#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#define rep(i,a,b) for(int i=a; i<=b; ++i)

int world_rank;
int world_size;

const unsigned char IGNOREME = 0;
int N;
int* VALUES;
int* ACCSUMS;
int TARGET_SUM;
enum Tag { CHUNK_TAG=1, SOLUTION_FOUND_TAG=2 };

// ------ UTILS --------------
int calc_bytes(int n_bits) {
    return (n_bits >> 3) + ((n_bits & 7) > 0);
}

// ------- BitSet ---------
typedef struct BitSet {
    unsigned char* bits;
    int size;
    int n_bytes;
} BitSet;
BitSet* BitSet_init(int size) {
    BitSet* bs = malloc(sizeof(BitSet));
    bs->size = size;
    bs->n_bytes = calc_bytes(size);
    bs->bits = calloc(bs->n_bytes, 1); // all bits set to 0 by default
}
void BitSet_free(BitSet* bs) {
    free(bs->bits);
    free(bs);
}
void BitSet_setbit(BitSet* bs, int i, bool val) {
    assert (i < bs->size);
    if (val) bs->bits[i >> 3] |= 1 << (i & 7);
    else bs->bits[i >> 3] &= ~(1 << (i & 7));
}
bool BitSet_getbit(BitSet* bs, int i) {
    assert (i < bs->size);
    return bs->bits[i >> 3] & (1 << (i & 7));
}
BitSet* BitSet_copy(BitSet* bs) {
    BitSet* copy = malloc(sizeof(BitSet));
    memcpy(copy, bs, sizeof(BitSet));
    copy->bits = malloc(copy->n_bytes);
    memcpy(copy->bits, bs->bits, copy->n_bytes);
    return copy;
}
int BitSet_countSetBits(BitSet* bs) {
    int count = 0;
    rep(i,0,bs->size-1) if (BitSet_getbit(bs,i)) count++;
    return count;
}

// ------- PartialSolution --------
typedef struct PartialSolution {
    BitSet* bitset;
    int index;
    int partial_sum;
} PartialSolution;
PartialSolution* PartialSolution_init() {
    PartialSolution* ps = malloc(sizeof(PartialSolution));
    ps->bitset = BitSet_init(N);
    ps->index = 0;
    ps->partial_sum = 0;
    return ps;
}
void PartialSolution_free(PartialSolution* ps) {
    BitSet_free(ps->bitset);
    free(ps);
}
PartialSolution* PartialSolution_copy(PartialSolution* ps) {
    PartialSolution* copy = malloc(sizeof(PartialSolution));
    memcpy(copy, ps, sizeof(PartialSolution));
    copy->bitset = BitSet_copy(ps->bitset);
    return copy;
}
bool PartialSolution_complete(PartialSolution* ps) {
    return ps->index == N;
}
void PartialSolution_setnextbit(PartialSolution* ps, bool val) {
    assert (ps->index < N);
    int i = ps->index++;
    BitSet_setbit(ps->bitset, i, val);
    if (val) ps->partial_sum += VALUES[i];
}

// ------- Queue -------------
typedef struct Queue {
    void** items;
    int front;
    int rear;
    int cap; // capacity
} Queue;
Queue* Queue_init(int cap) {
    assert (cap > 0);
    Queue* q = malloc(sizeof(Queue));
    q->items = malloc(cap * sizeof(void*));
    q->front = 0;
    q->rear = -1;
    q->cap = cap;
    return q;
}
void Queue_free(Queue* q) {
    free(q->items);
    free(q);
}
int Queue_size(Queue* q) {
    return q->rear - q->front + 1;
}
void Queue_push(Queue* q, void* item) {
    if (q->rear + 1 == q->cap) { // duplicate capacity if we run out of capacity
        printf("duplicating queue!");
        int new_cap = q->cap * 2;
        void** new_items = malloc(new_cap * sizeof(void*));
        int size = Queue_size(q);
        memcpy(new_items, &q->items[q->front], size * sizeof(void*));
        free(q->items);
        q->items = new_items;
        q->front = 0;
        q->rear = size - 1;
        q->cap = new_cap;
    }
    q->items[++(q->rear)] = item;
}
void* Queue_pop(Queue *q) {
    return q->items[(q->front)++];
}
void* Queue_item(Queue *q, int i) {
    return q->items[q->front + i];
}

// -------- BFS ------------
/* We explore the solution tree using BFS until we find enough
    nodes to distribute among parallel workers.
   In case we find a solution during BFS, we return it.
*/
PartialSolution* bfs(Queue** q_ptr) {
    assert (N > 0);
    assert (TARGET_SUM > 0);
    assert (world_size > 1);
    
    PartialSolution* full_sol = NULL;
    PartialSolution* ps;
    ps = PartialSolution_init();

    int min_n_jobs = (world_size-1) * 4;
    Queue* q = Queue_init(min_n_jobs * 4);

    Queue_push(q, ps);
    while (0 < Queue_size(q) && Queue_size(q) < min_n_jobs) {
        ps = Queue_pop(q);
        assert(PartialSolution_complete(ps) == false);
        
        // option 1: pick next item
        PartialSolution* ps1 = PartialSolution_copy(ps);
        PartialSolution_setnextbit(ps1, true);
        if (ps1->partial_sum <= TARGET_SUM) {
            if (ps1->partial_sum == TARGET_SUM) { // problem already solved!
                full_sol = ps1;
                goto answer;
            } else if (PartialSolution_complete(ps1)) {
                PartialSolution_free(ps1);
            } else {
                Queue_push(q, ps1);
            }
        } else {
            PartialSolution_free(ps1);
        }
        
        // option 2: skip next item
        PartialSolution* ps2 = PartialSolution_copy(ps);
        PartialSolution_setnextbit(ps2, false);
        assert (ps2->partial_sum == ps->partial_sum);
        assert (ps2->partial_sum < TARGET_SUM);
        if (PartialSolution_complete(ps2)) {
            PartialSolution_free(ps2);
        } else {
            Queue_push(q, ps2);
        }

        // free ps's allocated memory
        PartialSolution_free(ps);
    }
    
    answer:
    *q_ptr = q;
    return full_sol;
}

void print_solution(BitSet* sol) {
    if (sol == NULL) puts("0");
    else {
        printf("%d", BitSet_countSetBits(sol));
        rep(i,0,N-1) if (BitSet_getbit(sol, i)) printf(" %d", VALUES[i]);
        puts("");
    }
}

// -------- Backtracking (SERIAL) ------
bool search(BitSet* sol, int i, int partial_sum) {
    if (i == N) return partial_sum == TARGET_SUM;
    assert (i < N);    
    if (partial_sum + ACCSUMS[i] < TARGET_SUM) // pruning
        return false;
    // option 1: include i-th value
    int tmp = partial_sum + VALUES[i];
    if (tmp <= TARGET_SUM) {
        BitSet_setbit(sol, i, true);
        if (search(sol, i+1, tmp)) return true;
    }
    // option 2: skip i-th value
    BitSet_setbit(sol, i, false);
    if (search(sol, i+1, partial_sum)) return true;
    // default: nothing worked
    return false;
}

// -------- Backtracking (PARALLEL) ------
int search_flag;
MPI_Request search_request;
MPI_Status search_status;
bool search_interrupted;
unsigned int polling_count = 0;

bool search_parallel(BitSet* sol, int i, int partial_sum) {
    if (i == N) return partial_sum == TARGET_SUM;
    assert (i < N);
    if (partial_sum + ACCSUMS[i] < TARGET_SUM) // pruning
        return false;
    // check if someone else already solved the problem so we can stop
    if (!(++polling_count & 1023)) { // test every 1024 iterations
        MPI_Test(&search_request, &search_flag, &search_status);
        if (search_flag) {
            // printf("worker %d: getting interrupted!\n", world_rank);
            search_interrupted = true;
            return true;
        }
    }
    // option 1: include i-th value
    int tmp = partial_sum + VALUES[i];
    if (tmp <= TARGET_SUM) {
        BitSet_setbit(sol, i, true);
        if (search_parallel(sol, i+1, tmp)) return true;
    }
    // option 2: skip i-th value
    BitSet_setbit(sol, i, false);
    if (search_parallel(sol, i+1, partial_sum)) return true;
    // default: nothing worked
    return false;
}

// ACCSUMS[i] = the sum of all VALUES[j] for j in range i ... N-1
// This can be used later on for pruning during backtracking search
void precompute_ACCSUMS() {
    ACCSUMS = malloc((N+1) * sizeof(int));
    ACCSUMS[N] = 0;
    for (int i = N-1; i >= 0; i--) {
        ACCSUMS[i] = ACCSUMS[i+1] + VALUES[i];
    }
}

void read_input() {
    scanf("%d", &N);
    VALUES = malloc(N * sizeof(int));
    rep(i,0,N-1) scanf("%d", VALUES+i);
    scanf("%d", &TARGET_SUM);
    precompute_ACCSUMS();
}

int main(int argc, char** argv) {
    // initialize the MPI environment
    MPI_Init(&argc, &argv);

    // find out world rank and size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // printf("world_rank=%d, world_size=%d\n", world_rank, world_size);

    // ==============================
    // special case: a single process
    if (world_size == 1) {
        assert (world_rank == 0);
        read_input(); // read input from stdin
        // execute a serial backtracking
        BitSet* sol = BitSet_init(N);
        if (search(sol, 0, 0)) {
            print_solution(sol);
        } else {
            puts("0");
        }
        MPI_Finalize();
        return 0;
    }

    // ==============================
    // general case: 2+ processes
    bool done;
    if (world_rank == 0) {
        read_input(); // read input from stdin
        // execute initial BFS        
        Queue* q;
        PartialSolution* full_sol;
        if (N > 0 && TARGET_SUM > 0) {
            full_sol = bfs(&q);
            if (full_sol != NULL) { // we found a solution during BFS
                print_solution(full_sol->bitset);
                done = true;
            } else if (Queue_size(q) == 0) { // there are no solutions
                puts("0");
                done = true;
            } else { // there are nodes to explore -> not done yet
                done = false;
            }
        } else { // trivially 0
            puts("0");
            done = true;
        }
        // broadcast 'done' status to workers
        MPI_Bcast(&done,1,MPI_C_BOOL,0,MPI_COMM_WORLD);
        // if not done, then distribute jobs to workers
        if (!done) {
            // 1) Broadcast inputs to workers
            MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
            MPI_Bcast(VALUES,N,MPI_INT,0,MPI_COMM_WORLD);
            MPI_Bcast(&TARGET_SUM,1,MPI_INT,0,MPI_COMM_WORLD);
            // 2) Distribute nodes collected in BFS queue among workers
            int n_workers = world_size-1;
            int n_jobs = Queue_size(q);
            int ref_chunk_size = n_jobs / n_workers;
            int remainder = n_jobs % n_workers;
            int offset = 0;
            rep(rank,1,world_size-1) { // send chunk of work from queue to each worker                
                int chunk_size = ref_chunk_size + (rank <= remainder);
                MPI_Send(&chunk_size, 1, MPI_INT, rank, CHUNK_TAG, MPI_COMM_WORLD);
                rep(i,0,chunk_size-1) {
                    PartialSolution* ps = Queue_item(q, offset+i);
                    MPI_Send(&ps->index, 1, MPI_INT, rank, CHUNK_TAG, MPI_COMM_WORLD);
                    MPI_Send(&ps->partial_sum, 1, MPI_INT, rank, CHUNK_TAG, MPI_COMM_WORLD);
                    MPI_Send(ps->bitset->bits, ps->bitset->n_bytes,
                            MPI_UNSIGNED_CHAR, rank, CHUNK_TAG, MPI_COMM_WORLD);
                }
                offset += chunk_size;
            }
            assert (offset == n_jobs);            
            // 3) Blockingly wait for all workers to finish in a first-come first-served basis.
            // When the first solution is received, then we NON-blockingly broadcast that
            // we are done, so other workers can stop working immediately.
            BitSet* tmp = BitSet_init(N);
            BitSet* first_sol = NULL;
            int count = n_workers;
            while(count-- > 0) {
                MPI_Status status;
                int length;
                MPI_Probe(MPI_ANY_SOURCE, SOLUTION_FOUND_TAG, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &length);
                MPI_Recv(tmp->bits, length, MPI_UNSIGNED_CHAR,
                        status.MPI_SOURCE, SOLUTION_FOUND_TAG,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (length > 0 && first_sol == NULL) {
                    first_sol = BitSet_copy(tmp);
                    done = true;
                    // printf("  ***** ROOT::solution found by worker %d! MPI_Ibcast-ing!\n", status.MPI_SOURCE);
                    MPI_Ibcast(&done, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD, &search_request);
                }
            }
            print_solution(first_sol);
        }
    } else { // non-ROOT
        // wait for 'done'
        MPI_Bcast(&done,1,MPI_C_BOOL,0,MPI_COMM_WORLD);
        // if not done, then receive jobs from root
        if (!done) {
            // 1) receive inputs
            MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
            VALUES = malloc(N * sizeof(int));
            MPI_Bcast(VALUES,N,MPI_INT,0,MPI_COMM_WORLD);
            MPI_Bcast(&TARGET_SUM,1,MPI_INT,0,MPI_COMM_WORLD);
            precompute_ACCSUMS();
            // 2) receive chunk of work
            int chunk_size;
            MPI_Status status;
            MPI_Recv(&chunk_size, 1, MPI_INT, 0, CHUNK_TAG, MPI_COMM_WORLD, &status);
            BitSet** sols = malloc(chunk_size * sizeof(BitSet));
            int* indexes = malloc(chunk_size * sizeof(int));
            int* partial_sums = malloc(chunk_size * sizeof(int));
            rep(i,0,chunk_size-1) {
                sols[i] = BitSet_init(N);
                MPI_Recv(&indexes[i], 1, MPI_INT, 0, CHUNK_TAG, MPI_COMM_WORLD, &status);
                MPI_Recv(&partial_sums[i], 1, MPI_INT, 0, CHUNK_TAG, MPI_COMM_WORLD, &status);
                MPI_Recv(sols[i]->bits, sols[i]->n_bytes,
                        MPI_UNSIGNED_CHAR, 0, CHUNK_TAG, MPI_COMM_WORLD, &status);
            }
            // 3) search for solutions in the solution subtrees received
            search_interrupted = false;
            MPI_Ibcast(&done, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD, &search_request);
            BitSet* found_sol = NULL;
            rep(i,0,chunk_size-1) {
                if (search_parallel(sols[i], indexes[i], partial_sums[i])) {
                    found_sol = sols[i]; break;
                }
            }
            // printf("worker %d: search_interrupted = %d\n", world_rank, search_interrupted);
            // 4) send solution found to root
            // printf("worker %d: BEFORE sending\n", world_rank);
            if (found_sol == NULL || search_interrupted) {                
                MPI_Send(&IGNOREME, 0,
                        MPI_UNSIGNED_CHAR, 0, SOLUTION_FOUND_TAG, MPI_COMM_WORLD);
            } else {
                MPI_Send(found_sol->bits, found_sol->n_bytes,
                        MPI_UNSIGNED_CHAR, 0, SOLUTION_FOUND_TAG, MPI_COMM_WORLD);
            }
            // printf("worker %d: AFTER sending\n", world_rank);
        }
    }
    
    // printf("worker %d: BEFORE MPI_Finalize()\n", world_rank);
    MPI_Finalize();
    // printf("worker %d: AFTER MPI_Finalize()\n", world_rank);
    return 0;
}
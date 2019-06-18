To compile, run
> make

To execute, run
> mpiexec -n {num_process} ./subsetsum-mpi.o < {input_file}

for example:
> mpiexec -n 20 ./subsetsum-mpi.o < ./testcases/N\=30_nosolution.txt
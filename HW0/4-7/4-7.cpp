#include <iostream>
#include <mpi.h>

int main(int argc, char** argv)
{
    int rank, size, sum, global_sum;
    double elapsed_time;
    MPI_Init(&argc, &argv);
    MPI_Barrier (MPI_COMM_WORLD);
    elapsed_time = - MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    sum = rank + 1;
    MPI_Reduce (&sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(!rank){
        std::cout << "Calulated Global sum: " << global_sum << std::endl;
        std::cout << "Correct Global sum: " << size * (size + 1) / 2 << std::endl;
    }

    elapsed_time += MPI_Wtime();
    std::cout << "elapsed_time: " << elapsed_time << std::endl;
    MPI_Finalize();
    return 0;
}

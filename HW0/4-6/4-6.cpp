#include <iostream>
#include <mpi.h>

int main(int argc, char** argv)
{
    int rank;
    double elapsed_time;
    MPI_Init(&argc, &argv);
    MPI_Barrier (MPI_COMM_WORLD);
    elapsed_time = - MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout << "hello, world, from process " << rank << std::endl;
    elapsed_time += MPI_Wtime();
    std::cout << "elapsed_time: " << elapsed_time << std::endl;
    MPI_Finalize();
    return 0;
}

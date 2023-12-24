#include <iostream>
#include <mpi.h>

int main(int argc, char** argv)
{
    int rank, size;
    double elapsed_time;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Barrier (MPI_COMM_WORLD);
    int* data = new int[5242880]; // 20M bytes
    elapsed_time = - MPI_Wtime();

    if(rank == 1){
        MPI_Recv (data, 5242880, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send (data, 5242880, MPI_INT, 0, 2, MPI_COMM_WORLD);
    }

    if(rank == 0){
        MPI_Send (data, 5242880, MPI_INT, 1, 1, MPI_COMM_WORLD);
        MPI_Recv (data, 5242880, MPI_INT, 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    elapsed_time += MPI_Wtime();
    if(rank == 0){
        std::cout << "elapsed_time: " << elapsed_time << std::endl;
    }

    delete [] data;
    MPI_Finalize();

    return 0;
}

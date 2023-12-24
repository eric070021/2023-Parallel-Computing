#include <iostream>
#include <mpi.h>

#define n 50

double f(int i){
    double x;
    x = static_cast<double>(i) / static_cast<double>(n);
    return 4.0 / (1.0 + x * x);
}

int main(int argc, char** argv)
{
    int id, p;
    double elapsed_time;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier (MPI_COMM_WORLD);
    elapsed_time = - MPI_Wtime();

    double area = 0.0, global_area;
    for (int i = id + 1; i <= n/2; i += p)
        area += (4.0*f(2*i - 1) + 2*f(2*i));
    
    MPI_Reduce (&area, &global_area, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(id == 0){
        global_area += (f(0) - f(n));
        global_area /= (3.0 * n);
    }
   
    elapsed_time += MPI_Wtime();
    if(id == 0){
        std::cout.precision(15);
        std::cout << "Approximatation of pi: " << global_area << std::endl;
        std::cout << "elapsed_time: " << elapsed_time << std::endl;
    }

    MPI_Finalize();

    return 0;
}

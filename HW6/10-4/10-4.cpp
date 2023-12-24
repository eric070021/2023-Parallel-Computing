#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <mpi.h>

using namespace std;

#define s 2
#define d 0.3
#define num_iter 10000000

bool check_distance(double x, double y, double z){
    double norm2 = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    double distance = (x + y + z) / (sqrt(3) * norm2);
    distance = sin(acos(distance)) * norm2;
    if(distance > (d / 2)) return true;
    else return false;
}

int main(int argc, char** argv){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(-s/2, s/2);
    double x, y, z;
    int point_cnt = 0;
    int global_point_cnt = 0;
    int id, p;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /*************************Start Timer***************************/
    MPI_Barrier (MPI_COMM_WORLD);
    double elapsed_time;
    elapsed_time = - MPI_Wtime();

    for(int i = 0; i < (num_iter / p); i++){
        x = distribution(gen);
        y = distribution(gen);
        z = distribution(gen);

        if(check_distance(x, y, z)) point_cnt++;
    }

    MPI_Reduce(&point_cnt, &global_point_cnt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(!id) {
        cout << "area: " << (static_cast<double>(global_point_cnt) / ((num_iter / p) * p)) * pow(s, 3) << endl;
    }

    /*************************Stop Timer***************************/
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    if(!id){
        cout << "Processes = " << p << ", Time = " << elapsed_time << std::endl;
    }

    MPI_Finalize();

    return 0;
}
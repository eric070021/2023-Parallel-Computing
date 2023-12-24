#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <mpi.h>

using namespace std;

#define PROC_READY 0
#define PERFECT_NUM 1
#define TERMINATE 2
#define TEST_NUM 3

typedef unsigned long long int uint64_t;

bool is_prime(uint64_t n) {
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0) return false; // Exclude even numbers except 2

    // Check divisibility for odd numbers starting from 3 up to sqrt(n)
    for (int i = 3; i <= sqrt(n); i += 2) {
        if (n % i == 0) {
            return false; // If divisible, it's not a prime number
        }
    }
    return true;
}

int main(int argc, char** argv){
    int id, p;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    if(p < 2){
        cerr << "Need at least 2 processors\n";
        exit(1);
    }
    
    /*************************Start Timer***************************/
    MPI_Barrier (MPI_COMM_WORLD);
    double elapsed_time;
    elapsed_time = - MPI_Wtime();

    if(!id){ // manager
        uint64_t send_buf;          // buffer used to send data to worker
        uint64_t recv_buf;          // buffer used to receive data from worker
        uint64_t worker;            // worker id
        uint64_t tag;               // message tag
        uint64_t number = 1;        // number used for prime test
        MPI_Status status;          // MPI status
        vector<uint64_t> primeVec;  // Vector used to store perfect number

        while(primeVec.size() < 6){
            MPI_Recv(&recv_buf, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, 
                MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            worker = status.MPI_SOURCE;
            tag = status.MPI_TAG;
            if(tag == PROC_READY){
                send_buf = number;
                MPI_Send(&send_buf, 1, MPI_UNSIGNED_LONG_LONG, worker, TEST_NUM, MPI_COMM_WORLD);
                number++;
            }
            else if(tag == PERFECT_NUM){
                primeVec.push_back(recv_buf);
            }
        }

        for(int i = 1; i < p; i++){
            MPI_Send(&send_buf, 1, MPI_UNSIGNED_LONG_LONG, i, TERMINATE, MPI_COMM_WORLD);
        }

        sort(primeVec.begin(), primeVec.end());
        cout << "Perfect Number: ";
        for(int i = 0; i < primeVec.size(); i++){
            cout << primeVec[i] << " ";
        }
    }
    else{ // worker
        uint64_t send_buf;      // buffer used to send data to manager
        uint64_t recv_buf;      // buffer used to receive data from manager
        uint64_t tag;           // message tag
        bool prime;             // number is prime or not
        MPI_Status status;      // MPI status

        while (true) {
            MPI_Send(&send_buf, 1, MPI_UNSIGNED_LONG_LONG, 0, PROC_READY, MPI_COMM_WORLD);
            MPI_Recv(&recv_buf, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            tag = status.MPI_TAG;
            if(tag == TERMINATE){
                break;
            }
            else if(tag == TEST_NUM){
                if(is_prime(static_cast<uint64_t>(pow(2, recv_buf) - 1))){
                    send_buf = static_cast<uint64_t>((pow(2, recv_buf) - 1) * pow(2, recv_buf - 1));
                    MPI_Send(&send_buf, 1, MPI_UNSIGNED_LONG_LONG, 0, PERFECT_NUM, MPI_COMM_WORLD);
                }
            }
        }
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

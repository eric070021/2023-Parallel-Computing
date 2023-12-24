#include <iostream>
#include <vector>
#include <math.h>
#include <mpi.h>

#define MIN(a,b)   ((a)<(b)?(a):(b))
#define MAX(a,b)   ((a)>(b)?(a):(b))

#define n 1000000
#define BLOCK_SIZE 16384

int seq_sieve (bool *small_primes, int small_prime_array_size, int sqrt_n)
{
   int i, j;
   int prime_index;
   int prime_value;
   int count;

   /* small_primes[i] represents integer 2i+3 */

   for (i = 0; i < small_prime_array_size; i++) small_primes[i] = true;
   prime_index = 0;
   prime_value = 3;
   while (prime_value * prime_value <= sqrt_n) {
      j = prime_value * prime_value / 2 - 1;
      while (j < small_prime_array_size) {
         small_primes[j] = false;
         j += prime_value;
      }
      while (small_primes[++prime_index] == false);
      prime_value = 2*prime_index + 3;
   }
   count = 0;
   for (i = 0; i < small_prime_array_size; i++)
        if (small_primes[i] == true)
            count++;

   return count;
}

int main(int argc, char** argv)
{
    int i, j, k;
    int global_count;           /* Total number of primes up through n */
    int id;                     /* Process ID number */
    int p;                      /* Number of processes */
    double elapsed_time;        /* Elapsed wall clock time */
    int sqrt_n;                 /* Square root of n, rounded down */
    bool *small_primes;         /* Used to sieve odd primes up to sqrt(n) */
    int small_prime_count;      /* Number of odd primes through sqrt(n) */
    int *small_prime_values;    /* List of odd primes up to sqrt(n) */
    int size;                   /* Number of elements in 'primes' */
    bool *primes;               /* Process's portion of integers 3,5,...,n */
    int blocks_num;             /* Number of blocks in subarray */
    int low_proc_value;         /* Integer represented by 1st el in 'primes' */
    int high_proc_value;        /* Integer represented by last el in 'primes' */

    MPI_Init(&argc, &argv);
    MPI_Barrier (MPI_COMM_WORLD);

    /* Start timer */

    elapsed_time = - MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);
    sqrt_n = static_cast<int>(sqrt(static_cast<double>(n)));
    int small_prime_array_size = (sqrt_n - 1)/2;
    small_primes = new bool [small_prime_array_size];
    small_prime_count = seq_sieve (small_primes, small_prime_array_size, sqrt_n);
    small_prime_values = new int [small_prime_count];
    for (i = 0, j = 0; i < small_prime_array_size; i++)
        if (small_primes[i])
            small_prime_values[j++] = 2*i+3;
    int els = (n-1) / 2;
    int smaller_size = els / p;
    int larger_size = smaller_size + 1;
    int num_larger_blocks = els % p;
    if (id < num_larger_blocks) size = larger_size;
    else size = smaller_size;
    low_proc_value = 2*(id*smaller_size + MIN(id, num_larger_blocks)) + 3;
    high_proc_value = low_proc_value + 2 * (size-1);
    primes = new bool [size];

    blocks_num = static_cast<int>(ceil(static_cast<double>(size) / BLOCK_SIZE));
    for (i = 0; i < blocks_num; i++) {
        int low_block_index = i * BLOCK_SIZE;
        int high_block_index = ((i + 1) * BLOCK_SIZE) - 1;
        if (high_block_index > size-1) high_block_index = size-1;
        for (j = low_block_index; j <= high_block_index; j++) primes[j] = true;
        int low_block_value = low_proc_value + 2*low_block_index;
        int high_block_value = low_proc_value + 2*high_block_index;
        for (j = 0; j < small_prime_count; j++) {
            int index;
            int current_prime = small_prime_values[j];
            if (current_prime * current_prime > low_block_value)
                index = low_block_index + (current_prime * current_prime - low_block_value)/2;
            else {
                int r = low_block_value % current_prime;
                if (!r) index = low_block_index;
                else if ((current_prime - r) & 1)
                index = low_block_index + (2*current_prime - r)/2;
                else index = low_block_index + (current_prime - r)/2;
            }
            if (index > high_block_index) break;
            for (k = index; k <= high_block_index; k+= current_prime)
                primes[k] = false;
        }
    }
    /* Count the max gap of two consecutive prime numbers */
    int gap = -1, global_gap = -1;

    std::vector<int> prime_vec;
    for (i = 0; i < size; i++)
        if (primes[i])
            prime_vec.push_back(2*i+low_proc_value);

    for (i = 0; i < (prime_vec.size() - 1); i++) {
        int mid_gap = prime_vec[i+1] - prime_vec[i];
        gap = MAX(mid_gap, gap);
    }

    int* f_prime_arr = new int [p];
    int* l_prime_arr = new int [p];
    MPI_Gather(&prime_vec.front(), 1, MPI_INT, f_prime_arr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&prime_vec.back(), 1, MPI_INT, l_prime_arr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!id) {
        for (i = 0; i < (p - 1); i++) {
            int mid_gap =  f_prime_arr[i+1] - l_prime_arr[i];
            gap = MAX(mid_gap, gap);
        }
    }

    MPI_Reduce (&gap, &global_gap, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    if (!id) {
        std::cout << "Max gap of two consecutive prime numbers is " << global_gap << std::endl;
        std::cout << "elapsed_time: " << elapsed_time << std::endl;
    }

    /* Delete heap elements */
    delete [] small_primes;
    delete [] small_prime_values;
    delete [] primes;
    delete [] f_prime_arr;
    delete [] l_prime_arr;
    MPI_Finalize();

    return 0;
}
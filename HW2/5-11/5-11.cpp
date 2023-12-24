#include <iostream>
#include <mpi.h>
#include <mpfr.h>
#include <string.h>

#define n 1000000
#define d 100
#define max_precision (d+10)
#define max_string_length (max_precision+2)

int main(int argc, char** argv)
{
    int id, p;
    double elapsed_time;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier (MPI_COMM_WORLD);
    elapsed_time = - MPI_Wtime();

    mpfr_t sum;
    mpfr_init2 (sum, max_precision);
    mpfr_set_d(sum, 0.0, MPFR_RNDN);

    /* Calculate partial sum */
    mpfr_t numerator, denominator;
    mpfr_init2 (numerator, max_precision);
    mpfr_init2 (denominator, max_precision);
    mpfr_set_d(numerator, 1.0, MPFR_RNDN);
    for (int i = id + 1; i <= n; i += p) {
        mpfr_set_d(denominator, i, MPFR_RNDN);
        mpfr_div (denominator, numerator, denominator, MPFR_RNDN);
        mpfr_add (sum, sum, denominator, MPFR_RNDN);
    }
    mpfr_clear(numerator);
    mpfr_clear(denominator);

    /* Serialize partial sum to char array */
    long exponent;
    char* sum_get_str;
    sum_get_str = mpfr_get_str (NULL, &exponent, 10, max_precision, sum, MPFR_RNDN);
    int currentLength = strlen(sum_get_str);
    char *sum_send = new char [max_string_length];
    for (int i = 0; i < exponent; i++) {
        sum_send[i] = sum_get_str[i];
    }
    sum_send[exponent] = '.';
    for (int i = exponent; i < currentLength; i++) {
        sum_send[i+1] = sum_get_str[i];
    }
    sum_send[currentLength + 1] = '\0'; // null character

    char* sum_recv = NULL;
    if (id == 0) {
        // Allocate memory for the receive buffer
        sum_recv = new char[max_string_length * p];
    }

    // Gather the strings into the receive buffer on id 0
    MPI_Gather(sum_send, max_string_length, MPI_CHAR, sum_recv, max_string_length, MPI_CHAR, 0, MPI_COMM_WORLD);

    mpfr_t sum_global;
    mpfr_init2 (sum_global, max_precision);
    mpfr_set_d(sum_global, 0.0, MPFR_RNDN);
    if (id == 0) {
        mpfr_t sum_partial;
        mpfr_init2 (sum_partial, max_precision);
        for(int proc = 0; proc < p; proc++) {
            mpfr_set_str(sum_partial, &sum_recv[proc * max_string_length], 10, MPFR_RNDN);
            mpfr_add (sum_global, sum_global, sum_partial, MPFR_RNDN);
        }

        /* Serialize global sum to char array */
        sum_get_str = mpfr_get_str (NULL, &exponent, 10, max_precision, sum_global, MPFR_RNDN);
        for (int i = 0; i < exponent; i++) {
            sum_send[i] = sum_get_str[i];
        }
        sum_send[exponent] = '.';
        for (int i = exponent; i < (d + exponent); i++) {
            sum_send[i+1] = sum_get_str[i];
        }
        sum_send[d + exponent + 1] = '\0'; // null character

        delete [] sum_recv;
    }
   
    elapsed_time += MPI_Wtime();

    if(id == 0){
        std::cout << "S1,000,000: " << sum_send <<std::endl;
        std::cout << "elapsed_time: " << elapsed_time << std::endl;
    }

    /* Delete heap elements */
    delete [] sum_send;
    mpfr_free_str(sum_get_str);

    MPI_Finalize();

    return 0;
}

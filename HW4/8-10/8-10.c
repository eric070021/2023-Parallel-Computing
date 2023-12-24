#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

/* Change these two definitions when the matrix and vector
   element types change */

typedef double dtype;
#define mpitype MPI_DOUBLE

/************************* MACROS **************************/

#define DATA_MSG           0
#define PROMPT_MSG         1
#define RESPONSE_MSG       2

#define OPEN_FILE_ERROR    -1
#define MALLOC_ERROR       -2
#define TYPE_ERROR         -3

#define MIN(a,b)           ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) \
                     (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(j,p,n) (((p)*((j)+1)-1)/(n))
#define PTR_SIZE           (sizeof(void*))
#define CEILING(i,j)       (((i)+(j)-1)/(j))

/***************** MISCELLANEOUS FUNCTIONS *****************/

/*
 *   Function 'my_malloc' is called when a process wants
 *   to allocate some space from the heap. If the memory
 *   allocation fails, the process prints an error message
 *   and then aborts execution of the program.
 */

void *my_malloc (
   int id,     /* IN - Process rank */
   int bytes)  /* IN - Bytes to allocate */
{
   void *buffer;
   if ((buffer = malloc ((size_t) bytes)) == NULL) {
      printf ("Error: Malloc failed for process %d\n", id);
      fflush (stdout);
      MPI_Abort (MPI_COMM_WORLD, MALLOC_ERROR);
   }
   return buffer;
}

/*
 *   Given MPI_Datatype 't', function 'get_size' returns the
 *   size of a single datum of that data type.
 */

int get_size (MPI_Datatype t) {
   if (t == MPI_BYTE) return sizeof(char);
   if (t == MPI_DOUBLE) return sizeof(double);
   if (t == MPI_FLOAT) return sizeof(float);
   if (t == MPI_INT) return sizeof(int);
   printf ("Error: Unrecognized argument to 'get_size'\n");
   fflush (stdout);
   MPI_Abort (MPI_COMM_WORLD, TYPE_ERROR);
}

/*
 *   Print elements of a singly-subscripted array.
 */

void print_subvector (
   void        *a,       /* IN - Array pointer */
   MPI_Datatype dtype,   /* IN - Array type */
   int          n)       /* IN - Array size */
{
   int i;

   for (i = 0; i < n; i++) {
      if (dtype == MPI_DOUBLE)
         printf ("%6.3f ", ((double *)a)[i]);
      else {
         if (dtype == MPI_FLOAT)
            printf ("%6.3f ", ((float *)a)[i]);
         else if (dtype == MPI_INT)
            printf ("%6d ", ((int *)a)[i]);
      }
   }
}

/*
 *   Print a vector that is block distributed among the
 *   processes in a communicator.
 */

void print_block_vector (
   void        *v,       /* IN - Address of vector */
   MPI_Datatype dtype,   /* IN - Vector element type */
   int          n,       /* IN - Elements in vector */
   MPI_Comm     comm)    /* IN - Communicator */
{
   int        datum_size; /* Bytes per vector element */
   int        i;
   int        prompt;     /* Dummy variable */
   MPI_Status status;     /* Result of receive */
   void       *tmp;       /* Other process's subvector */
   int        id;         /* Process rank */
   int        p;          /* Number of processes */

   MPI_Comm_size (comm, &p);
   MPI_Comm_rank (comm, &id);
   datum_size = get_size (dtype);

   if (!id) {
      print_subvector (v, dtype, BLOCK_SIZE(id,p,n));
      if (p > 1) {
         tmp = my_malloc (id,BLOCK_SIZE(p-1,p,n)*datum_size);
         for (i = 1; i < p; i++) {
            MPI_Send (&prompt, 1, MPI_INT, i, PROMPT_MSG,
               comm);
            MPI_Recv (tmp, BLOCK_SIZE(i,p,n), dtype, i,
               RESPONSE_MSG, comm, &status);
            print_subvector (tmp, dtype, BLOCK_SIZE(i,p,n));
         }
         free (tmp);
      }
      printf ("\n\n");
   } else {
      MPI_Recv (&prompt, 1, MPI_INT, 0, PROMPT_MSG, comm,
         &status);
      MPI_Send (v, BLOCK_SIZE(id,p,n), dtype, 0,
         RESPONSE_MSG, comm);
   }
}

/*
 *   Function 'read_checkerboard_matrix' reads a matrix from
 *   a file. The first two elements of the file are integers
 *   whose values are the dimensions of the matrix ('m' rows
 *   and 'n' columns). What follows are 'm'*'n' values
 *   representing the matrix elements stored in row-major
 *   order.  This function allocates blocks of the matrix to
 *   the MPI processes.
 *
 *   The number of processes must be a square number.
 */

void read_checkerboard_matrix (
   char *s,              /* IN - File name */
   void ***subs,         /* OUT - 2D array */
   void **storage,       /* OUT - Array elements */
   MPI_Datatype dtype,   /* IN - Element type */
   int *m,               /* OUT - Array rows */
   int *n,               /* OUT - Array cols */
   MPI_Comm grid_comm)   /* IN - Communicator */
{
   void      *buffer;         /* File buffer */
   int        coords[2];      /* Coords of proc receiving
                                 next row of matrix */
   int        datum_size;     /* Bytes per elements */
   int        dest_id;        /* Rank of receiving proc */
   int        grid_coord[2];  /* Process coords */
   int        grid_id;        /* Process rank */
   int        grid_period[2]; /* Wraparound */
   int        grid_size[2];   /* Dimensions of grid */
   int        i, j, k;
   FILE      *infileptr;      /* Input file pointer */
   void      *laddr;          /* Used when proc 0 gets row */
   int        local_cols;     /* Matrix cols on this proc */
   int        local_rows;     /* Matrix rows on this proc */
   void     **lptr;           /* Pointer into 'subs' */
   int        p;              /* Number of processes */
   void      *raddr;          /* Address of first element
                                 to send */
   void      *rptr;           /* Pointer into 'storage' */
   MPI_Status status;         /* Results of read */

   MPI_Comm_rank (grid_comm, &grid_id);
   MPI_Comm_size (grid_comm, &p);
   datum_size = get_size (dtype);

   /* Process 0 opens file, gets number of rows and
      number of cols, and broadcasts this information
      to the other processes. */

   if (grid_id == 0) {
      infileptr = fopen (s, "r");
      if (infileptr == NULL) *m = 0;
      else {
         fread (m, sizeof(int), 1, infileptr);
         fread (n, sizeof(int), 1, infileptr);
      }
   }
   MPI_Bcast (m, 1, MPI_INT, 0, grid_comm);

   if (!(*m)) MPI_Abort (MPI_COMM_WORLD, OPEN_FILE_ERROR);

   MPI_Bcast (n, 1, MPI_INT, 0, grid_comm);

   /* Each process determines the size of the submatrix
      it is responsible for. */

   MPI_Cart_get (grid_comm, 2, grid_size, grid_period,
      grid_coord);
   local_rows = BLOCK_SIZE(grid_coord[0],grid_size[0],*m);
   local_cols = BLOCK_SIZE(grid_coord[1],grid_size[1],*n);

   /* Dynamically allocate two-dimensional matrix 'subs' */

   *storage = my_malloc (grid_id,
      local_rows * local_cols * datum_size);
   *subs = (void **) my_malloc (grid_id,local_rows*PTR_SIZE);
   lptr = (void *) *subs;
   rptr = (void *) *storage;
   for (i = 0; i < local_rows; i++) {
      *(lptr++) = (void *) rptr;
      rptr += local_cols * datum_size;
   }

   /* Grid process 0 reads in the matrix one row at a time
      and distributes each row among the MPI processes. */

   if (grid_id == 0)
      buffer = my_malloc (grid_id, *n * datum_size);

   /* For each row of processes in the process grid... */
   for (i = 0; i < grid_size[0]; i++) {
      coords[0] = i;

      /* For each matrix row controlled by this proc row...*/
      for (j = 0; j < BLOCK_SIZE(i,grid_size[0],*m); j++) {

         /* Read in a row of the matrix */

         if (grid_id == 0) {
            fread (buffer, datum_size, *n, infileptr);
         }

         /* Distribute it among process in the grid row */

         for (k = 0; k < grid_size[1]; k++) {
            coords[1] = k;

            /* Find address of first element to send */
            raddr = buffer +
               BLOCK_LOW(k,grid_size[1],*n) * datum_size;

            /* Determine the grid ID of the process getting
               the subrow */
            MPI_Cart_rank (grid_comm, coords, &dest_id);

            /* Process 0 is responsible for sending...*/
            if (grid_id == 0) {

               /* It is sending (copying) to itself */
               if (dest_id == 0) {
                  laddr = (*subs)[j];
                  memcpy (laddr, raddr, local_cols * datum_size);

               /* It is sending to another process */
               } else {
                  MPI_Send (raddr,
                     BLOCK_SIZE(k,grid_size[1],*n), dtype,
                  dest_id, 0, grid_comm);
               }

            /* Process 'dest_id' is responsible for
               receiving... */
            } else if (grid_id == dest_id) {
               MPI_Recv ((*subs)[j], local_cols, dtype, 0,
                  0, grid_comm,&status);
            }
         }
      }
   }
   if (grid_id == 0) free (buffer);
}

/*
 *   Open a file containing a vector, read its contents,
 *   and distributed the elements by block among the
 *   processes in a communicator.
 */

void read_block_vector (
    char        *s,      /* IN - File name */
    void       **v,      /* OUT - Subvector */
    MPI_Datatype dtype,  /* IN - Element type */
    int         *n,      /* OUT - Vector length */
    MPI_Comm     comm)   /* IN - Communicator */
{
   int        datum_size;   /* Bytes per element */
   int        i;
   FILE      *infileptr;    /* Input file pointer */
   int        local_els;    /* Elements on this proc */
   MPI_Status status;       /* Result of receive */
   int        id;           /* Process rank */
   int        p;            /* Number of processes */
   int        x;            /* Result of read */

   datum_size = get_size (dtype);
   MPI_Comm_size(comm, &p);
   MPI_Comm_rank(comm, &id);

   /* Process p-1 opens file, determines number of vector
      elements, and broadcasts this value to the other
      processes. */

   if (id == (p-1)) {
      infileptr = fopen (s, "r");
      if (infileptr == NULL) *n = 0;
      else fread (n, sizeof(int), 1, infileptr);
   }
   MPI_Bcast (n, 1, MPI_INT, p-1, comm);
   if (! *n) {
      if (!id) {
         printf ("Input file '%s' cannot be opened\n", s);
         fflush (stdout);
      }
   }

   /* Block mapping of vector elements to processes */

   local_els = BLOCK_SIZE(id,p,*n);

   /* Dynamically allocate vector. */

   *v = my_malloc (id, local_els * datum_size);
   if (id == (p-1)) {
      for (i = 0; i < p-1; i++) {
         x = fread (*v, datum_size, BLOCK_SIZE(i,p,*n),
            infileptr);
         MPI_Send (*v, BLOCK_SIZE(i,p,*n), dtype, i, DATA_MSG,
            comm);
      }
      x = fread (*v, datum_size, BLOCK_SIZE(id,p,*n),
             infileptr);
      fclose (infileptr);
   } else {
      MPI_Recv (*v, BLOCK_SIZE(id,p,*n), dtype, p-1, DATA_MSG,
         comm, &status);
   }
}

/***************** MAIN FUNCTIONS *****************/

int main (int argc, char *argv[]) {
   dtype **a;       /* First factor, a matrix */
   dtype *b;        /* Second factor, a vector */
   int    base;
   dtype *c_block;  /* Partial product vector */
   dtype *c_sums;
   int       cols;
   double    max_seconds;
   double    elapsed_time;    /* Elapsed time for matrix-vector multiply */
   dtype *storage;  /* Matrix elements stored here */
   int       grid_id;
   int    grid_size[2]; /* Number of procs in each grid dimension */
   MPI_Comm grid_comm;
   MPI_Comm row_comm;
   MPI_Comm col_comm;
   dtype *tmp;
   int    i, j;     /* Loop indices */
   int    id;       /* Process ID number */
   int    m;        /* Rows in matrix */
   int    n;        /* Columns in matrix */
   int    nprime;   /* Elements in vector */
   int    p;        /* Number of processes */
   int    sqrt_p;   /* Square Root of Number of processes */
   int    rows;     /* Number of rows on this process */
   int   *recv_cnt;
   int   *recv_disp;
   int   *send_cnt;
   int   *send_disp;
   int    grid_coords[2];
   int    periodic[2];
   int    is_square;
   dtype *btrans;
   int    send_length;
   int    b_length;
   int    src;
   int    dest;
   int    coords[2];
   MPI_Status status;

   MPI_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);

   /* Create 2D communicators */
   grid_size[0] = grid_size[1] = 0;
   MPI_Dims_create (p, 2, grid_size);
   periodic[0] = periodic[1] = 0;
   MPI_Cart_create (MPI_COMM_WORLD, 2, grid_size, periodic, 1, &grid_comm);
   MPI_Comm_rank (grid_comm, &grid_id);
   MPI_Cart_coords (grid_comm, grid_id, 2, grid_coords);
   MPI_Comm_split (grid_comm, grid_coords[0], grid_coords[1], &row_comm);
   MPI_Comm_split (grid_comm, grid_coords[1], grid_coords[0], &col_comm);
  
   /* Read Matrix and Vector from file */
   read_checkerboard_matrix (argv[1], (void *) &a,
      (void *) &storage, mpitype, &m, &n, grid_comm);
   rows = BLOCK_SIZE(grid_coords[0],grid_size[0],m);
   cols = BLOCK_SIZE(grid_coords[1],grid_size[1],n);
   if (grid_coords[1] == 0) {
      read_block_vector (argv[2], (void *) &b, mpitype, &nprime, col_comm);
   }

   /* Malloc needed buffers */
   c_block = (dtype *) malloc (rows * sizeof(dtype));
   c_sums = (dtype *) malloc (rows * sizeof(dtype));
   b_length = BLOCK_SIZE(grid_coords[1], grid_size[1], n);
   btrans = (dtype *) malloc (b_length * sizeof(dtype));
   tmp = (dtype *) malloc (n * sizeof(dtype)); 

   /*************************Start Timer***************************/
   MPI_Barrier (MPI_COMM_WORLD);
   elapsed_time = - MPI_Wtime();

   /* Dispatch b within processors */
   sqrt_p = (int)sqrt((double)p);
   is_square = (sqrt_p * sqrt_p == p);
   if (is_square) {
      /* Proc at (i,0) sends subvector to proc at (0,i) */
      /* Proc at (0,0) does a copy. */
      if ((grid_coords[0] == 0) && (grid_coords[1] == 0)) {
         for (i = 0; i < b_length; i++) btrans[i] = b[i];
      } 
      else if ((grid_coords[0] > 0) && (grid_coords[1] == 0)) {
         send_length = BLOCK_SIZE(grid_coords[0],grid_size[0],n);
         coords[0] = 0;
         coords[1] = grid_coords[0];
         MPI_Cart_rank(grid_comm, coords, &dest);
         MPI_Send (b, send_length, mpitype, dest, 0, grid_comm);
      } 
      else if ((grid_coords[1] > 0) && (grid_coords[0] == 0)) {
         coords[0] = grid_coords[1];
         coords[1] = 0;
         MPI_Cart_rank(grid_comm,coords, &src);
         MPI_Recv (btrans, b_length, mpitype, src, 0, grid_comm, &status);
      }
   } 
   else {
      /* Process at (0,0) gathers vector elements from procs in column 0 */
      if (grid_coords[1] == 0) {
         recv_cnt = (int *) malloc (grid_size[0] * sizeof(int));
         recv_disp = (int *) malloc (grid_size[0] * sizeof(int));
         for (i = 0; i < grid_size[0]; i++)
            recv_cnt[i] = BLOCK_SIZE(i,grid_size[0],n);
         recv_disp[0] = 0;
         for (i = 1; i < grid_size[0]; i++)
            recv_disp[i] = recv_disp[i-1] + recv_cnt[i-1];
         MPI_Gatherv (b, BLOCK_SIZE(grid_coords[0], grid_size[0], n),
            mpitype, tmp, recv_cnt, recv_disp, mpitype, 0, col_comm);
      }

      /* Process at (0,0) scatters vector elements to row 0 procs */
      if (grid_coords[0] == 0) {
         if (grid_size[1] > 1) {
            send_cnt = (int *) malloc (grid_size[1] * sizeof(int));
            send_disp = (int *) malloc (grid_size[1] * sizeof(int));
            for (i = 0; i < grid_size[1]; i++) {
               send_cnt[i] = BLOCK_SIZE(i,grid_size[1],n);
            }
            send_disp[0] = 0;
            for (i = 1; i < grid_size[1]; i++) {
               send_disp[i] = send_disp[i-1] + send_cnt[i-1];
            }
            MPI_Scatterv (tmp, send_cnt, send_disp, mpitype, btrans,
               b_length, mpitype, 0, row_comm);
         } 
         else { // Process (0, 0)
            for (i = 0; i < n; i++) btrans[i] = tmp[i];
         }
      }
   }
   /* Row 0 processors broadcast their subvectors to processors in same column */
   MPI_Bcast (btrans, b_length, mpitype, 0, col_comm);

   /* Matrix times vector */
   for (i = 0; i < rows; i++) {
      c_block[i] = 0.0;
      for (j = 0; j < cols; j++) {
         c_block[i] += a[i][j] * btrans[j];
      }
   }
   MPI_Reduce(c_block, c_sums, rows, mpitype, MPI_SUM, 0, row_comm);
   if (grid_coords[1] == 0) {
      print_block_vector (c_sums, mpitype, n, col_comm);
   }

   /*************************Stop Timer***************************/
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time += MPI_Wtime();
   if(!id){
      printf ("Matrix = %dx%d\nProcesses = %d\nTime = %.6f sec\n",
         n, n, p, elapsed_time);
   }

   MPI_Finalize();

   return 0;
}

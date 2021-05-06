#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
//gcc Rep_Sum_demo.c -o Rep_Sum_demo -lm -lmpfr -lgmp -O3

#define PRECISION            (700)
#define EPS_DOUBLE 2.0e-52


//计算时间
static struct timeval start;
static struct timeval end;
void tic(void){
  gettimeofday( &start, NULL );
}
double toc(void){
  gettimeofday( &end, NULL );
  return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

void MPI_Reduce_Reprosum(double *in, double *out, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm communicator, int k){
  double max_in[count];
  int my_rank, my_size;
  MPI_Comm_rank(communicator, &my_rank);
  MPI_Comm_size(communicator, &my_size);

  MPI_Allreduce(in, max_in, count, MPI_DOUBLE, MPI_MAX, communicator);

  double q1[count], q2[count], M1[count], M2[count];
  double local_high[count], local_low[count], x_low[count];
  double res_high[count], res_low[count];
  for(int i=0; i < count; i++){
    q1[i] = my_size*(max_in[i])/(1 - 2 * my_size * EPS_DOUBLE);
    M1[i] = pow(2, ceil(log2(q1[i])));
    q2[i] = my_size*( 2*EPS_DOUBLE *M1[i]) / (1 - 2*my_size*EPS_DOUBLE);
    M2[i] = pow(2, ceil(log2(q2[i+1])));
    local_high[i] = M1[i] + in[i] - M1[i];
    x_low[i] = in[i] - local_high[i];
    local_low[i] = M2[i] + x_low[i] - M2[i];
  }

  MPI_Reduce(local_high, res_high, count, datatype, op, root, communicator);
  MPI_Reduce(local_low,  res_low, count, datatype, op, root, communicator);
  for(int i=0; i < count; i++){
    out[i] = res_high[i] + res_low[i];
  }
}

int main(int argc, char** argv){
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int n = 200;
  //n = n + size - (n % size);
  int local_n = n / size;
  double elapsed_time;
  double *x = NULL;
  double *x_shuffled = NULL;
  double *local_x = malloc(local_n * sizeof(double));
  double *local_x_shuffled = malloc(local_n * sizeof(double));

  if(rank == 0){
    x = malloc(n * sizeof(double));
    x_shuffled = malloc(n * sizeof(double));

    // Set x to be a sine wave
    for(int i = 0; i < n; i++){
      x[i] = sin(2 * M_PI * (i / (double)n - 0.5));
    }
    // Shuffle x into x_shuffled
    for(int i = 0; i < n; i++){
      x_shuffled[i] = x[i];
    }
    double t;
    int r;
    for(int i = 0; i < n; i++){
      r = rand();
      t = x_shuffled[i];
      x_shuffled[i] = x_shuffled[i + (r % (n - i))];
      x_shuffled[i + (r % (n - i))] = t;
    }
  }

  MPI_Scatter(x, local_n, MPI_DOUBLE, local_x, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(x_shuffled, local_n, MPI_DOUBLE, local_x_shuffled, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(rank == 0){
    printf("Sum of sin(2* M_PI * (i / (double)n - 0.5)).  n = %d.\n\n", n);
  }

  double sum, local_sum, repro_sum;

  for(int i = 0; i < local_n; i++){
    local_sum += local_x[i];
  }
  MPI_Reduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce_Reprosum(&local_sum, &repro_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, 2);


  double sum_shuffled, repro_sum_shuffled;
  for(int i = 0; i < local_n; i++){
    local_sum += local_x_shuffled[i];
  }

  MPI_Reduce(&local_sum, &sum_shuffled, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(rank == 0){
    printf("MPI_Reduce: |%.17e - %.17e| = %g\n", sum, sum_shuffled, fabs(sum - sum_shuffled));
    printf("MPI_Reduce_Reprosum: |%.17e - %.17e| = %g\n",  repro_sum, repro_sum_shuffled, fabs(repro_sum - repro_sum_shuffled));
  }

  if(rank == 0){
    free(x);
    free(x_shuffled);
  }
  free(local_x);
  free(local_x_shuffled);
  MPI_Finalize();
}

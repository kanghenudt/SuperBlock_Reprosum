#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include <sys/time.h>
#include "MPI_ACCU_REDUCE.h"

//mpicc -o test_MPI_ACCU_REDUCE test_MPI_ACCU_REDUCE.c -L. -lMPI_ACCU_REDUCE -lm
//mpicc -c MPI_ACCU_REDUCE.c
//ar -crv libMPI_ACCU_REDUCE.a MPI_ACCU_REDUCE.o

/*MPI_ACCU_REDUCE函数定义为:double MPI_ACCU_REDUCE(void *sendbuf,void *recvbuf,int count,int optype, int root, MPI_Comm comm);
 IN   sendbuf    发送消息缓冲区的起始地址
 OUT  recvbuf    接收消息缓冲区中的地址
 IN   count      发送消息缓冲区中的数据个数（整形）
 IN   optype     归约操作符，在MPI_ACC_REDUCE.c文件中定义了三种高精度操作，ddsum,ddprod,norm，其optype分别为0，1，2（整形）
 IN   root       根进程序列号（整形）
 IN   comm       通信域（句柄）*/

static struct timeval start;
static struct timeval end;

void tic(void){
  gettimeofday( &start, NULL );
}

double toc(void){
  gettimeofday( &end, NULL );

  return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

int main(int argc, char ** argv)
{
  int size,rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int n = 25;
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
      printf("x[%d] is %.17e\n",i,x[i]);
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
    printf("%15s : Time (s) : |Sum - Sum of Shuffled| = ?\n", "Sum Method");
  }

  double local_sum = 0;
  double sum = 0;
  double repro_sum=0;
  double Kfoldrepro_sum=0;
  tic();
  for(int i = 0; i < local_n; i++){
    local_sum += local_x[i];
  }
  MPI_Reduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_ACCU_REDUCE(&local_sum, &repro_sum, 1, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD, 0);
  MPI_ACCU_REDUCE(&local_sum, &Kfoldrepro_sum, 1, MPI_DOUBLE, 4, 0, MPI_COMM_WORLD, 5);  

  elapsed_time = toc();

  double sum_shuffled = 0;
  double repro_sum_shuffled = 0;
  double Kfoldrepro_sum_shuffled = 0; 
  local_sum = 0;
  for(int i = 0; i < local_n; i++){
    local_sum += local_x_shuffled[i];
  }

  MPI_Reduce(&local_sum, &sum_shuffled, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_ACCU_REDUCE(&local_sum, &repro_sum_shuffled, 1, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD, 0);
  MPI_ACCU_REDUCE(&local_sum, &Kfoldrepro_sum_shuffled, 1, MPI_DOUBLE, 4, 0, MPI_COMM_WORLD, 5);  

  if(rank == 0){
    printf("MPI_Reduce: |%.17e - %.17e| = %g\n", sum, sum_shuffled, fabs(sum - sum_shuffled));
    printf("MPI_REPRO_REDUCE: |%.17e - %.17e| = %g\n",  repro_sum, repro_sum_shuffled, fabs(repro_sum - repro_sum_shuffled));
    printf("MPI_REPRO_ACCU_REDUCE: |%.17e - %.17e| = %g\n",  Kfoldrepro_sum, Kfoldrepro_sum_shuffled, fabs(Kfoldrepro_sum - Kfoldrepro_sum_shuffled));

  }

  if(rank == 0){
    free(x);
    free(x_shuffled);
  }
  free(local_x);
  free(local_x_shuffled);
  MPI_Finalize();
}

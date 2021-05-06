#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <time.h>

#define EPS_DOUBLE 2.0e-52
#define min(a,b) ((a) < (b) ? (a) : (b))

static struct timeval start;
static struct timeval end;
void tic(void){
  gettimeofday( &start, NULL );
}
double toc(void){
  gettimeofday( &end, NULL );
  return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

double OMP_Kfold_ReproSum(double *x, int n, int K, int thread_count){
  double  res = 0.0;
  double q1, M1, q2, M2;
  double max_num = fabs(x[0]);
  #pragma omp parallel for
  for(int i =1; i < n; i++){
    if(max_num < fabs(x[i])){
      max_num = fabs(x[i]);
    }
  }
  q1  = n*max_num/(1 - 2 * n * EPS_DOUBLE);
  M1 = pow(2, ceil(log2(q1)));
  q2 = n*( 2*EPS_DOUBLE *M1) / (1 - 2*n*EPS_DOUBLE);
  M2 = pow(2, ceil(log2(q2)));
  double res_high, res_low, temp_T_high, temp_T_low, temp_x;
  res_high = 0.0;
  res_low  = 0.0;
  #pragma omp parallel for private(temp_T_high, temp_T_low, temp_x) reduction(+:res_high, res_low)
  for(int j = 0; j < n; j++){
    temp_T_high = M1 + x[j] - M1;
    temp_x = x[j] - temp_T_high;
    res_high = res_high + temp_T_high;
    temp_T_low = M2 + temp_x - M2;
    res_low = res_low + temp_T_low;
  }
  res = res_high + res_low;
  return res;
}

double OMP_SuperBlock_reprosum(double *x , int n, int block, int thread_count){
  double max_num = fabs(x[0]);
  double q1, M1, q2, M2;
  #pragma omp parallel for
  for(int i =1; i < n; i++){
    if(max_num < fabs(x[i])){
      max_num = fabs(x[i]);
    }
  }
  q1  = n*max_num/(1 - 2 * n * EPS_DOUBLE);
  M1 = pow(2, ceil(log2(q1)));
  q2 = n*( 2*EPS_DOUBLE *M1) / (1 - 2*n*EPS_DOUBLE);
  M2 = pow(2, ceil(log2(q2)));

  double s_high, s_low, sres_high, sres_low, res_high, res_low, s;
  double temp_T_high, temp_T_low, temp_x;
  int ibeg,iend;
  int nblks = ceil((double)n/block);
  int nsblks = ceil(sqrt((double)nblks));
  int blksInSblk = ceil((double)nblks/nsblks);
  s_high=0.0;
  s_low=0.0;
  for(int k = 0; k < nsblks; k++){
    sres_high = 0.0;
    sres_low=0.0;
    for(int b = 0; b < blksInSblk; b++){
      ibeg = b*block+k*blksInSblk*block;
      iend = min((b+1)*block+k*blksInSblk*block, n);
      res_high = 0.0;
      res_low  = 0.0;
      #pragma omp parallel for private(temp_T_high, temp_T_low, temp_x) reduction(+:res_high, res_low)
      for(int i = ibeg; i < iend; i++){
        temp_T_high = M1 + x[i] - M1;
        temp_x = x[i] - temp_T_high;
        res_high = res_high + temp_T_high;
        temp_T_low = M2 + temp_x - M2;
        res_low = res_low + temp_T_low;
      }
      sres_high = sres_high + res_high;
      sres_low = sres_low + res_low;
    }
    s_high = s_high + sres_high;
    s_low = s_low + sres_low;
  }
  s = s_high + s_low;
  return s;
}

int main(int argc, char *argv[]){
  int n = 700000;
  int nb = 2048;
  int thread_count = strtol(argv[1], NULL, 10);
  double *x = malloc(n * sizeof(double));
  double *x_shuffled = malloc(n * sizeof(double));

  srand(time(0));
  for(int i=0; i<size; i++){
    x[i] = 0+1.0*(rand()%RAND_MAX)/RAND_MAX *(1-0);	//设为RAND_MAX,随机效果更好
  }

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

  printf("----------n is %d----------\n", n);
  double sum = 0.0, sum_shuffled = 0.0;
  #pragma omp parallel for num_threads(thread_count) reduction(+:sum)
  for(int i=0; i <n; i++){
    sum += x[i];
  }
  #pragma omp parallel for num_threads(thread_count) reduction(+:sum_shuffled)
  for(int i=0; i <n; i++){
    sum_shuffled += x_shuffled[i];
  }
  printf("OpenMP reduction sum: |%.17e - %.17e| = %g\n",  sum, sum_shuffled, fabs(sum - sum_shuffled));

  sum = OMP_SuperBlock_reprosum(x, n, nb, thread_count);
  sum_shuffled = OMP_SuperBlock_reprosum(x_shuffled, n, nb, thread_count);
  printf("OpenMP SuperBlock_reprosum: |%.17e - %.17e| = %g\n",  sum, sum_shuffled, fabs(sum - sum_shuffled));

  sum = OMP_Kfold_ReproSum(x, n, 2, thread_count);
  sum_shuffled = OMP_Kfold_ReproSum(x_shuffled, n, 2, thread_count);
  printf("OpenMP OMP_Kfold_ReproSum: |%.17e - %.17e| = %g\n",  res, res_shuffled, fabs(res - res_shuffled));

  free(x);
  free(x_shuffled);
  return 0;
}


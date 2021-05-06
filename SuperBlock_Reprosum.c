#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
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

double Kfold_ReproSum(double *x, int n, int K){
    double  res = 0.0;
    double q1, M1, q2, M2;
    double max_num = fabs(x[0]);
    for(int i =1; i < n; i++){
        if(max_num < fabs(x[i])){
        max_num = fabs(x[i]);
        }
    }
    q1  = n*max_num/(1 - 2 * n * EPS_DOUBLE);
    M1 = pow(2, ceil(log2(q1)));
    q2 = n*( 2*EPS_DOUBLE *M1) / (1 - 2*n*EPS_DOUBLE);
    M2 = pow(2, ceil(log2(q2)));
    double res_high, res_low, temp_T_high, temp_T_low, temp_x_low;
    res_high = 0.0;
    res_low  = 0.0;
    for(int j = 0; j < n; j++){
        temp_T_high = M1 + x[j] - M1;
        temp_x_low = x[j] - temp_T_high;
        res_high = res_high + temp_T_high;
        temp_T_low = M2 + temp_x_low - M2;
        res_low = res_low + temp_T_low;
    }
    res = res_high + res_low;
    return res;
}

double SuperBlock_Reprosum(double *x , int n, int b){
    double q1, M1, q2, M2;
    double max_num = fabs(x[0]);
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
    double temp_T_high, temp_T_low, temp_x_low;
    int ibeg,iend;
    int nblks = ceil((float)n/b);
    int nsblks = ceil(sqrt(nblks));
    int blksInSblk = ceil((float)nblks/nsblks);
    s_high=0.0;
    s_low=0.0;
    for(int s = 0; s < nsblks; s++){
        sres_high = 0.0;
        sres_low=0.0;
        for(int k = 0; k < blksInSblk; k++){
            ibeg = k*b+s*blksInSblk*b;
            iend = min((k+1)*b+s*blksInSblk*b,n);
            res_high = 0.0;
            res_low  = 0.0;
            for(int j = ibeg; j < iend; j++){
                temp_T_high = M1 + x[j] - M1;
                temp_x_low = x[j] - temp_T_high;
                res_high = res_high + temp_T_high;
                temp_T_low = M2 + temp_x_low - M2;
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


int main(){
    int size = 1200000;
    int block = 2048;

    double *x = malloc(size * sizeof(double));
    double *x_shuffled = malloc(size* sizeof(double));
    double *y = malloc(size * sizeof(double));
    double *y_shuffled = malloc(size* sizeof(double));

    srand(time(0));
    for(int i=0; i<size; i++){
        x[i] = 0+1.0*(rand()%RAND_MAX)/RAND_MAX *(1-0);	//设为RAND_MAX,随机效果更好
        y[i] = x[i];
    }
    for(int i = 0; i < size; i++){
        x_shuffled[i] = x[i];
    }
    double t;
    int r;
    for(int i = 0; i < size; i++){
        r = rand();
        t = x_shuffled[i];
        x_shuffled[i] = x_shuffled[i + (r % (size - i))];
        x_shuffled[i + (r % (size - i))] = t;
        y_shuffled[i] = x_shuffled[i];
    }
    printf("%18s : |Sum - Sum of Shuffled| = ?\n", "Sum Method");

    double sum = 0.0,sum_shuffled = 0.0;
    double elapsed_time;
    tic();
    for(int i = 0; i < size; i++){
        sum += x[i];
    }
    elapsed_time = toc();
    for(int i = 0; i < size; i++){
        sum_shuffled += x_shuffled[i];
    }
    printf("%15s : %-8g : |%.17e - %.17e| = %g\n", "sum", elapsed_time, sum, sum_shuffled, fabs(sum - sum_shuffled));

    tic();
    sum = Kfold_ReproSum(x, size, 2);
    elapsed_time = toc();
    sum_shuffled = Kfold_ReproSum(x_shuffled, size, 2);
    printf("%15s : %-8g : |%.10e - %.10e| = %g\n", "Kfold_ReproSum", elapsed_time, sum, sum_shuffled, fabs(sum - sum_shuffled));

    tic();
    sum = SuperBlock_Reprosum(y, size, block);
    elapsed_time = toc();
    sum_shuffled = SuperBlock_Reprosum(y_shuffled, size, block);
    printf("%15s : %-8g : |%.10e - %.10e| = %g\n", "SuperBlock_Reprosum", elapsed_time, sum, sum_shuffled, fabs(sum - sum_shuffled));




    free(x);
    free(x_shuffled);
    free(y);
    free(y_shuffled);

    return 0;

}


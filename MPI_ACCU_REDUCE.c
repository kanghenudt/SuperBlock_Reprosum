#include "MPI_ACCU_REDUCE.h"

dd_real  quick_two_sum(double a, double b){
  dd_real res;
  res.hi = a + b;
  res.lo = b - (res.hi - a);
  return res;
}

dd_real two_sum(double a, double b){
  dd_real res;
  res.hi = a + b;
  double bb = res.hi - a;
  res.lo = (a - (res.hi - bb)) + (b -bb);
  return res;
}

//将浮点数分为两部分的无误差变换算法，将一个精度为p的浮点数分成两个精度最多为s-1位的浮点数
dd_real split(double a) {
  double temp;
  dd_real res;
  temp = SPLITTER * a;//SPLITTER是一个即常数
  res.hi = temp - (temp - a);
  res.lo = a - res.hi;
  return res;
}

//TwoProd()算法,两个浮点数乘法的无误差变换
dd_real two_prod(double a, double b) {
  dd_real res,aa,bb;
  res.hi = a * b;
  aa=split(a);
  bb=split(b);
  res.lo = ((aa.hi * bb.hi - res.hi) + aa.hi * bb.lo + aa.lo * bb.hi) + aa.lo * bb.lo;
  return res;
}

//两个double-double格式数 a=(ah,al)和b=(bh,bl)的乘法运算
dd_real prod_dd_dd(const dd_real aa, const dd_real bb) {
  dd_real res,temp;
  temp=two_prod(aa.hi,bb.hi);
  temp.lo += (aa.lo * bb.hi + aa.hi * bb.lo);
  res=quick_two_sum(temp.hi, temp.lo);
  return res;
}

//两个double-double格式数 a=(ah,al)和b=(bh,bl)的加法运算
dd_real add_dd_dd(dd_real aa, dd_real bb){
  dd_real res,temp;
  temp=two_sum(aa.hi,bb.hi);
  temp.lo+=(aa.lo+bb.lo);
  res=quick_two_sum(temp.hi, temp.lo);
  return res;
}

void myddSum(dd_real *in, dd_real *inout, int * len, MPI_Datatype *dptr){
  int i;
  dd_real temp;
  for(i = 0; i < *len ; i++){
    temp = add_dd_dd(*inout, *in);
    *inout = temp;
    in++;
    inout++;
  }
}

void myddprod(dd_real *in, dd_real *inout, int * len, MPI_Datatype *dptr){
  int i;
  dd_real temp;
  for(i = 0; i < *len; i++){
    temp = prod_dd_dd(*inout, *in);
    *inout = temp;
    in++;
    inout++;
  }
}
double normNumber(double x){
  dd_real temp;
  double res;
  temp = two_prod(x,x);
  res = temp.hi+temp.lo;
  return res;
}
/*double *normNumber(double *x, unsigned int n)
{
    dd_real temp1,temp2,temp3;
    double p,h,l,res;
    double *r;
    int i;
    int m = n*2;
    r = (double *)malloc(sizeof(double)* m);
    temp1 = two_prod(x[0],x[0]);
    p = temp1.hi;
    r[0] = temp1.lo;
    for(i=1 ;i < n ; i++)
    {
        temp2 = two_prod(x[i],x[i]);
        h = temp2.hi;
        r[i] = temp2.lo;
        temp3 = two_sum(p,h);
        p = temp3.hi;
        r[n+i-1] = temp3.lo;
    }
    r[m-1] = p;
    return r;
	
}*/

void mynorm(double *in, double *inout, int *len, MPI_Datatype *dptr){
  int i;
  dd_real temp1;   
  for(i=0; i < *len; i++) {
    temp1 = two_sum(*inout,*in);
    *inout = temp1.hi+temp1.lo;
    in++;
    inout++;
  }
}

void MPI_REPRO_REDUCE(double *in,double *inout, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm){  
  double  *max_in = malloc(count * sizeof(double));
  int my_rank, my_size;
  MPI_Comm_rank(comm, &my_rank);
  MPI_Comm_size(comm, &my_size);
  double *abs_in = malloc(count * sizeof(double));
  for(int i=0; i< count;i++){
     abs_in[i] = fabs(in[i]);
  } 
  MPI_Allreduce(abs_in, max_in, count , MPI_DOUBLE, MPI_MAX,  comm);
  double dp = 2.0e-52;
  double  *oq = malloc(count * sizeof(double));
  double  *M = malloc(count * sizeof(double));
  double *q = malloc(count * sizeof(double));
  for(int i=0; i< count;i++){
    oq[i] = my_size*(max_in[i])/(1-2*my_size*dp);
    M[i] = pow(2, ceil(log(oq[i])/log(2)));
    q[i] = M[i] + in[i] - M[i];
  } 
  MPI_Reduce(q, inout, count,  datatype, op, root, comm );
  free(oq);
  free(M);
  free(q);
  free(max_in);
  free(abs_in);  
}

dd_real ExtractVector(double M, double v){
    dd_real temp;
    double q;
    q = M + v - M;
    temp.lo = v - q;
    temp.hi = q;
    return temp;
}

void MPI_ACCU_REPRO_REDUCE(double *in, double *out, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm communicator, int k){  
  double  *max_in = malloc(count * sizeof(double));
  int my_rank, my_size;
  MPI_Comm_rank(communicator, &my_rank);
  MPI_Comm_size(communicator, &my_size);

  MPI_Allreduce(in, max_in, count, MPI_DOUBLE, MPI_MAX, communicator);

  double dp = 2.0e-52;
  double **oq = (double **)malloc(count *sizeof(double *));
  double **M  = (double **)malloc(count *sizeof(double *));
  double **T  = (double **)malloc(count *sizeof(double *));
  for(int i=0; i < count; i++){
    oq[i] =  (double *)malloc(k *sizeof(double));
    M[i]  =  (double *)malloc(k *sizeof(double));
    T[i]  =  (double *)malloc(k *sizeof(double));    
  }

  for(int i=0; i < count; i++){
    oq[i][0] = my_size*(max_in[i])/(1 - 2 * my_size * dp);
    M[i][0] = pow(2, ceil(log(oq[i][0])/log(2)));
  }
  dd_real temp;
  for(int j=0; j < count; j++){
    for(int i=0; i < k-1; i++){
      temp = ExtractVector(M[j][i], in[j]);
      T[j][i] = temp.hi;
      in[j] = temp.lo;
      oq[j][i+1] = my_size * (2 * dp * M[j][i]) /(1 - 2 * my_size * dp);
      M[j][i+1] = pow(2, ceil(log(oq[j][i+1]) / log(2)));  
    }
  }
  for(int j=0; j < count; j++){ 
    T[j][k-1] = 0;
  }
  double *a = (double *)malloc(count * sizeof(double));
  for(int j=0; j < count; j++){ 
    a[j] = (M[j][k-1] + in[j]) - M[j][k-1];
    T[j][k-1] = a[j];
  }
  double *res        = (double *)malloc(count *sizeof(double));
  double **local_res = (double **)malloc(count *sizeof(double *));
  for(int i=0; i < count; i++){
    local_res[i] =  (double *)malloc(k *sizeof(double));
  }
  for(int j=0; j < count; j++){ 
    MPI_Reduce(T[j], local_res[j], k, datatype, op, root, communicator);
  }
  for(int j=0; j < count; j++){ 
    for(int i=0;i<k;i++){
      res[j] = res[j] +local_res[j][i];
    }
  }
  *out = *res;
  free(max_in);
  free(res); 
  free(a);
  for(int i=0; i < count; i++){
    free(oq[i]);
    free(M[i]);
    free(T[i]);
    free(local_res[i]);
  }
  free(oq);
  free(M);
  free(T);
  free(local_res);
}
void MPI_ACCU_REDUCE(void *in, void *inout, int count, MPI_Datatype datatype, int optype, int root, MPI_Comm comm, int k){
  dd_real temp;
  MPI_Op DDSUM,DDPROD,NORM;
  MPI_Datatype ctype;

  MPI_Type_contiguous(2, MPI_DOUBLE, &ctype);
  MPI_Type_commit(&ctype);


  MPI_Op_create((MPI_User_function*)myddSum, 1, &DDSUM);
  MPI_Op_create((MPI_User_function*)myddprod, 1, &DDPROD);
  MPI_Op_create((MPI_User_function*)mynorm, 1, &NORM);

  switch(optype){
    case 0:
    MPI_Reduce(in, inout, count, ctype, DDSUM, root, comm);
    break;
    case 1:
    MPI_Reduce(in, inout, count, ctype, DDPROD, root, comm);
    break;
    case 2:
    MPI_Reduce(in, inout, count, datatype, NORM, root, comm);
    break;
    case 3:
    MPI_REPRO_REDUCE(in, inout, count, datatype, MPI_SUM, root, comm);
    break;
    case 4:
    MPI_ACCU_REPRO_REDUCE(in, inout, count, datatype, MPI_SUM, root, comm, k);
    break;
    default:
    printf("error\n");
    break;
  }
  MPI_Op_free(&DDSUM);
  MPI_Op_free(&DDPROD);
  MPI_Op_free(&NORM);
}


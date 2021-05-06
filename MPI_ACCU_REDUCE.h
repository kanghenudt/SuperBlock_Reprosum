
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#ifndef _MPI_ACCU_REDUCE
#define _MPI_ACCU_REDUCE

#ifndef SPLITTER
	#define SPLITTER 134217729.0               // = 2^27 + 1
#endif

#ifndef unit
	#define unit  1.110223024625157e-16
#endif

typedef struct doubledouble{
    double hi;
    double lo;
}dd_real;

dd_real  quick_two_sum(double a, double b);
dd_real two_sum(double a, double b);
dd_real split(double a);
dd_real two_prod(double a, double b);
dd_real ExtractVector(double M, double v);

dd_real add_dd_dd(const dd_real aa, const dd_real bb);
dd_real prod_dd_dd(const dd_real aa, const dd_real bb);
double normNumber(double x);

void myddSum(dd_real *in, dd_real *inout, int * len, MPI_Datatype *dptr);
void myddprod(dd_real *in, dd_real *inout, int * len, MPI_Datatype *dptr);
void mynorm(double *in, double *inout, int * len, MPI_Datatype *dptr);

void MPI_REPRO_REDUCE(double *in,double *inout, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

void MPI_ACCU_REPRO_REDUCE(double *in, double *out, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, int k);

void MPI_ACCU_REDUCE(void *in, void *inout, int count, MPI_Datatype datatype, int optype, int root, MPI_Comm comm, int k);
#endif

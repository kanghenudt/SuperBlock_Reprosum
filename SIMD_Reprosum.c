#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//------------------------------------------------------------------------------------------//
#ifndef EPS_DOUBLE
#define EPS_DOUBLE 2.0e-52
#endif

#define N		"x0"	/* vector length */
#define X		"x1"	/* "X" vector address */
#define INC_X		"x2"	/* "X" stride */
#define M		"x3"	/* vector length */
#define J		"x5"	/* loop variable */
#define REG0		"xzr"
#define SUMR            "d0"

#define MAXF	"d0"
#define TMPF	"d1"
#define TMPVF	"{v1.d}[0]"
#define SZ	8

#define MAXINIT_F1                                               \
    "ldr         "MAXF", ["X"], #SZ                              \n"

#define MAXKERNEL_F1                                             \
    "ldr         "TMPF", ["X"], #SZ                              \n"   \
    "fcmp	 "MAXF", "TMPF"                                  \n"   \
    "fcsel	 "MAXF", "MAXF", "TMPF", ge                    \n"

#define MAXINIT_F4                                                \
    "ld2	{v0.2d,v1.2d}, ["X"], #32                         \n"   \
    "fmax	v0.2d, v0.2d, v1.2d                               \n"   \
    "fmaxp	"MAXF", v0.2d                                     \n"

#define MAXKERNEL_F4                                               \
    "ld2	{v1.2d,v2.2d}, ["X"], #32                            \n"   \
    "fmax	v1.2d, v1.2d, v2.2d                                \n"   \
    "fmaxp	"TMPF", v1.2d                                        \n"   \
    "fcmp	"MAXF", "TMPF"                                        \n"   \
    "fcsel	"MAXF", "MAXF", "TMPF", ge                            \n"

#define MAXINIT_S                                                   \
    "lsl	"INC_X", "INC_X", #3                                   \n"    \
    "ld1	{v0.d}[0], ["X"], "INC_X"                              \n"

#define MAXKERNEL_S1                                                \
    "ld1	"TMPVF", ["X"], "INC_X"                                   \n"   \
    "fcmp	"MAXF", "TMPF"                                         \n"   \
    "fcsel	"MAXF", "MAXF", "TMPF", ge                              \n"

#define SUMR_F8                                                                    \
	"mov	 "M", %[M_]			                                   \n"   \
    "ld1     {v2.2d, v3.2d, v4.2d, v5.2d}, ["X"]    \n"   \
    "PRFM    PLDL1KEEP, ["X", #1024]                         \n"    \
    "PRFM    PLDL1KEEP, ["X", #1024+64]                         \n"    \
	"fadd    v6.2d, v2.2d, "M"                     \n"     \
	"fsub    v7.2d, v6.2d, "M"                     \n"     \
	"fsub    v8.2d, v2.2d, v7.2d                   \n"     \
    "str        v8.2d, ["X"]                                   \n"  \
	"fadd    v9.2d, v3.2d, "M"                     \n"     \
	"fsub    v10.2d, v9.2d, "M"                     \n"     \
	"fsub    v11.2d, v3.2d, v10.2d                   \n"     \
    "add	"X", "X", #16			\n"	\
    "str        v11.2d, ["X"]                                  \n"  \
	"fadd    v12.2d, v4.2d, "M"                     \n"     \
	"fsub    v13.2d, v12.2d, "M"                     \n"     \
	"fsub    v14.2d, v4.2d, v13.2d                   \n"     \
    "add	"X", "X", #16			\n"	\
    "str        v14.2d, ["X"]                                  \n"  \
	"fadd    v15.2d, v5.2d, "M"                     \n"     \
	"fsub    v16.2d, v15.2d, "M"                     \n"     \
	"fsub    v17.2d, v5.2d, v16.2d                   \n"     \
    "add	"X", "X", #16			\n"	\
    "str        v17.2d, ["X"]                                  \n"  \
    "add	"X", "X", #16			\n"	\
	"fadd    v15.2d, v5.2d, "M"                     \n"     \
	"fsub    v16.2d, v15.2d, "M"                     \n"     \
	"fsub    v17.2d, v5.2d, v16.2d                   \n"     \
    "add	"X", "X", #8			\n"	\
    "str        v17.2d, ["X"]                                  \n"  \
    "add	"X", "X", #8			\n"	\
	"fadd    v7.2d, v10.2d, v7.2d                 \n"     \
	"fadd    v13.2d, v16.2d, v13.2d                 \n"     \
    "fadd    v31.2d, v13.2d, v7.2d                 \n"

#define SUMR_F8_FINALIZE                                       \
        "fadd   v0.2d, v0.2d, v31.2d                    \n"    \
        "faddp  "SUMR", v0.2d                               \n"

#define SUMR_F1                                       \
	"mov	 "M", %[M_]			               \n"  \
    "ldr     "TMPF", ["X"]                       \n"  \
	"fadd    d2, "TMPF", "M"               \n"  \
	"fsub    d2, d2, "M"                          \n"  \
	"fsub    d3, "TMPF", d2                  \n"  \
    "str        d3, ["X"]                               \n"  \
	"add	"X", "X", #8		       	           \n"	\
    "fadd    d31, d31, d2                          \n"

#define SUMR_F1_FINALIZE                			\
    "fadd    "SUMR" , "SUMR" , d31         \n"

double ReproSum(double *x, int n, int K, double maxf){
    double sumr = 0.0;
    double T, q, M;

    q = n * maxf/(1 - 2 * n * EPS_DOUBLE);
    M = pow(2, ceil(log(q) / log(2)));
    for(int i=0; i < K; i++ ){
        __asm__ __volatile__ (
	        "	            mov	            "N",  %[N_]			                       \n"
            "	            mov	            "X",  %[X_]			                        \n"
            "	            mov	            "INC_X", %[INCX_]		              \n"
            "	            fmov	        "SUMR", "REG0"			               \n"
            "	            fmov	        d1, "REG0"			\n"
            "	            fmov	        d2, "REG0"			\n"
            "	            fmov	        d3, "REG0"			\n"
            "	            fmov	        d4, "REG0"			\n"
            "	            fmov	        d5, "REG0"			\n"
            "	            fmov	        d6, "REG0"			\n"
            "	            fmov	        d7, "REG0"			\n"
            "	            cmp              "N", xzr                                          \n"
            "	            ble              Lsumr_L999                                \n"
            "Lsum_F_BEGIN:                                                           \n"
            "	            asr             "J", "N", #3                                   \n"
            "	            cmp         "J", xzr                                            \n"
            "	            beq          Lsumr_F1                                       \n"

            "Lsumr_F8:                                                                      \n"
            "	            "SUMR_F8"                                                     \n"
            "	            subs    "J", "J", #1                                         \n"
            "	            bne     Lsumr_F8                                            \n"
            "	            "SUMR_F8_FINALIZE"                                \n"

            "Lsumr_F1:                                                                       \n"
            "	            ands    "J", "N", #7                                        \n"
            "	            ble     Lsumr_L999                                         \n"
            "	            fmov    d31, xzr                                               \n"

            "Lsumr_F1_plus:                                                           \n"
            "	            "SUMR_F1"                                                     \n"
            "	            subs    "J", "J", #1                                        \n"
            "	            bne     Lsumr_F1_plus                                \n"
            "	            "SUMR_F1_FINALIZE"                                \n"

            "Lsumr_L999:				                                                  \n"
	        "	            fmov	%[SUMR_], "SUMR"		              \n"
        :
        : [SUMR_]  "r"  (&sumr),		//%0
          [N_]        "r"  (n),		          //%1
          [X_]         "r"  (x),		       //%2
          [INCX_] "r"  (inc_x),		  //%3
          [M] "r"  (M),		  //%3
        : "cc",
          "memory",
	        "x0", "x1", "x2", "x3", "x4", "x5",
	        "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7"
        );
        q = n * (2 * EPS_DOUBLE * M) / (1 - 2 * n * EPS_DOUBLE);
      	M = pow(2, ceil(log(q) / log(2)));
    }
    return sumr;
}

double SIMD_ReproSum(int n, double *x, int inc_x)
{
    int max_i=0;
    int max_ix=0;
    double maxf=0.0;

    double sumf = 0.0;
    if (n <= 0 || inc_x <= 0) return(sumf);

	__asm__ __volatile__ (
	"	        mov	    "N", %[N_]				\n"
	"	        mov	    "X", %[X_]				\n"
	"	        mov	    "INC_X", %[INCX_]			\n"

	"	        cmp	    "INC_X", #1				\n"
	"	        bne	      .Lmax_kernel_S_BEGIN		\n"

	".Lmax_kernel_F_BEGIN:			\n"
    "           asr	        "J", "N", #2                                     \n"
    "           cmp	     "J", xzr                                     \n"
    "           beq       .Lmax_kernel_F1_INIT    \n"
    "           "MAXINIT_F4"                                \n"
    "           subs	  "J", "J", #1                                     \n"
    "           beq	      .Lmax_kernel_F1          \n"

    ".Lmax_kernel_F4:                        \n"
    "           "MAXKERNEL_F4"                      \n"
    "           subs	    "J", "J", #1                  \n"
    "           bne	        .Lmax_kernel_F4              \n"

    ".Lmax_kernel_F1:                        \n"
    "           ands	  "J", "N", #3                  \n"
    "           ble	       .Lmax_kernel_L999              \n"

    ".Lmax_kernel_F10:         \n"
    "           "MAXKERNEL_F1"             \n"
    "           subs	    "J", "J", #1                  \n"
    "           bne         .Lmax_kernel_F10          \n"
    "	        mov	        %[MAXF_], "MAXF"			\n"

    ".Lmax_kernel_F1_INIT:              \n"
    "           "MAXKERNEL_F1"            \n"
    "           subs	    "N", "N", #1        \n"
    "           b	            .Lmax_kernel_F1     \n"

    ".Lmax_kernel_S_BEGIN:   \n"
    "           "MAXINIT_S"                   \n"
    "           subs	    "N", "N", #1        \n"
    "           ble	          .Lmax_kernel_L999    \n"
    "           asr            "J", "N", #2                    \n"
    "           cmp	        "J", xzr          \n"
    "           ble	         .Lmax_kernel_S1     \n"

    ".Lmax_kernel_S4:                     \n"
    "           "MAXKERNEL_S1"         \n"
    "           "MAXKERNEL_S1"         \n"
    "           "MAXKERNEL_S1"         \n"
    "           "MAXKERNEL_S1"         \n"
    "           subs	    "J", "J", #1                  \n"
    "           bne	        .Lmax_kernel_S4            \n"

    ".Lmax_kernel_S1:            \n"
    "           ands	    "J", "N", #3          \n"
    "           ble	           .Lmax_kernel_L999  \n"

    ".Lmax_kernel_S10:  \n"
    "           "MAXKERNEL_S1"         \n"
    "           subs	    "J", "J", #1                  \n"
    "           bne	        .Lmax_kernel_S10            \n"

    ".Lmax_kernel_L999:  \n"
    "	        mov	       %[MAXF_], "MAXF"			\n"

	: [MAXF_] "=r" (maxf)		//%0
	: [N_]    "r"  (n),		//%1
	  [X_]    "r"  (x),		//%2
	  [INCX_] "r"  (inc_x)		//%3
	: "cc",
	  "memory",
	  "x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7",
	  "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7"
	);

    int K = 2;
    sumf = ReproSum(x, n, K, maxf);
    return(sumf);
}

int main(){
    int size = 40000;

    double *x = malloc(size * sizeof(double));
    double *x_shuffled = malloc(size* sizeof(double));

    printf("Sum of sin(2* M_PI * (i / (double)size - 0.5)). size = %d\n", size);

    srand(time(0));
    for(int i=0; i<size; i++){
        x[i] = 0+1.0*(rand()%RAND_MAX)/RAND_MAX *(1-0);	//设为RAND_MAX,随机效果更好
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
    printf("%15s : %-8g : |%.17e - %.17e| = %g\n", "sum ", elapsed_time, sum, sum_shuffled, fabs(sum - sum_shuffled));

    tic();
    sum = SIMD_ReproSum(size , x , 1);
    elapsed_time = toc();
    sum_shuffled = SIMD_ReproSum(size , x_shuffled , 1);
    printf("%15s : %-8g : |%.17e - %.17e| = %g\n", "SIMD_ReproSum ", elapsed_time, sum, sum_shuffled, fabs(sum - sum_shuffled));

    free(x);
    free(x_shuffled);
    return 0;

}

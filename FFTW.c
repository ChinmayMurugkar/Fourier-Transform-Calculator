{\rtf1\ansi\ansicpg1252\cocoartf1265\cocoasubrtf190
{\fonttbl\f0\froman\fcharset0 Times-Roman;}
{\colortbl;\red255\green255\blue255;}
\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\deftab720
\pard\pardeftab720

\f0\fs24 \cf0 /*Parallel implementation using MPI and FFTW.\
\'a0* to compile on RCC using:\
\'a0* mpicc -o fftw fftw.c -lfftw_mpi -lfftw\
-L/usr/local/fftw-2.1.5_mpich_pgi72/lib\
-I/usr/local/fftw-2.1.5_mpich_pgi72/include\
\
\'a0* to RUN:\
\'a0* mpirun -np 2 fftw\
\
\'a0*/\
\
#include "mpi.h"\
#include "fftw_mpi.h"\
#include <stdio.h>\
#include <math.h>\
#include <time.h>\
\
#define pi 3.14\
#define h 3.14/1024\
\
int main ( int argc, char *argv[] ) \{\
\
\'a0 int myid, numprocs;\
\'a0 fftw_mpi_plan plan1, plan2; /*plans for forward & reverse transformations*/\
\
\'a0 const int N = 2048; /* size of vector to transform */\
\'a0 fftw_complex u2[N]; /* stores initial vector and transposed vector*/\
\'a0 fftw_complex A[N]; /* work space */\
\'a0 fftw_complex u222[N];\
\'a0 int i, j, k, n1;\
\
\'a0 int local_n, local_start;\
\'a0 int local_n_after, local_start_after, total_local_size;\
\
\'a0 double startwtime, endwtime;\
\'a0 time_t t1;\
\
\'a0 MPI_Init ( &argc, &argv );\
\'a0 MPI_Comm_size ( MPI_COMM_WORLD, &numprocs );\
\'a0 MPI_Comm_rank ( MPI_COMM_WORLD, &myid );\
\
\'a0 if ( myid == 0 )\
\'a0 \'a0 \'a0 \'a0 startwtime = MPI_Wtime();\
\
\'a0 /* create plans; estimate the best approach rather than measuring it in full*/\
\'a0 plan1 = fftw_mpi_create_plan ( MPI_COMM_WORLD, N, FFTW_FORWARD,\
\'a0 \'a0 FFTW_ESTIMATE );\
\'a0 plan2 = fftw_mpi_create_plan ( MPI_COMM_WORLD, N, FFTW_BACKWARD,\
\'a0 \'a0 FFTW_ESTIMATE );\
\
\'a0 /* figure out how this plan breaks up the data */\
\'a0 fftw_mpi_local_sizes ( plan1, &local_n, &local_start,\
\'a0 \'a0 &local_n_after, &local_start_after, &total_local_size );\
\
\'a0 /* get system time and use it to seed the random number generator */\
\'a0 (void) time(&t1);\
\
\
\'a0 /* initialize data to complex array with values in RHS of the equation ie F */\
\'a0 for ( j = 0; j < local_n; j++ )\
\'a0 \'a0\{\
\'a0 \'a0 c_re(u2[j]) = (2*(h*h)*cos(j));\
\'a0 \'a0// c_im(u2[j]) = 0;\
\'a0 \}\
\
\'a0 /*output current contents of this array */\
\'a0/*printf("*********************Initial array:********************\\n");\
\'a0 for ( j = 0; j < local_n; j++) \{\
\'a0 \'a0 printf("%i (on %i): (%f, %f)\\n", j+local_start, myid,\
\'a0 \'a0 \'a0 c_re(u2[j]), c_im(u2[j]));\
\
\'a0 \}\
\'a0 printf("\\n");\
*/\
\'a0 /*apply distributed FFTW*/\
\'a0 fftw_mpi( plan1, 1, u2, A);\
\
\'a0 /*output transformed array */\
\'a0/* printf("******************After transformation:********************\\n");\
\'a0 for ( j = 0; j < local_n_after; j++) \{\
\'a0 \'a0 printf("%i (on %i): (%f, %f)\\n", j+local_start_after, myid,\
\'a0 \'a0 \'a0 c_re(u2[j]), c_im(u2[j]));\
\'a0 \}\
\'a0 printf("\\n");\
*/\
\'a0 \'a0//calculating the equation U(k)\
\
\'a0 \'a0double temp;\
\
\'a0 \'a0for(j=0;j<local_n_after;j++)\
\'a0 \'a0 \{\
\
\'a0 \'a0 \'a0 \'a0temp= -2*cos((2*pi*j)/N)+2+(h*h);\
\'a0 \'a0 \'a0 \'a0c_re(u222[j])=c_re(u2[j])/temp;\
\'a0 \'a0 \}\
\
\
\'a0 /* divide real and imaginary parts by N to normalize before retransforming */\
\'a0 for ( i = 0; i < local_n_after; i++ ) \{\
\'a0 \'a0 c_re(u222[i]) /= (double)N;\
\'a0 \'a0 c_im(u222[i]) /= (double)N;\
\'a0 \}\
\
\'a0 /* Finally, compute the inverse transform;*/\
\'a0 fftw_mpi ( plan2, 1, u222, A);\
\
\'a0 /*print retransformed array */\
\'a0 printf("$$$$$$$$$$$$$$$444After inverse\
transformation:$$$$$$$$$$$$$$$$$$$$4\\n");\
\'a0 for ( j = 0; j < local_n_after; j++) \{\
\'a0 \'a0 printf("%i (on %i): (%f, %f)\\n", j+local_start, myid,\
\'a0 \'a0 \'a0 c_re(u222[j]), c_im(u222[j]));\
\'a0 \}\
\
\
\
\'a0 /* clean up!*/\
\'a0 fftw_mpi_destroy_plan ( plan1 );\
\'a0 fftw_mpi_destroy_plan ( plan2 );\
\
\'a0 if ( myid == 0 ) \{\
\'a0 \'a0 endwtime = MPI_Wtime();\
\'a0 \'a0 fprintf(stdout, "\\nWall clock time = %.6f seconds\\n", endwtime-startwtime);\
\'a0 \}\
\
\'a0 MPI_Finalize ();\
\
\'a0 return 0;\
\
\}}
/*
AUTHOR: ADRIAN DOBRA (ADOBRA@STAT.WASHINGTON.EDU)
*/

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "util.h"

#ifndef DBL_MAX
#define DBL_MAX 99999999999999999.0
#endif


//LAPACK FORTRAN FUNCTIONS

extern "C"
{
  void dgelsd_(int*,int*,int*,double*,int*,double*,int*,double*,double*,int*,double*,int*,int*,int*);  
  void dgetrf_(int*,int*,double*,int*,int*,int*);
  void dgetri_(int*,double*,int*,int*,double*,int*,int*);
  void dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
}

//some constants//////////////////
#define EPSILON 0.000001
//////////////////////////////////
const int randomseed = 16384;
//computes n choose m
int choose(int m, int n)
{
   if(m > n) return 0;
   if(m == 0) return 1;
   if(m == 1) return n;
   if(m == n) return 1;
   if(m == (n-1)) return n;

   return(choose(m, n-1) + choose(m-1, n-1));
}

///////////////////////////////
int IsValid(double s)
{
   if(fabs(s)<DBL_MAX)
      return 1;
   return 0;
}

//adjusts numbers that are too small to be significant
double AdjustDouble(double a)
{
   double b = fabs(a);
   if((b>EPSILON)&&(b<DBL_MAX)) return a;
   return 0.0;
}


void CheckPointer(void* pointer)
{
   if(pointer==NULL)
   {
      printf("Memory allocation error.\n");
      exit(1);
   }
   return;
}

int Solve(int M, double* A)
{
	int n = M;
	int ipiv[M];
	double work[5*M];
	int LWORK=5*M;
	int info=0;
	dgetrf_(&n,&n,&A[0],&n,&ipiv[0],&info);
	if(info!=0) return 0;
	dgetri_(&n,&A[0],&n,&ipiv[0],&work[0],&LWORK,&info);	
	if(info!=0) return 0;
	return 1;
}

//inverts a MxM matrix A
/*
int Solve(int M, double** A)
{
   int i,j;
   int n = M;
   double X[M][M];
   int ipiv[M];
   double work[5*M];
   int LWORK=5*M;
   int info=0;
 
   for(i=0;i<n;i++)
   {
      for(j=0;j<n;j++)
      {
	 X[i][j] = A[i][j];
      }
   }

   dgetrf_(&n,&n,&X[0][0],&n,&ipiv[0],&info);
   if(info!=0) return(0);
   dgetri_(&n,&X[0][0],&n,&ipiv[0],&work[0],&LWORK,&info);
   if(info!=0) return(0);

   for(i=0;i<n;i++)
   {
      for(j=0;j<n;j++)
      {
	 A[i][j] = X[i][j];
      }
   }

   return(1);
}
*/
double getlogdet(int nsize, double* A) {
	
	int N = nsize;
	char jobz = 'N';
	char uplo = 'L';
	int lda = nsize;
	double* eigvals = new double[N];
	double* work = new double[5*N];
	int LWORK=5*N;
	int info;

	double* X = new double[nsize * nsize];
	memcpy(X, A, nsize * nsize * sizeof(double));
 
	dsyev_(&jobz,&uplo,&N, X, &lda,eigvals, work, &LWORK,&info);
	double logdet = 0.0;
	if(0==info) {
		for(int i=0;i<N;i++) {
			logdet += log(fabs(eigvals[i]));
		}
	}
	delete[] eigvals;
	delete[] work;
	delete[] X;
	return logdet;

}
/*
double getlogdet(int nsize,double** A)
{
   int i,j;
   double logdet = 0.0;
   int N = nsize;
   char jobz = 'N';
   char uplo = 'L';
   int lda = N;
   double X[N][N];
   double eigvals[N];
   double work[5*N];
   int LWORK=5*N;
   int info;
 
   for(i=0;i<N;i++)
   {
      for(j=0;j<N;j++)
      {
         X[i][j] = A[i][j];
      }
   }

   dsyev_(&jobz,&uplo,&N,&X[0][0],&lda,&eigvals[0],&work[0],
          &LWORK,&info);

   if(0==info)
   {
      for(i=0;i<N;i++)
      {
         logdet += log(fabs(eigvals[i]));
      }
   }

   return(logdet);
}
*/






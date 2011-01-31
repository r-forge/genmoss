#include "f2c.h"
#include "clapack.h"
#include "blaswrap.h"

void copymatrix(int n,int p,double** source,double** dest)
{
	int i,j;
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			dest[i][j] = source[i][j];
		}
	}
	return;
}

double ** allocmatrix(int n,int p)
{
	int i;
	double** m;
	
	m = new double*[n];
	for(i=0;i<n;i++)
	{
		m[i] = new double[p];
		memset(m[i],0,p*sizeof(double));
	}
	return(m);
}

void freematrix(int n,double** m)
{
	int i;
	
	for(i=0;i<n;i++)
	{
		delete[] m[i]; m[i] = NULL;
	}
	delete[] m; m = NULL;
}

void meanmatrix(int n,int p,double** m,double** v)
{
	int i,j;
	double s;
	
	for(j=0;j<p;j++)
	{
		s = 0.0;
		for(i=0;i<n;i++)
		{
			s += m[i][j];
		}
		v[0][j] = s/((double)n);
	}
	return;
}

void samplemean(int n,int p,double** m,double* v)
{
	int i,j;
	double s;
	
	for(j=0;j<p;j++)
	{
		s = 0.0;
		for(i=0;i<n;i++)
		{
			s += m[i][j];
		}
		v[j] = s/((double)n);
	}
	return;
}

void matrixproduct(int n,int p,int l,double** m1,double** m2,double** m)
{
	int i,j,k;
	double s;
	
	for(i=0;i<n;i++)
	{
		for(k=0;k<l;k++)
		{
			s = 0;
			for(j=0;j<p;j++)
			{
				s += m1[i][j]*m2[j][k];
			}
			m[i][k] = s;
		}
	}
	return;
}

void readmatrix(char* filename,int n,int p,double* m[])
{
	int i,j;
	double s;
	FILE* in = fopen(filename,"r");
	
	if(NULL==in)
	{
		printf("Cannot open input file [%s]\n",filename);
		exit(1);
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
		{
			fscanf(in,"%lf",&s);
			m[i][j] = s;
		}
	}
	fclose(in);
	return;
}

void printmatrix(char* filename,int n,int p,double** m)
{
	int i,j;
	double s;
	FILE* out = fopen(filename,"w");
	
	if(NULL==out)
	{
		printf("Cannot open output file [%s]\n",filename);
		exit(1);
	}
	for(i=0;i<n;i++)
	{
		fprintf(out,"%.25lf",m[i][0]);
		for(j=1;j<p;j++)
		{
			fprintf(out,"\t%.10lf",m[i][j]);
		}
		fprintf(out,"\n");
	}
	fclose(out);
	return;
}

//input matrix of n samples and p variables
//output: sample covariance matrix
void samplecov(int n,int p,double** d,double** m)
{
	int i,j;
	double** X = allocmatrix(n,p);
	double** tX = allocmatrix(p,n);
	double** dmean = allocmatrix(1,p);
	double** ones = allocmatrix(n,1);
	double** onesmean = allocmatrix(n,p);

	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
		{
			X[i][j] = d[i][j];
		}
	}
	for(i=0;i<n;i++) ones[i][0] = 1;
	meanmatrix(n,p,X,dmean);
	matrixproduct(n,1,p,ones,dmean,onesmean);
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
		{
			X[i][j] -= onesmean[i][j];
			tX[j][i] = X[i][j];
		}
	}

	matrixproduct(p,n,p,(double**)tX,(double**)X,m);
	for(i=0;i<p;i++)
	{
		for(j=0;j<p;j++)
		{
			m[i][j] /= ((double)(n-1));
		}
	}

	//clean memory
	freematrix(1,dmean);
	freematrix(n,ones);
	freematrix(n,X);
	freematrix(p,tX);
	freematrix(n,onesmean);
	return;
}

//computes the Cholesky decomposition of
//a positive definite matrix
void chol(int p,double** m)
{
	int i,j;
	char uplo = 'L';
	double a[p][p];
	integer nn = (integer)p;
	integer lda = (integer)p;
	integer info;
	for(i=0;i<p;i++)
	{
		for(j=0;j<p;j++)
		{
			a[i][j] = m[i][j];
		}
	}
	dpotrf_(&uplo,&nn,(double*)a,&lda,&info);
	if(0!=info)
	{
		printf("Smth wrong in the call of 'dpotrf' error is [info = %d]\n",info);
		exit(1);
	}

	for(i=0;i<p;i++)
	{
		for(j=0;j<p;j++)
		{
			if(j<=i)
			{
				m[i][j] = a[j][i];
			}
			else
			{
				m[i][j] = a[i][j];
			}
		}
	}
}

//computes the inverse of a matrix
void inverse(int p,double** m)
{
	int i,j;
	char uplo = 'L';
	double a[p][p];
	integer nn = (integer)p;
	integer lda = (integer)p;
	integer info;
	for(i=0;i<p;i++)
	{
		for(j=0;j<p;j++)
		{
			a[i][j] = m[i][j];
		}
	}
	dpotri_(&uplo,&nn,(double*)a,&lda,&info);
	if(0!=info)
	{
		printf("Smth wrong in the call of 'dpotri' error is [info = %d]\n",info);
		exit(1);
	}

	for(i=0;i<p;i++)
	{
		for(j=0;j<p;j++)
		{
			if(j<=i)
			{
				m[i][j] = a[j][i];
			}
			else
			{
				m[i][j] = a[i][j];
			}
			printf("%d\t%d\t%.5lf\n",i,j,m[i][j]);
		}
	}
	exit(1);
	return;
}

double logdet(int p,double** m)
{
	int i,j;
	char jobvl = 'N';
	char jobvr = 'N';
	integer lda = (integer)p;
	double wr[2*p];
	double wi[2*p];
	double vl[p][p];
	integer ldvl = (integer)(p*p);
	double vr[p][p];
	integer ldvr = (integer)(p*p);
	double work[p*p];
	integer lwork = (integer)(p*p);
	double a[p][p];
	integer info;
	integer pp = (integer)p;

	for(i=0;i<p;i++)
	{
		for(j=0;j<p;j++)
		{
			a[i][j] = m[i][j];
		}
	}
	dgeev_(&jobvl,&jobvr,&pp,(double*)a,&lda,(double*)wr,(double*)wi,(double*)vl, 
		   &ldvl,(double*)vr,&ldvr,(double*)work,&lwork,&info);
	if(0!=info)
	{
		printf("Smth wrong in the call of 'dgeev' error is [info = %d]\n",info);
		exit(1);
	}	   
	
	double logdet = 0;
	for(i=0;i<p;i++) logdet+=log(wr[i]);	
	return(logdet);
}

double mvnpdf(double* r,int nParams,double thetaCovLogdet,double* thetaMean,double** thetaInverse)
{
	int i;
	double s = -(nParams*0.9189385332046726695409689)-(thetaCovLogdet/2.0);
	double** dmeant = allocmatrix(1,nParams);
	double** dmean = allocmatrix(nParams,1);
	double** M1 = allocmatrix(nParams,1);
	double** M2 = allocmatrix(1,1);
	
	for(i=0;i<nParams;i++)
	{
		dmeant[0][i] = r[i]-thetaMean[i];
		dmean[i][0] = dmeant[0][i];
	}
	matrixproduct(nParams,nParams,1,thetaInverse,dmean,M1);
	matrixproduct(1,nParams,1,dmeant,M1,M2);
	
	s -= M2[0][0]/2.0;
	
	//clean memory
	freematrix(1,dmeant);
	freematrix(nParams,dmean);
	freematrix(nParams,M1);
	freematrix(1,M2);
	return s;
}

/*
AUTHOR: ADRIAN DOBRA (ADOBRA@STAT.DUKE.EDU)
		QUANLI WANG (QUANLI@STAT.DUKE.EDU)
*/

// Data.cpp: implementation of the CData class.
//
//////////////////////////////////////////////////////////////////////

#include "Data.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CData::CData()
{
	NumberOfGenes = 0; 
	SampleSize = 0; 
	DataFile = "";

	Y = NULL;
	X = NULL;
	data = NULL;
}

CData::~CData()
{
	Cleanup();
}

void CData::ReadData()
{
	int i,j;
	double s;
   
	data = new double*[SampleSize];
	for(i=0;i<SampleSize;i++) {
		data[i] = new double[NumberOfGenes];
	}
	FILE* in = fopen(DataFile.c_str(),"r");
	if(NULL==in) {
		printf("Cannot open data file %s.\n",DataFile.c_str());
		exit(1);
	}
	for(i=0;i<SampleSize;i++) {
		for(j=0;j<NumberOfGenes;j++) {
			fscanf(in,"%lf",&s);
			data[i][j] = s;
		}
	}
	fclose(in);
}

void CData::WriteData(int* goodvars)
{
	int i,j;
	
	FILE* out = fopen(DataFile.c_str(),"w");
	if(NULL==out)
	{
		printf("Cannot open data file %s.\n",DataFile.c_str());
		exit(1);
	}
	for(i=0;i<SampleSize;i++)
	{
		for(j=0;j<NumberOfGenes;j++)
		{
			if(goodvars[j]) break;
		}
		fprintf(out,"%.5lf",data[i][j]);
		j++;
		for(;j<NumberOfGenes;j++)
		{
			if(goodvars[j])
			{
				fprintf(out,"\t%.5lf",data[i][j]);
			}
		}
		fprintf(out,"\n");
	}
	fclose(out);

	
	return;
}

void CData::Allocate(int ncolsx)
{
	ncolsX = ncolsx;
	///////////////////////////////
	//alloc the Y and X variables//
	///////////////////////////////
	Y = new double[SampleSize];
	X = new double*[SampleSize];
	// set the rest of the X matrix
	for(int i=0;i<SampleSize; i++) {
		X[i] = new double[ncolsX];
		memset(X[i],0,ncolsX*sizeof(double));
	}
}

void CData::FreeData()
{
	if (data != NULL) {
		int i;
		for(i=0;i<SampleSize;i++) {
			delete[] data[i]; 
			data[i] = NULL;
		}
		delete[] data;
		data = NULL;
	}
}	

void CData::Cleanup()
{
	//clean memory
	FreeData();
	if (X != NULL) {
		for(int i=0;i<SampleSize;i++) {  
			delete[] X[i]; X[i] = NULL;
		}
		delete[] X; X = NULL;
	}
	if (Y != NULL) {
		delete[] Y; Y = NULL;
	}
}

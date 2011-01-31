/*
AUTHOR: ADRIAN DOBRA (ADOBRA@STAT.WASHINGTON.EDU)
*/

#if !defined(CMODEL_H)
#define CMODEL_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

//#include "mkl_vsl.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class CModel
{
public:
	CModel();
	virtual ~CModel();
	
	double mdCutoffMax;
	double mdCutoffMin;
	void SetCutoffs(double vmax,double vmin);
	
	int notStudied;
	int Counter;
	int Length;
	int* vars;
	int* generators;
	double laplace;
	CModel* Next;
	
	void AddLaplace(int len,int* amodel,int* gen,double alaplace);
	CModel* Find(int len,int* amodel);
	void DeleteAll();
	void DeleteNotRelevant();
	void SaveAll(FILE* file);
	void SaveBestModel(FILE* file);
	void NormalizeWeights();
	void MakeCovProbs(int ncovprob,double* covprob); //to be always called after NormalizeWeights()
	
	double* NormalizeWeightsShotgun(double* w,int nmax);
	int WeightedSamplingShotgun(int maxk,double* weights, gsl_rng * stream);
	CModel* Shotgun(gsl_rng * stream,FILE* logfile);
	
	CModel* getBestModel();
	
	void Diagnostics(double& maxprob,int& nrelmodels);
	void SaveRelevant(FILE* file);
};

//returns -1 if m1 > m2, 0 if m1 == m2, 1 if m1 < m2
int LexCompare(int len,int* m1,int* m2);
#endif

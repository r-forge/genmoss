/*
AUTHOR: ADRIAN DOBRA (ADOBRA@STAT.WASHINGTON.EDU)
*/

#if !defined(REGMODELS_H)
#define REGMODELS_H

//#include "mkl_vsl.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <set>
using namespace std;

class CRegression
{
public:
	CRegression();
	virtual ~CRegression();
	
	int notStudied;
	double mdCutoffMax;
	double mdCutoffMin;
	int nIterationNumber;
	
	double Weight;
	double normWeight;
	set<int> Vars;
	CRegression* Next;
	int nmodels;
	
	void SetCutoffs(double vmax,double vmin);
	
	void ResetStudied();
	double getMaxProb();
	void AddList(CRegression* reglist);
	int NumberRegressions();
	int existReg(double weight,set<int> vars);
	
	void DeleteAll();
	void Delete();
	void DeleteNotRelevant(double maxweight);
	int AddReg(double weight,set<int> vars);
	double NormalizeWeights();
	void MakeCovProbs(int ncov,double* probs);
	void LoadAll(int nMaxRegressors,FILE* file);
	void SaveAll(int nMaxRegressors,FILE* file);
	
	double* NormalizeWeightsShotgun(double* w,int nmax);
	int WeightedSamplingShotgun(int maxk,double* weights, gsl_rng * stream);
	CRegression* Shotgun(gsl_rng * stream,int& nPmodels,int& nPtotalmodels);
	CRegression* FirstModel(gsl_rng * stream,int& nPmodels,int& nPtotalmodels);
};
#endif

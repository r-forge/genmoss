/*
AUTHOR: ADRIAN DOBRA (ADOBRA@STAT.WASHINGTON.EDU)
*/
// gibbssampler.h: interface for the CGibbsSampler class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(GIBBSSAMPLER_H)
#define GIBBSSAMPLER_H

//#include "mkl_vsl.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#include "f2c_.h"
#include "clapack.h"
#include "blaswrap.h"

//#include "cblas.h"


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <set>
#include <vector>
using namespace std;

#if !defined(SEED)
#define SEED 7777777
#endif

#if !defined(BRNG)
#define BRNG    gsl_rng_mt19937
#endif

#if !defined(DBL_MAX)
#define DBL_MAX 99999999999999999.0
#endif

#include "Data.h"

#ifndef _TABLE_H
#include "newtable.h"
#endif

#include "equations.h"
#include "util.h"
#include "regmodels.h"

class CGibbsSampler  
{
public:
	CGibbsSampler();
	virtual ~CGibbsSampler();

	//VSLStreamStatePtr stream;
	gsl_rng * stream;
	double tauPrior;

	//GIBBS sampling parameters
	double myPI;
	int nMaxRegressors; //this gives the maximum number of accepted regressors
	//int nMaxRegressorsInteraction;
	//int nDoInteractions;
	//int nDoSpaceRatio;
	//int nDoMetropolisH;
	int nChainIterations;
	//int nChainIterationsBurnin;
	int nChainReplicates;
	double mdCutoffMax;
	double mdCutoffMin;
	double mdProbMax;

	vector<int> mGeneList;
	void Cleanup();

	double calculateLogPost(int TargetGene, set<int>& ExplanatoryGenes, CData& Data);
	int ChooseNextGene(vector<int>& PossiblePredictors);
	int ChooseNextGeneSet(set<int>& PossiblePredictors);
	set<int> RandomRegression(set<int> Predictors);

	double metropolisReg(int targetGene,CData& Data,CRegression* hybridreglist);
		
	void exhaustiveReg(int targetGene,CData& Data,double* probs,FILE* modelsfile);

	CRegression* modelSelection(int targetGene,CData& Data,double* probs,FILE* modelsfile,FILE* logfile,int* countmodels);
	
	CRegression* modelSelectionMH(int targetGene,CData& Data,double* probs,FILE* modelsfile,FILE* logfile,int* countmodels);
	
	void modelSelectionInteraction(int targetGene,CData& Data,double* probs,FILE* modelsfile,FILE* logfile,CRegression* allregs);
	
	int shotgunReg(int targetGene,CData& Data,CRegression* reglist,int astartpoint,FILE* logfile);
	
	void mhReg(int targetGene,CData& Data,CRegression* reglist,int astartpoint,FILE* logfile,int nMhIterations);
	
	double interactReg(int targetGene,set<int>& MainEffects,CData& Data,CRegression* reglist);
	
	double* NormalizeWeights(double* w,int nmax);
	int WeightedSampling(int maxk,double* weights);
	int SelectEquation(int NumberOfGenes,CEquations* eq);
	
	void FirstIndex(int n,int* index,set<int>& VarInModel);
	int GetNextIndex(int n,int* index,set<int>& VarInModel);
	
	int MakeInteractionIndex(int p,int i,int j);
	void GetInteractionIndex(int ind,int p,int& vi,int& vj);
	
	double logitnr(int p, CData& Data);
	double loglikprior(double** beta,int p,CData& Data);
	double** gvector(double** beta,int p,CData& Data);
	double** gmatrix(double** beta,int p,CData& Data);
};

#endif 

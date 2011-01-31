#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iomanip>
#include <time.h>
#include <iostream>
#include <fstream>
using namespace std;

// includes for forking
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>

//#include "mkl_vsl.h"
//#include "mkl_lapack.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define SEED    7777777
#define BRNG    gsl_rng_mt19937
//#define METHOD  0

#ifndef PARAMS_H
#include "Params.h"
#endif

#ifndef _GRAPH_H
#include "graph.h"
#endif

#ifndef _TABLE_H
#include "newtable.h"
#endif

#ifndef DATA_H
#include "Data.h"
#endif

#ifndef MODEL_H
#include "model.h"
#endif

#include<R.h>
#include<Rmath.h>
#include<Rdefines.h>

/**** INPUT PARAMETERS ****/
double cPrior = -1; //prior for the table

/**** OUTPUT PARAMETERS ****/

/**** GLOBAL VARIABLES ****/
Params params;
//VSLStreamStatePtr stream;
gsl_rng * stream;
int n = -1; //number of dimensions in the table
CData TrainDiscrete;
CData TestDiscrete;

/*
int** VarSets = NULL;
int* lenVarSets = NULL;
int nVarSets = -1;
*/

int** DownLinks = NULL;
int* nDownLinks = NULL;
int** UpLinks = NULL;
int* nUpLinks = NULL;

double* compPost = NULL;
double* modelProb = NULL;
double* modelStd = NULL;
/**** GLOBAL VARIABLES ****/

/**** FUNCTIONS ****/
void InitVarSets(int n,int**& VarSets,int*& lenVarSets,int& nVarSets);
void DeleteVarSets(int**& VarSets,int*& lenVarSets,int& nVarSets);
int subsetVarSets(int n,int set1,int set2); //returns 1 if set1 is a subset of set2
int subset(int n,int* set1,int* set2);
extern "C" void rtraintest(char**);

/**** FUNCTIONS ****/



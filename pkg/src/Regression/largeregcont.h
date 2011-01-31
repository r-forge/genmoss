/*
 *  largeregcont.h
 *  
 *
 *  Created by Adrian Dobra on 11/22/06.
 *  Copyright 2006 __University of Washington__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <iomanip>
#include <time.h>

// necessary for forking:
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>

using namespace std;

//#include "mkl_vsl.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "f2c_.h"
#include "clapack.h"
#include "blaswrap.h"

#define SEED    777777
#define BRNG    gsl_rng_mt19937

#ifndef _GRAPH_H
#include "graph.h"
#endif

#ifndef _TABLE_H
#include "newtable.h"
#endif

#include "Model.h"
#include "util.h"
#include "Data.h"
#include "gibbssampler.h"

#include<R.h>
#include<Rmath.h>
#include<Rdefines.h>



/**** GLOBAL VARIABLES ****/
Model model;
//VSLStreamStatePtr mystream;
gsl_rng * mystream;
CData Data;
CGibbsSampler GS;
double* probs = NULL;

double rndModelBern = 0.3; 
//that's the probability of inclusion of an extra model term 
//in a random model used as the starting point for the sampler

double M = 0.01;
/**** GLOBAL VARIABLES ****/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <iomanip>
#include <time.h>
#include <iostream>
#include <fstream>

// necessary for forking:
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>


using namespace std;

//#include "mkl_vsl.h"
#include <gsl/gsl_rng.h>

#define SEED    777777
//#define BRNG    VSL_BRNG_MT19937
#define BRNG gsl_rng_mt19937
#define METHOD  0

#include "Model.h"

#ifndef _TABLE_H
#include "newtable.h"
#endif

#include "Data.h"



#include<R.h>
#include<Rmath.h>
#include<Rdefines.h>



/**** GLOBAL VARIABLES ****/
Model model;
//VSLStreamStatePtr mystream;
gsl_rng * mystream;
CData TrainData;
CData TrainDiscrete;
CData TestData;
CData TestDiscrete;
/**** GLOBAL VARIABLES ****/

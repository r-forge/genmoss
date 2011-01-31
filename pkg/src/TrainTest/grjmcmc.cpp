#include "grjmcmc.h"
#include "lattice.functions.cpp"
#include "varsets.functions.cpp"
#include "posterior.functions.cpp"
#include "process.models.cpp"

void traintest_main(char* fname)
{
	if (fname == NULL)
	{
		cout << "USAGE: [parametersfile]" << endl;
		exit(1);
	}
	if (!params.Load(fname)) {
	  cout << params.GetErrorMessage() << endl;
	  exit(1);
	}

    //int i, j;//, errcode;
 
    /***** Initialize the RND library*****/
    //errcode = vslNewStream( &stream,BRNG,SEED);
 	stream = gsl_rng_alloc(BRNG);
	gsl_rng_set(stream, SEED);

    //read the train data
	TrainDiscrete.NumberOfGenes = params.mnNumberOfGenes;
	TrainDiscrete.SampleSize = params.mnTrainSampleSize;   
	TrainDiscrete.DataFile = params.mstrTrainDataFile;
    TrainDiscrete.ReadData();
	
	//read the test data
	TestDiscrete.NumberOfGenes = params.mnNumberOfGenes;
	TestDiscrete.SampleSize = params.mnTestSampleSize;   
	TestDiscrete.DataFile = params.mstrTestDataFile;
    TestDiscrete.ReadData();
	
	n  = params.mnNumberOfGenes;
   		
	doRJMCMC(n,TrainDiscrete,TestDiscrete,cPrior);

    //clear data structures
	TrainDiscrete.Cleanup();
	TestDiscrete.Cleanup();
	
    /***** Deinitialize the RND library*****/
    //errcode = vslDeleteStream(&stream);
	gsl_rng_free(stream);
}

extern "C" void rtraintest(char **fname)
{
        pid_t pID = fork();

        if (pID == 0) {                 // child
		char *fname2 = fname[0];
		traintest_main(fname2);
        	exit(0);
        }
        else if (pID <0) {
                printf("Failed to fork to process the train/test request.\n");
        }
        else {                          // parent
                int childstatus;
                waitpid(pID, &childstatus, 0);
        }

}


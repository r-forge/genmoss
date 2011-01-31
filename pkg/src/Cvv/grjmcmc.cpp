#include "grjmcmc.h"
#include "lattice.functions.cpp"
#include "varsets.functions.cpp"
#include "posterior.functions.cpp"
#include "process.models.cpp"

void cvv_main(char* fname)
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

   // int i, j;//, errcode;
 
    /***** Initialize the RND library*****/
    //errcode = vslNewStream( &stream,BRNG,SEED);
	stream = gsl_rng_alloc(BRNG);
	gsl_rng_set(stream, SEED); 

    //read the train data
	Discrete.NumberOfGenes = params.mnNumberOfGenes;
	Discrete.SampleSize = params.mnSampleSize;   
	Discrete.DataFile = params.mstrDataFile;
    Discrete.ReadData();
	
	n  = params.mnNumberOfGenes;
   		
	doRJMCMC(n,Discrete,cPrior);

    //clear data structures
	Discrete.Cleanup();
	
    /***** Deinitialize the RND library*****/
    //errcode = vslDeleteStream(&stream);
	gsl_rng_free(stream);
}

extern "C" void rcvv(char **fname)
{
        pid_t pID = fork();

        if (pID == 0) {                 // child
		char * fname2 = fname[0];
		cvv_main(fname2);
                exit(0);
        }
        else if (pID <0) {
                printf("Failed to fork to process the cvv request.\n");
        }
        else {                          // parent
                int childstatus;
                waitpid(pID, &childstatus, 0);
        }

}

/*
 *  largeregcont.cpp
 *  
 *
 *  Created by Adrian Dobra on 11/22/06.
 *  Copyright 2006 __University of Washington__. All rights reserved.
 *
 */

#include "largeregcont.h"
#include "cmodel.h"

#include "lattice.functions.cpp"
#include "varsets.functions.cpp"
#include "posterior.functions.cpp"

int qsortcompare(const void * a, const void * b)
{
	int ia = *(int*)a;
	int ib = *(int*)b;
	
	if(probs[ia]>probs[ib]) return(-1);
	if(probs[ia]<probs[ib]) return(1);
	return(0);
}

void TableShotgunSearch(FILE* out,LPTable dataTable,
					    int ShotgunChainReplicates,double ShotgunCutoffMax,double ShotgunCutoffMin,double ShotgunProbMax, int nconfs)
{
	int i;
	
	if(dataTable->nDimens<=(2+nconfs))
	{
		fprintf(out,"0");
		for(i=1;i<dataTable->Total;i++)
		{
			fprintf(out," 1");
		}
		fprintf(out," 1.0\n");
		return;
	}
	
	int** VarSets = NULL;
	int* lenVarSets = NULL;
	int nVarSets = -1;

	int** DownLinks = NULL;
	int* nDownLinks = NULL;
	int** UpLinks = NULL;
	int* nUpLinks = NULL;
	
	LPTable priorTable = new Table;
	InitPriorTable(priorTable,dataTable,-1);

	InitVarSets(dataTable->nDimens,
	            VarSets,lenVarSets,nVarSets);
	InitLattice(dataTable->nDimens,
				VarSets,lenVarSets,nVarSets,
				DownLinks,nDownLinks, 
				UpLinks,nUpLinks);

	LPTable datapriorTable = new Table;
	if(!datapriorTable->Alloc(dataTable->Dimens,dataTable->nDimens))
	{
		printf("Failed to allocate data+prior table.\n"); exit(1);
	}
	for(i=0;i<datapriorTable->Total;i++)
	{
		datapriorTable->Data[i] = dataTable->Data[i]+priorTable->Data[i];
	}

	//this is where we store all the models we identify
	CModel* models = new CModel;
	models->SetCutoffs(ShotgunCutoffMax,ShotgunCutoffMax);
	
	for(int astartpoint=1;astartpoint<=ShotgunChainReplicates;astartpoint++)
	{
		CModel* localmodels = new CModel;
		localmodels->SetCutoffs(ShotgunCutoffMax,ShotgunCutoffMin);
		
		doRJMCMCstart(models,localmodels,astartpoint,
		              dataTable,priorTable,datapriorTable,
		              VarSets,lenVarSets,nVarSets,
				      DownLinks,nDownLinks,UpLinks,nUpLinks,
					  mystream);
						
		localmodels->DeleteAll();
		delete localmodels; localmodels = NULL;
	}

	//save the best model identified
	models->NormalizeWeights();
	models->SaveBestModel(out);

	//clean memory
	models->DeleteAll();
	delete models; models = NULL;
	datapriorTable->Reset(); delete datapriorTable; datapriorTable = NULL;
	priorTable->Reset(); delete priorTable; priorTable = NULL;
	DeleteLattice(VarSets,lenVarSets,nVarSets,
				  DownLinks,nDownLinks, 
				  UpLinks,nUpLinks);
	DeleteVarSets(VarSets,lenVarSets,nVarSets);			  
	return;
}

void FindHierarchicalModels(CRegression* allregs,CData& Data,char* DataFileName,int target,int nMaxRegressors)
{
	int i,j;
	char buffer[2048];
	FILE* out = NULL;
	int responseVariable = target-1;
	int lenmodel = 1+nMaxRegressors + Data.NumOfConfoundingVars; // + all confounding vars
	int* amodel = new int[lenmodel];
	CRegression* p = allregs->Next;
	
	int modelid = 1;
	while(NULL!=p)
	{
		set<int>::iterator it;
		lenmodel = 0;
		for(it=(p->Vars).begin();it!=(p->Vars).end();it++)
		{
			amodel[lenmodel] = *it;
			lenmodel++;
		}

		// Next insert all the confounding vars, followed by the response var
	        j = Data.NumOfConfoundingVars;
        	for(i=0; i<=j; i++)
        	{
                	amodel[lenmodel] = responseVariable-j+i;
			lenmodel++;
        	}

		//amodel[lenmodel] = responseVariable;
		//lenmodel++;
		
		printf("Model [%d] ::",modelid);
		for(i=0;i<lenmodel;i++) printf(" %d",amodel[i]);
		printf("\n");
		
		//get the marginal table for these variables
		LPTable dataTable = new Table;
		int* index = new int[lenmodel];
		for(j=0;j<lenmodel;j++)
		{
			index[j] = 2;
		}
		if(!dataTable->Alloc(index,lenmodel))
		{
			printf("Error allocating memory!\n");
			exit(1);
		}
		for(i=0;i<Data.SampleSize;i++)
		{
			for(j=0;j<lenmodel;j++)
			{
				index[j] = (int) (Data.data[i][amodel[j]]);
			}
			dataTable->SetIndex(index);
			dataTable->Set(dataTable->Get()+1);
		}
		
		//file where the best model will be saved
		sprintf(buffer,"%s.shotgun.%d.%d.reg.model%d.txt",DataFileName,target,nMaxRegressors,modelid);
		if(NULL==(out=fopen(buffer,"w"))) 
		{
			printf("Cannot open file [%s]\n",buffer);
			return;
		}
		
		//do the shotgun search
		TableShotgunSearch(out, dataTable, model.mnShotgunChainReplicates,
			model.mdShotgunCutoffMax, model.mdShotgunCutoffMin, model.mdShotgunProbMax, Data.NumOfConfoundingVars);
		fclose(out);
		
		//clean memory
		dataTable->Reset();
		delete dataTable; dataTable = NULL;
		
		//go to the next model
		p = p->Next;
		modelid++;
	}
	//clean memory
	delete[] amodel; amodel = NULL;
	return;
}

void searchRespReg(int target)
{
	double epsilon = 0.0001;
	int i;
	double s;
	char buffer[1024];
	FILE* out = NULL;
	FILE* logfile = NULL;
	//FILE* spaceratiofile = NULL;
	FILE* countmodelsfile = NULL;
	
	int* countmodels = new int[GS.nChainReplicates];
	memset(countmodels,0,GS.nChainReplicates*sizeof(int));
	
	//MAIN EFFECTS SELECTION BEGINS
	sprintf(buffer,"%s.shotgun.%d.%d.reg",Data.DataFile.c_str(),target,GS.nMaxRegressors);
	if(NULL==(out=fopen(buffer,"w"))) 
	{
		printf("Cannot open file [%s]\n",buffer);
		return;
	}
	sprintf(buffer,"%s.shotgun.%d.%d.log",Data.DataFile.c_str(),target,GS.nMaxRegressors);
	if(NULL==(logfile=fopen(buffer,"w"))) 
	{
		printf("Cannot open file [%s]\n",buffer);
		return;
	}
	// sprintf(buffer,"%s.shotgun.%d.%d.spaceratio.txt",Data.DataFile.c_str(),target,GS.nMaxRegressors);
	//if(NULL==(spaceratiofile=fopen(buffer,"w"))) 
	//{
	//	printf("Cannot open file [%s]\n",buffer);
	//	return;
	//}
	CRegression* allregs = GS.modelSelection(target-1,Data,probs,out,logfile,countmodels);
	fclose(out);
	fclose(logfile);
	//fclose(spaceratiofile);

	//take each regression and determine its best hierarchical log-linear model
	FindHierarchicalModels(allregs,Data,(char*)Data.DataFile.c_str(),target,GS.nMaxRegressors);

	sprintf(buffer,"%s.shotgun.%d.%d.countmodels.txt",Data.DataFile.c_str(),target,GS.nMaxRegressors);
	if(NULL==(countmodelsfile=fopen(buffer,"w"))) 
	{
		printf("Cannot open file [%s]\n",buffer);
		return;
	}
	for(i=0;i<GS.nChainReplicates;i++)
	{
		fprintf(countmodelsfile,"%d\n",countmodels[i]);
	}
	fclose(countmodelsfile);

	int index[Data.NumberOfGenes];
	for(i=0;i<Data.NumberOfGenes;i++)
	{
		index[i] = i;
	}
	qsort(index,Data.NumberOfGenes,sizeof(int),qsortcompare);
	
	sprintf(buffer,"%s.shotgun.%d.%d.var",Data.DataFile.c_str(),target,GS.nMaxRegressors);
	if(NULL==(out=fopen(buffer,"w")))
	{
		printf("Cannot open file [%s]\n",buffer);
		return;
	}
		
	for(i=0;i<Data.NumberOfGenes;i++)
	{
		if(index[i]!=(target-1))
		{
			if(probs[index[i]]>=epsilon)
			{
				//fprintf(out,"%d\t", target);
				// print the column ids of confounding vars:
				int fromm = Data.NumberOfGenes-Data.NumOfConfoundingVars;
				for(fromm; fromm<=target; fromm++)
					fprintf(out,"%d\t", fromm);
				fprintf(out,"%d\t%.5lf\n",index[i]+1,probs[index[i]]);
				// fprintf(out,"%d\t%d\t%.5lf\n",target,index[i]+1,probs[index[i]]);
			}
		}
	}
	fclose(out);
	//MAIN EFFECTS SELECTION ENDS
	allregs->DeleteAll();
	delete allregs; allregs = NULL;
	delete[] countmodels; countmodels = NULL;
	return;
}

void shotgun_main(char* fname)  
{
	if (fname == NULL)
	{
		cout << "USAGE: [modelfile]" << endl;
		exit(1);
	}
	if (!model.Load(fname)) {
	  cout << model.GetErrorMessage() << endl;
	  exit(1);
	}

	int i,j;
	int nMaxCandidates = (int)((model.mnMaxRegressors/2.0)*((double)(model.mnMaxRegressors+3)))+1+model.mnNumberOfGenes;
	
	probs = new double[nMaxCandidates];
	mystream = gsl_rng_alloc(BRNG);
	gsl_rng_set(mystream, SEED);

	Data.NumberOfGenes = model.mnNumberOfGenes;
	Data.SampleSize = model.mnSampleSize;   
	Data.DataFile = model.mstrDataFile;
	Data.NumOfConfoundingVars = model.mnNumOfConfoundingVars;
    //read the dataset
    Data.ReadData();
    Data.Allocate(nMaxCandidates);


	//simulation parameters
	GS.tauPrior = 1;
	GS.nChainIterations = model.mnChainIterations;
	//GS.nChainIterationsBurnin = model.mnChainIterationsBurnin;
	GS.nChainReplicates = model.mnChainReplicates;
	GS.nMaxRegressors = model.mnMaxRegressors;
	//GS.nMaxRegressorsInteraction = model.mnMaxRegressorsInteraction;
	//GS.nDoInteractions = model.mnDoInteractions;
	//GS.nDoSpaceRatio = model.mnDoSpaceRatio;
	//GS.nDoMetropolisH = model.mnDoMetropolisH;
	GS.mdCutoffMax = model.mdCutoffMax;
	GS.mdCutoffMin = model.mdCutoffMin;
	GS.mdProbMax = model.mdProbMax;

	//it is assumed that the response variable is the last
	//variable in the dataset
	searchRespReg(Data.NumberOfGenes);

    //clean memory
    Data.Cleanup();

	gsl_rng_free(mystream);

	delete[] probs; probs = NULL;
}

extern "C" void rshotgun(char **fname)
{
	pid_t pID = fork();

        if (pID == 0) {			// child
		char *fname2 = fname[0];
		shotgun_main(fname2);
		exit(0);
	}
	else if (pID <0) {
                printf("Failed to fork to process the request.\n");
        }
	else {				// parent
		int childstatus;
                pid_t ws = waitpid(pID, &childstatus, 0);
	}
}



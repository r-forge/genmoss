void LoadModels(CModel* models,int nMaxRegressors,int responseVariable)
{
	FILE* modelsfile = NULL;
	int i,j;
	double s,laplace;
	char buffer[2048];

	if(NULL==(modelsfile=fopen(params.mstrModelsFile.c_str(),"r")))
	{
		fprintf(stderr,"Cannot open file [%s]\n",params.mstrModelsFile.c_str());
		exit(1);
	}
	
	int lenmodel = 1+nMaxRegressors;
	int* amodel = new int[lenmodel];

	printf("lenmodel = %d :: responseVariable = %d\n",lenmodel,responseVariable);

	int modelid = 1;
	while(1)
	{
		if(1!=fscanf(modelsfile,"%lf",&laplace))
		{
			break;
		}
		fscanf(modelsfile,"%lf",&s);
		
		for(j=0;j<nMaxRegressors;j++)
		{
			fscanf(modelsfile,"%d",&i);
			amodel[j] = i-1;
		}
		amodel[nMaxRegressors] = responseVariable;
		
		sprintf(buffer,"%s.model%d.txt",params.mstrModelsFile.c_str(),modelid);
		string modelPath(buffer);
		models->AddLaplace(lenmodel,amodel,modelPath,laplace);
		modelid++;
	}
	fclose(modelsfile);
	delete[] amodel; amodel = NULL;
	return;
}

void WriteTrueValues(double* predictions,int ResponseVariable,LPTable dataTable)
{
	int i,j,k;
	
	dataTable->GetFirst();
	i = 0;
	k = 0;
	while(1)
	{
		for(j=0;j<dataTable->Data[k];j++)
		{
			predictions[i] = dataTable->Index[ResponseVariable-1];
			i++;
		}
		k++;
		if(!dataTable->GetNext())
		{
			break;
		}
	}
	
	return;
}

void SavePredictions(LPTable P,LPTable dataTable,int* vars,int firstsample,int lastsample,double** testdata,double* predictions)
{
	int i,j,k;
	//int nSampleSize = (int) dataTable->GetGrandTotal();
	int ResponseVariable = dataTable->nDimens-1;
		
	k = 0;
	for(i=firstsample;i<lastsample;i++)
	{
		for(j=0;j<dataTable->nDimens;j++)
		{
			P->Index[j] = (int) testdata[k][vars[j]];
		}
		k++;
		
		P->Index[ResponseVariable] = 1;
		double p1 = P->Get();
		P->Index[ResponseVariable] = 0;
		double p0 = P->Get();		
		predictions[i] = p1/(p0+p1);
	}
	
	return;
}

void getMarginal(LPTable dataTable,int lenmodel,int* model,int nSampleSize,double** Discrete)
{
   int i,j;

   int nvars = 0;
   for(i=0;i<lenmodel;i++)
   {
      if(model[i]>=0) nvars++;
   }
   
   int* index = new int[nvars];
   for(j=0;j<nvars;j++)
   {
      index[j] = 2;
   }
   if(!dataTable->Alloc(index,nvars))
   {
      printf("Error allocating memory!\n");
      exit(1);
   }
   int* vars = new int[nvars];
   nvars = 0;
   for(i=0;i<lenmodel;i++)
   {
      if(model[i]>=0)
	  {
	     vars[nvars] = model[i];
	     nvars++;
	  }
   }	
    
   for(i=0;i<nSampleSize;i++)
   {
      for(j=0;j<nvars;j++)
      {
	     index[j] = (int) Discrete[i][vars[j]];
      }
      dataTable->SetIndex(index);
      dataTable->Set(dataTable->Get()+1);
   }

   delete[] index; index = NULL;
   delete[] vars; vars = NULL;
   return;
}


// void SampleModelParameters(string modelPath,FILE* predfile,int niterations,int nvars,int* vars,CData& Discrete)
void SampleModelParameters(string modelPath,const char* predfilename,int niterations,int nvars,int* vars,CData& Discrete)
{
	int i,j,k;
	double s;
	
	const int nSampleSize = Discrete.SampleSize;
	int blocksize = (int)floor(((double)nSampleSize)/((double)params.mnCvvFold));
	int residualsize = nSampleSize-params.mnCvvFold*blocksize;
	
	printf("blocksize = %d :: residualsize = %d\n",blocksize,residualsize);
	
	const int NITER = niterations;
	double allpredictions[NITER][nSampleSize];
	
	int firstsample = 0;
	int lastsample = firstsample;
	for(int asplit=0;asplit<params.mnCvvFold;asplit++)
	{
		firstsample = lastsample;
		lastsample += blocksize;
		if(residualsize)
		{
			lastsample++;
			residualsize--;
		}
		if(lastsample>nSampleSize)
		{
			lastsample = nSampleSize;
		}
		printf("firstsample = %d :: lastsample = %d\n",firstsample,lastsample);
		
		//initialize the train and the test data
		int nTestSampleSize = lastsample-firstsample;
		int nTrainSampleSize = nSampleSize-nTestSampleSize;
		//printf("nTestSampleSize = %d :: nTrainSampleSize = %d",nTestSampleSize,nTrainSampleSize);
		
		double** testdata = new double*[nTestSampleSize];
		for(i=0;i<nTestSampleSize;i++)
		{
			testdata[i] = new double[Discrete.NumberOfGenes];
		}
		double** traindata = new double*[nTrainSampleSize];
		for(i=0;i<nTrainSampleSize;i++)
		{
			traindata[i] = new double[Discrete.NumberOfGenes];
		}
		
		k=0;
		for(i=firstsample;i<lastsample;i++)
		{
			for(j=0;j<Discrete.NumberOfGenes;j++)
			{
				testdata[k][j] = Discrete.data[i][j];
			}
			k++;
		}
		k=0;
		for(i=0;i<firstsample;i++)
		{
			for(j=0;j<Discrete.NumberOfGenes;j++)
			{
				traindata[k][j] = Discrete.data[i][j];
			}
			k++;
		}
		for(i=lastsample;i<nSampleSize;i++)
		{
			for(j=0;j<Discrete.NumberOfGenes;j++)
			{
				traindata[k][j] = Discrete.data[i][j];
			}
			k++;
		}
		
		//create the tables
		LPTable dataTable = new Table;
		getMarginal(dataTable,nvars,vars,nTrainSampleSize,traindata);
		LPTable mytable = new Table;
		if(!mytable->Alloc(dataTable->Dimens,dataTable->nDimens))
		{
			printf("Failed to allocate data+prior table.\n"); exit(1);
		}
		s = cPrior;
		if(s<0)
		{
			s = 1.0/((double)mytable->Total);
		}
		for(j=0;j<mytable->Total;j++)
		{
			mytable->Data[j] = dataTable->Data[j]+s;
		}
	
		LPTable Y = new Table;
		MakeMarginals(Y,mytable);
	
		LPTable Theta = new Table;
		if(!Theta->Alloc(mytable->Dimens,mytable->nDimens))
		{
			printf("Failed to allocate the theta's table.\n"); exit(1);
		}
		
		LPTable junkTheta = new Table;
		if(!junkTheta->Alloc(mytable->Dimens,mytable->nDimens))
		{
			printf("Failed to allocate the theta's table.\n"); exit(1);
		}
		
		LPTable P = new Table;
		if(!P->Alloc(mytable->Dimens,mytable->nDimens))
		{
			printf("Failed to allocate the probs' table.\n"); exit(1);
		}

		//load the best hierarchical model for this marginal
		int* localmodel = new int[mytable->Total];
		memset(localmodel,0,mytable->Total*sizeof(int));
		int* localmodelgenerators = new int[mytable->Total];
		memset(localmodelgenerators,0,mytable->Total*sizeof(int));
	
		FILE* modelfile = fopen(modelPath.c_str(),"r");
		if(NULL==modelfile)
		{
			printf("Cannot open model file [%s]\n",modelPath.c_str());
			exit(1);
		}
		for(i=0;i<mytable->Total;i++)
		{
			if(1!=fscanf(modelfile,"%d",&j))
			{
				printf("Model file [%s] is incorrect.\n",modelPath.c_str());
				exit(1);
			}
			localmodel[i] = j;
		}
		fclose(modelfile);
	
		int** VarSetsLocal = NULL;
		int* lenVarSetsLocal = NULL;
		int nVarSetsLocal = -1;
		InitVarSets(mytable->nDimens,VarSetsLocal,lenVarSetsLocal,nVarSetsLocal);
	
		int** DownLinks = NULL;
		int* nDownLinks = NULL;
		int** UpLinks = NULL;
		int* nUpLinks = NULL;
	
		InitLattice(mytable->nDimens,
					VarSetsLocal,lenVarSetsLocal,nVarSetsLocal,
					DownLinks,nDownLinks, 
					UpLinks,nUpLinks);
	
		//determine the generators of this model
		findGenerators(0,localmodel,nVarSetsLocal,UpLinks,nUpLinks,localmodelgenerators);
	
		//this is the mode of the distribution
		ipf(VarSetsLocal,lenVarSetsLocal,nVarSetsLocal,
			P,mytable,localmodelgenerators,localmodel);
		SavePredictions(P,dataTable,vars,firstsample,lastsample,testdata,allpredictions[0]);	
		ThetafromP(Theta,P,localmodel);	
		for(int iter=1;iter<NITER;iter++)
		{
			sampleNewThetas(VarSetsLocal,lenVarSetsLocal,nVarSetsLocal,
							localmodel,localmodelgenerators,mytable,Y,Theta,junkTheta);
			PfromTheta(P,Theta);
			SavePredictions(P,dataTable,vars,firstsample,lastsample,testdata,allpredictions[iter]);	
		}
		
		//clean memory
		DeleteLattice(VarSetsLocal,lenVarSetsLocal,nVarSetsLocal,
						DownLinks,nDownLinks, 
						UpLinks,nUpLinks);
		DeleteVarSets(VarSetsLocal,lenVarSetsLocal,nVarSetsLocal);
		delete[] localmodel; localmodel = NULL;
		delete[] localmodelgenerators; localmodelgenerators = NULL;
		Y->Reset(); delete Y; Y = NULL;
		P->Reset(); delete P; P = NULL;
		Theta->Reset(); delete Theta; Theta = NULL;
		junkTheta->Reset(); delete junkTheta; junkTheta = NULL;
		mytable->Reset(); delete mytable; mytable = NULL;					
		dataTable->Reset(); delete dataTable; dataTable = NULL;
		for(i=0;i<nTestSampleSize;i++)
		{
			delete[] testdata[i]; testdata[i] = NULL;
		}
		delete[] testdata; testdata = NULL;
		for(i=0;i<nTrainSampleSize;i++)
		{
			delete[] traindata[i]; traindata[i] = NULL;
		}
		delete[] traindata; traindata = NULL;
	}
	
	for(i=0;i<NITER;i++)
	{
		ofstream predfile;

		predfile.open(predfilename, ios::out | ios::app);
                
		if(!predfile.is_open()) {
			fprintf(stderr,"Cannot re-open file [%s]\n",predfilename);
			return;
		}
		predfile << allpredictions[i][0]; //fprintf(predfile,"%.3lf",allpredictions[i][0]);
		for(j=1;j<nSampleSize;j++)
		{
			predfile << "\t" << allpredictions[i][j]; //fprintf(predfile,"\t%.3lf",allpredictions[i][j]);
		}
		predfile << "\n"; //fprintf(predfile,"\n");
		predfile.close();
	}
	
	return;
}

void doRJMCMC(int n,CData& Discrete,double cPrior)
{
	int j;
	//double s;
	char predfilename[2048];
	//char buffer[2048];
	//FILE* predfile = NULL;
	double maxprob;
	int nrelmodels;
			
	//the response variable is assumed to be the last column in the data matrix		
	int responseVariable = Discrete.NumberOfGenes-1;			
						
	sprintf(predfilename,"%s.%d.cvv.txt",params.mstrDataFile.c_str(),responseVariable);

	ofstream predfile(predfilename);
	//if(NULL==(predfile=fopen(predfilename,"w")))
	if(!predfile.is_open())
	{
		fprintf(stderr,"Cannot open file [%s]\n",predfilename);
		return;
	}
					
	//this is where we store all the models we identify
	CModel* models = new CModel;
	models->SetCutoffs(0,0);
	LoadModels(models,params.mnMaxRegressors,responseVariable);
	models->NormalizeWeights();
	models->Diagnostics(maxprob,nrelmodels);
	printf("nRelevantModels = %d\n",nrelmodels);
	
	predfile << Discrete.data[0][responseVariable]; //fprintf(predfile,"%.3lf",Discrete.data[0][responseVariable]);
	for(j=1;j<Discrete.SampleSize;j++)
	{
		predfile << "\t" << Discrete.data[j][responseVariable]; //fprintf(predfile,"\t%.3lf",Discrete.data[j][responseVariable]);
	}
	predfile << "\n"; //fprintf(predfile,"\n");

	//params.mnChainIterations
	CModel* amodel = models->Next;
	int modelid = 1;
	int totaliterations = 0;
	predfile.close(); 
	while(NULL!=amodel)
	{
		int niterations = (int)rint(params.mnChainIterations*amodel->laplace);
		printf("Sampling model [%d] for [%d] iterations.\n",modelid,niterations);
		totaliterations += niterations;
		modelid++;
		if(niterations>=1)
		{			
			SampleModelParameters(amodel->modelPath,predfilename,niterations,amodel->Length,amodel->vars,Discrete);
		}
		amodel = amodel->Next;
	}
	printf("Total iterations = %d\n",totaliterations);
	
	models->DeleteAll();	
	//predfile.close(); //fclose(predfile);
	//clean memory
	delete models; models = NULL;
	return;
}

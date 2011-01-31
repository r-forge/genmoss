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

void SavePredictions(LPTable P,LPTable dataTable,const char* predfilename,int* vars,CData& Discrete)
{
	int i,j;
	const int nSampleSize = Discrete.SampleSize;
	double predictions[nSampleSize];
	int ResponseVariable = dataTable->nDimens-1;
	
	for(i=0;i<Discrete.SampleSize;i++)
	{
		for(j=0;j<dataTable->nDimens;j++)
		{
			//printf("vars[%d] = %d\n",j,vars[j]);
			P->Index[j] = (int) Discrete.data[i][vars[j]];
		}
		P->Index[ResponseVariable] = 1;
		double p1 = P->Get();
		P->Index[ResponseVariable] = 0;
		double p0 = P->Get();		
		predictions[i] = p1/(p0+p1);
	}

	ofstream predfile;
	predfile.open(predfilename, ios::out | ios::app);

	if(!predfile.is_open()) {
		fprintf(stderr,"Cannot re-open file [%s]\n",predfilename);
		return;
	}

	predfile << predictions[0]; //fprintf(predfile,"%.3lf",predictions[0]);
	for(j=1;j<nSampleSize;j++)
	{
		predfile << "\t" << predictions[j]; //fprintf(predfile,"\t%.3lf",predictions[j]);
	}
	predfile << "\n"; //fprintf(predfile,"\n");
	predfile.close();

	return;
}

//void SaveCoef(const char* filename,LPTable mytable)
void SaveCoef(LPTable mytable)
{
	int ResponseVariable = mytable->nDimens;
	
	mytable->GetFirst();
	int first = 1;

	/*ofstream predfile;
	predfile.open(filename, ios::out | ios::app);
	if(!predfile.is_open()) {
		fprintf(stderr,"Cannot re-open file [%s]\n",filename);
		return;
	}*/
	while(mytable->GetNext())
	{
		if(mytable->Index[ResponseVariable-1])
		{
			if(first)
			{
				first = 0;
				cout << mytable->Get(); //fprintf(file,"%.5lf",mytable->Get());
			}
			else
			{
				cout << "\t" << mytable->Get(); //fprintf(file,"\t%.5lf",mytable->Get());
			}
		}
	}
	cout << "\n"; // fprintf(file,"\n");
	//predfile.close();
	return;
}

void SampleModelParameters(string modelPath,
		                   int niterations,LPTable dataTable,LPTable mytable,const char* predfilename,
                           int* vars,CData& TestDiscrete)
{
//,const char* coeffile)
	int i,j;
	
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
	ThetafromP(Theta,P,localmodel);
	SaveCoef(Theta); //SaveCoef(coeffile,Theta);
	
	for(int iter=1;iter<niterations;iter++)
	{
		sampleNewThetas(VarSetsLocal,lenVarSetsLocal,nVarSetsLocal,
		                localmodel,localmodelgenerators,mytable,Y,Theta,junkTheta);
		SaveCoef(Theta);//SaveCoef(coeffile,Theta);				
		PfromTheta(P,Theta);
		SavePredictions(P,dataTable,predfilename,vars,TestDiscrete);
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
	return;
}

void getMarginal(LPTable dataTable,int lenmodel,int* model,CData& Discrete)
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
    
   for(i=0;i<Discrete.SampleSize;i++)
   {
      for(j=0;j<nvars;j++)
      {
	     index[j] = (int) Discrete.data[i][vars[j]];
      }
      dataTable->SetIndex(index);
      dataTable->Set(dataTable->Get()+1);
   }

   delete[] index; index = NULL;
   delete[] vars; vars = NULL;
   return;
}

void PrintCoefDescription(FILE* file,int ResponseVariable,LPTable mytable)
{
	int i,j;
	mytable->GetFirst();
	j = 0;
	while(mytable->GetNext())
	{
		if(mytable->Index[ResponseVariable-1])
		{
			j++;
			fprintf(file,"%d [",j);
			int first = 1;
			for(i=0;i<mytable->nDimens;i++)
			{
				if(i!=ResponseVariable-1)
				{
					if(mytable->Index[i])
					{
						if(first)
						{
							first = 0;
							fprintf(file,"%d",i+1);
						}
						else
						{
							fprintf(file,",%d",i+1);
						}
					}
				}
			}
			fprintf(file,"]\n");
		}
	}
	return;
}

void doRJMCMC(int n,CData& TrainDiscrete,CData& TestDiscrete,double cPrior)
{
	int j;
	double s;
	char predfilename[2048];
	//char coeffilename[2048];
	//char buffer[2048];
	//FILE* predfile = NULL;
	//FILE* coeffile = NULL;
	double maxprob;
	int nrelmodels;
			
	//the response variable is assumed to be the last column in the data matrix		
	int responseVariable = TestDiscrete.NumberOfGenes-1;			
						
	sprintf(predfilename,"%s.%d.fitted",params.mstrTestDataFile.c_str(),responseVariable);

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
	
	int nSampleSize = TestDiscrete.SampleSize;
	//double* predictions = new double[nSampleSize];
	//int lastrecord = 0;

	predfile << TestDiscrete.data[0][responseVariable]; //fprintf(predfile,"%.3lf",TestDiscrete.data[0][responseVariable]);
	for(j=1;j<nSampleSize;j++)
	{
		predfile << "\t" << TestDiscrete.data[j][responseVariable]; //fprintf(predfile,"\t%.3lf",TestDiscrete.data[j][responseVariable]);
	}
	predfile << "\n"; //fprintf(predfile,"\n");

	predfile.close();	
	//params.mnChainIterations
	CModel* amodel = models->Next;
	int modelid = 1;
	int totaliterations = 0;

	while(NULL!=amodel)
	{
		int niterations = (int)rint(params.mnChainIterations*amodel->laplace);
		printf("Sampling model [%d] for [%d] iterations.\n",modelid,niterations);
		totaliterations += niterations;
		modelid++;
		if(niterations>=1)
		{
			LPTable dataTable = new Table;
			getMarginal(dataTable,amodel->Length,amodel->vars,TrainDiscrete);
			LPTable datapriorTable = new Table;
			if(!datapriorTable->Alloc(dataTable->Dimens,dataTable->nDimens))
			{
				printf("Failed to allocate data+prior table.\n"); exit(1);
			}
			s = cPrior;
			if(s<0)
			{
				s = 1.0/((double)datapriorTable->Total);
			}
			for(j=0;j<datapriorTable->Total;j++)
			{
				datapriorTable->Data[j] = dataTable->Data[j]+s;
			}
			
			SampleModelParameters(amodel->modelPath,
								  niterations,dataTable,datapriorTable,predfilename,
			                      amodel->vars,TestDiscrete);//,stdout);
								
			datapriorTable->Reset(); delete datapriorTable; datapriorTable = NULL;					
			dataTable->Reset(); delete dataTable; dataTable = NULL;
		}
		amodel = amodel->Next;
	}
	printf("Total iterations = %d\n",totaliterations);
	
	/*
	//just use the best model
	CModel* amodel = models->getBestModel();
	LPTable dataTable = new Table;
	getMarginal(dataTable,amodel->Length,amodel->vars,TrainDiscrete);
	LPTable datapriorTable = new Table;
	if(!datapriorTable->Alloc(dataTable->Dimens,dataTable->nDimens))
	{
		printf("Failed to allocate data+prior table.\n"); exit(1);
	}
	s = cPrior;
	if(s<0)
	{
		s = 1.0/((double)datapriorTable->Total);
	}
	for(j=0;j<datapriorTable->Total;j++)
	{
		datapriorTable->Data[j] = dataTable->Data[j]+s;
	}
	
	sprintf(coeffilename,"%s.%d.coefdescription.txt",params.mstrTestDataFile.c_str(),TestDiscrete.NumberOfGenes);
	if(NULL==(coeffile=fopen(coeffilename,"w")))
	{
		fprintf(stderr,"Cannot open file [%s]\n",coeffilename);
		return;
	}
	PrintCoefDescription(coeffile,dataTable->nDimens,dataTable);
	fclose(coeffile);
	
	sprintf(coeffilename,"%s.%d.coef.txt",params.mstrTestDataFile.c_str(),TestDiscrete.NumberOfGenes);
	if(NULL==(coeffile=fopen(coeffilename,"w")))
	{
		fprintf(stderr,"Cannot open file [%s]\n",coeffilename);
		return;
	}
	
	SampleModelParameters(amodel->modelPath,
						  params.mnChainIterations,dataTable,datapriorTable,predfilename,
	                      amodel->vars,TestDiscrete);//coeffile);
	fclose(coeffile);
															
	datapriorTable->Reset(); delete datapriorTable; datapriorTable = NULL;					
	dataTable->Reset(); delete dataTable; dataTable = NULL;
	*/
	
	models->DeleteAll();	
	//fclose(predfile);
	//clean memory
	//delete[] predictions; predictions = NULL;
	delete models; models = NULL;
	return;
}


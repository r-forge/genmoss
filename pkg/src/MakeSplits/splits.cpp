/*
	We assume that the last column of the data contains the response variable
	and that the response variable is binary.
*/

#include "splits.h"

int subset(int n,int* set1,int* set2)
{
	int i;
	
	for(i=0;i<n;i++)
	{
		if(set1[i]>set2[i])
		{
			return(0);
		}
	}
	return(1);
}

int qsortcompare(const void * a, const void * b)
{
	double ia = *(double*)a;
	double ib = *(double*)b;
	
	if(ia<ib) return(-1);
	if(ia>ib) return(1);
	return(0);
}

void PrintTable(LPTable p)
{
	int i,j;
	
	p->GetFirst();
	for(i=0;i<p->Total;i++)
	{
		fprintf(stdout,"[%d,",p->Index[0]+1);
		for(j=1;j<p->nDimens-1;j++) fprintf(stdout,"%d,",p->Index[j]+1);
		fprintf(stdout,"%d] :: ",p->Index[p->nDimens-1]+1);
		//fprintf(stdout,"%.10lf\n",p->Data[i]);
		cout << p->Data[i] << "\n";
		p->GetNext();
	}
	return;
}

void MakeMarginals(LPTable Y,LPTable dataTable)
{
	int i;
	double s;
	int lenC = dataTable->nDimens;
	int* C = new int[lenC];

	if(!Y->Alloc(dataTable->Dimens,dataTable->nDimens))
	{
		printf("Failed to allocate marginals table.\n"); exit(1);
	}

	Y->GetFirst();
	Y->Set(dataTable->GetGrandTotal());
	while(Y->GetNext())
	{
		lenC = 0;
		for(i=0;i<Y->nDimens;i++)
		{
			if(Y->Index[i])
			{
				C[lenC] = i;
				lenC++;
			}
		}
		/*
		LPNTable p = new NTable;
		p->Create(lenC,C,dataTable);
		Y->Set(p->Data[p->Total-1]);
		
		p->Reset(); delete p; p = NULL;
		*/
		dataTable->Table::SetIndex(Y->Index);
		s = dataTable->Get();
		while(dataTable->GetNextWithIndex(lenC,C))
		{
			s += dataTable->Get();
		}
		Y->Set(s);
	}
	delete[] C; C = NULL;
	return;
}

void CreateShape(LPTable S,LPTable shape)
{
	int i,k;
	int lenF;
	int lenD;
	double s;
	int okay;
		
	for(i=0;i<shape->Total;i++)
	{
		shape->Data[i] = 0;
	}
			
	shape->GetFirst();
	s = 0.0;
	S->GetFirst();
	okay = 1;
	while(okay)
	{
		lenF = 0;
		for(k=0;k<S->nDimens;k++)
		{
			lenF += S->Index[k];
		}
		s += pow(-1.0,lenF)*S->Get();
		okay = S->GetNext();
	}
	shape->Set(s);
	
	int mokay = shape->GetNext();
	while(mokay)
	{
		lenD = 0;
		for(i=0;i<shape->nDimens;i++)
		{
			lenD += shape->Index[i];
		}
	
		s = 0.0;
		S->GetFirst();
		okay = 1;
		while(okay)
		{
			if(subset(S->nDimens,shape->Index,S->Index))
			{
				lenF = 0;
				for(k=0;k<S->nDimens;k++)
				{
					lenF += S->Index[k];
				}
				s += pow(-1.0,lenF-lenD)*S->Get();
			}
			okay = S->GetNext();
		}
		shape->Set(s);
		mokay = shape->GetNext();
	}
			
	return;
}

double ExactPost(LPTable marginal)
{
	int i;
	double mypost = 0;
		
	LPTable Y = new Table;
	MakeMarginals(Y,marginal);
	LPTable Shape = new Table;
	if(!Shape->Alloc(marginal->Dimens,marginal->nDimens))
	{
		printf("Failed to allocate shape table.\n"); exit(1);
	}
	CreateShape(Y,Shape);
	
	for(i=0;i<Shape->Total;i++)
	{
		mypost += lgamma(Shape->Data[i]);
	}
	
	//clean memory
	Y->Reset(); delete Y; Y = NULL;
	Shape->Reset(); delete Shape; Shape = NULL;
	return(mypost);
}

void MakeDiscreteDataTrain(CData& Data,CData& Discrete,int* goodvars,double* splits, double** score, bool* flags)
{
	int i;
	double s;
	int yindex = Data.NumberOfGenes-1;

	//double* score = new double[Data.NumberOfGenes];
	//memset(score,0,Data.NumberOfGenes*sizeof(double));
	double* xsort = new double[Data.SampleSize];
        double* xunique = new double[Data.NumberOfGenes];
        memset(xunique,0,Data.NumberOfGenes*sizeof(double));

	//need a 2x2 table
	const int n = 2;
	int index[n];
	for(i=0;i<n;i++) index[i] = 2;
	
	LPTable pYX = new Table;
	pYX->Alloc(index,n);

	LPTable pX = new Table;
	pX->Alloc(index,1);
	
	LPTable pY = new Table;
	pY->Alloc(index,1);

	s = 1.0/((double)pYX->Total);
	for(i=0;i<pYX->Total;i++)
	{
		pYX->Data[i] = s;
	}
	double priorDoubleScore = ExactPost(pYX);

	cout << "priorDoubleScore = " << priorDoubleScore << "\n"; //printf("priorDoubleScore = %.5lf\n",priorDoubleScore);

	s = 1.0/((double)pY->Total);
	for(i=0;i<pY->Total;i++)
	{
		pY->Data[i] = s;
	}
	double priorSingleScore = ExactPost(pY);
	
	cout << "priorSingleScore = " << priorSingleScore << "\n"; //printf("priorSingleScore = %.5lf\n",priorSingleScore);
	for(i=0;i<Data.SampleSize;i++)
	{
		index[0] = (int) (Data.data[i][yindex]);
		pY->SetIndex(index);
		pY->Set(pY->Get()+1);
	}

	double Yscore = ExactPost(pY);

	cout << "Yscore = " << Yscore << "\n"; //printf("Yscore = %.5lf\n",Yscore);
	
	double emptySeparatorScore = lgamma((double)(1+Data.SampleSize));
	
	cout << "emptySeparatorScore = " << emptySeparatorScore << "\n"; //printf("emptySeparatorScore = %.5lf\n",emptySeparatorScore);

	for(i=0;i<Data.SampleSize;i++)
	{
		Data.Y[i] = Data.data[i][yindex];
	}

	for(int agene=0;agene<yindex;agene++)
	{
		for(i=0;i<Data.SampleSize;i++)
		{
			xsort[i] = Data.data[i][agene];
		}
		qsort(xsort,Data.SampleSize,sizeof(double),qsortcompare);

		// Check for discrete case with only 3 data values
		xunique[0] = xsort[0];	
		int countunique = 1;
		for(i=0; i<Data.SampleSize-1; i++) 
		{
			if(xsort[i] != xsort[i+1]) {
				xunique[countunique] = xsort[i+1];
				countunique ++;
			}
		}
		// Additional variables to allow the fexibility of dealing with 3 vals		
		bool flagThree = false;
		double mmin = model.mdSplitsMin;
		double mmax = model.mdSplitsMax;
		double minc = model.mdSplitsInc;

//printf("count: %d\n", countunique);
		if(countunique <= 3) {
			flagThree = true;
			mmin = 0;
			mmax = 2;
			minc = 1;
		}

		
		int firstsplit = 1;
		double bestscore = -999999999;
		bool bestflag = false;
//		for(double a = model.mdSplitsMin;a<=model.mdSplitsMax;a+=model.mdSplitsInc)

		for(double a = mmin;a<=mmax;a+=minc)
		{ 
			double asplit = 0;

			// Define asplit and Data.X for the 3-valued case
			if(flagThree) {
				asplit = xunique[(int)a];
				for(i=0;i<Data.SampleSize;i++)
                                {
					if(a != 2) {
						if(Data.data[i][agene]<=asplit)
						{
							Data.X[i][0] = 0;
						}
						else
						{
							Data.X[i][0] = 1;
						}	
					} 
					else // additional setting: 2=1 and everything else to 0
					{
						if(Data.data[i][agene]==xunique[1])
                                                {
                                                        Data.X[i][0] = 1;
                                                }
                                                else
                                                {
                                                        Data.X[i][0] = 0;
                                                }
					}
                                }			
			}
			else 
			{
				//printf("should not be here!!!\n");
				asplit = xsort[(int)(a*(Data.SampleSize-1))];
				//printf("agene = %d :: asplit = %.3lf\n",agene,asplit);

				//this prevents all the samples being classified in one category
				if(asplit==xsort[Data.SampleSize-1])
				{
					continue;
				}

				for(i=0;i<Data.SampleSize;i++)
				{
					if(Data.data[i][agene]<=asplit)
					{
						Data.X[i][0] = 0;
					}
					else
					{
						Data.X[i][0] = 1;
					}
				}

			}			

			s = 1.0/((double)pYX->Total);
			for(i=0;i<pYX->Total;i++)
			{
				pYX->Data[i] = s;
			}
			s = 1.0/((double)pX->Total);
			for(i=0;i<pX->Total;i++)
			{
				pX->Data[i] = s;
			}
			
			for(i=0;i<Data.SampleSize;i++)
			{
				index[0] = (int) Data.Y[i];
				index[1] = (int) Data.X[i][0];
				pYX->SetIndex(index);
				pYX->Set(pYX->Get()+1);
			}
			
			for(i=0;i<Data.SampleSize;i++)
			{
				index[0] = (int) Data.X[i][0];
				pX->SetIndex(index);
				pX->Set(pX->Get()+1);
			}
			
			/*
			printf("pYX table\n");
			PrintTable(pYX);
			printf("pX table\n");
			PrintTable(pX);
			*/
			
			double score1 = ExactPost(pYX)-ExactPost(pX);
			if(a < 3)
				score[agene][(int)a] = score1;

			//printf("score = %.3lf\n",score1);
			if(firstsplit)
			{
				firstsplit = 0;
				bestscore = score1;
				for(i=0;i<Data.SampleSize;i++)
				{
					Discrete.data[i][agene] = Data.X[i][0];
				}
				splits[agene] = asplit;
			}
			else
			{
				if(score1>bestscore)
				{
					bestscore = score1;
					for(i=0;i<Data.SampleSize;i++)
					{
						Discrete.data[i][agene] = Data.X[i][0];
					}
					splits[agene] = asplit;
					// If different order is encoded, save the second variable as the split. 
					if(flagThree && a==2) 
					{
						//printf("Reversed order, %d\n", agene);
						bestflag = true;
						splits[agene] = xunique[1];
					} else
					{
						bestflag = false; // not necessary, since flag can only be activated at the very end of the loop over 'a'
					}
				}
			}
		}
		score[agene][3] = bestscore;
		flags[agene] = bestflag;
	}

//printf("****************** After formula: ***********************");

/*	for(i=0;i<yindex;i++)
	{
		double s1 = score[i]-priorDoubleScore;
		double s2 = Yscore-emptySeparatorScore-2*priorSingleScore;
		score[i] = s1-s2;
		if(i < 20)
			printf("new = %.3lf\n", score[i]);
		if(score[i]>0)
		{
			printf("good: %d, %.3lf\n", i, score[i]);
			//goodvars[i] = 1;
		}
	}*/
	/*double maxscore = score[0];
	for(i=1;i<yindex;i++)
	{
		if(maxscore<score[i])
		{
			maxscore = score[i];
		}
	}*/


	/*for(i=0;i<yindex;i++)
	{

		if(goodvars[i])
		{
			if(exp(score[i]-maxscore)>=model.mdCutoffMax)
			{
				goodvars[i] = 1;
			}
			else
			{
				goodvars[i] = 0;
			}
		}
	}*/
/*	goodvars[yindex] = 1;
*/	
	for(i=0;i<yindex;i++)
	{
		goodvars[i] = 1;
	}
	goodvars[yindex] = 1;

	//clean memory
	//delete[] score; score = NULL;
	delete[] xsort; xsort = NULL;
	delete[] xunique; xunique = NULL;
	pYX->Reset(); delete pYX; pYX = NULL;
	pX->Reset(); delete pX; pX = NULL;
	pY->Reset(); delete pY; pY = NULL;
	return;
}

void MakeDiscreteDataTest(CData& Data,CData& Discrete,double* splits, bool* flags)
{
        int i;
        //double s;
        int yindex = Data.NumberOfGenes-1;

        for(int agene=0;agene<yindex;agene++)
        {
		//if(flags[agene])
		//	printf("FLAG: %d\n", agene);

                for(i=0;i<Data.SampleSize;i++)
                {
                       /* if(Data.data[i][agene]<=splits[agene])
                        {
                                Discrete.data[i][agene] = 0;
                        }
                        else
                        {
                                Discrete.data[i][agene] = 1;
                        }*/

			// In 3-valued discrete case, if the split is of different
			// ordering, then the split represents value which should be
			// set to 1, and everything else is 0. 
			if(flags[agene]) {
				if(Data.data[i][agene]==splits[agene])
				{
					Discrete.data[i][agene] = 1;
				}
				else
				{
					Discrete.data[i][agene] = 0;
				}
			}
			else
			{
				if(Data.data[i][agene]<=splits[agene])
				{
					Discrete.data[i][agene] = 0;
				}
				else
				{
					Discrete.data[i][agene] = 1;
				}
			}

                }
        }
        return;
}

LPTable getMarginal(int* goodvars,CData& Discrete)
{
   int i,j;

   LPTable mytable = new Table;
   
   int nvars = 0;
   for(j=0;j<Discrete.NumberOfGenes;j++)
   {
      if(goodvars[j])
	  {
		 nvars++;
	  }
   }
   int* index = new int[nvars];
   for(j=0;j<nvars;j++)
   {
      index[j] = 2;
   }
   
   printf("nvars = %d\n",nvars);
   if(!mytable->Alloc(index,nvars))
   {
      printf("Error allocating memory!\n");
      exit(1);
   }
   int* vars = new int[nvars];
   i = 0;
   for(j=0;j<Discrete.NumberOfGenes;j++)
   {
      if(goodvars[j])
	  {
		 vars[i] = j;
		 i++;
	  }
   }
   
   for(i=0;i<Discrete.SampleSize;i++)
   {
      for(j=0;j<nvars;j++)
      {
	     index[j] = (int) (Discrete.data[i][vars[j]]);
      }
      mytable->SetIndex(index);
      mytable->Set(mytable->Get()+1);
   }

   delete[] index; index = NULL;
   delete[] vars; vars = NULL;
   return(mytable);
}


int splits_main(char* fname)
{
	if (fname == NULL)
	{
		cout << "USAGE: " << "[modelfile]" << endl;
		exit(1);
	}
	if (!model.Load(fname)) {
	  cout << model.GetErrorMessage() << endl;
	  exit(1);
	}

	//vslNewStream(&mystream,BRNG,SEED);
	mystream = gsl_rng_alloc(BRNG);
	gsl_rng_set(mystream, SEED);

	//read the train data
	TrainData.NumberOfGenes = model.mnNumberOfVariables;
	TrainData.SampleSize = model.mnTrainSampleSize;   
	TrainData.DataFile = model.mstrTrainDataFile;
	TrainData.ReadData();
	TrainData.Allocate(2);

	//read the test data
	TestData.NumberOfGenes = model.mnNumberOfVariables;
	TestData.SampleSize = model.mnTestSampleSize;   
	TestData.DataFile = model.mstrTestDataFile;
	TestData.ReadData();
	TestData.Allocate(2);

	
	TrainDiscrete.NumberOfGenes = TrainData.NumberOfGenes;
	TrainDiscrete.SampleSize = TrainData.SampleSize;   
	TrainDiscrete.DataFile = TrainData.DataFile;
	TrainDiscrete.ReadData();

	TestDiscrete.NumberOfGenes = TestData.NumberOfGenes;
	TestDiscrete.SampleSize = TestData.SampleSize;   
	TestDiscrete.DataFile = TestData.DataFile;
	TestDiscrete.ReadData();

	int* goodvars = new int[TrainDiscrete.NumberOfGenes];
	memset(goodvars,0,TrainDiscrete.NumberOfGenes*sizeof(int));
	double* trainsplits = new double[TrainDiscrete.NumberOfGenes];
	memset(trainsplits,0,TrainDiscrete.NumberOfGenes*sizeof(double));
	//double* testsplits = new double[TestDiscrete.NumberOfGenes];
	//memset(testsplits,0,TestDiscrete.NumberOfGenes*sizeof(double));

	// Define score estimates to output them to file 
	//double* score = new double[TrainData.NumberOfGenes];
	//memset(score,0,TrainData.NumberOfGenes*sizeof(double));

	int scorecols = 4;
	double** score = new double*[TrainData.NumberOfGenes];
	memset(score,0,TrainData.NumberOfGenes*sizeof(double));

	for (int i = 0; i < TrainData.NumberOfGenes; i++) {
		score[i] = new double[scorecols];
		memset(score[i], 0, scorecols*sizeof(double));
	}



	// Define flags: True at genes where 3 discrete data values are used  
	// AND the unusual ordering is taking place: 
	// where second value is 1, and everything else is 0
	bool* flags = new bool[TrainData.NumberOfGenes];
	memset(flags,false,TrainData.NumberOfGenes*sizeof(bool));

	//MakeDiscreteData(TrainData,TrainDiscrete,goodvars,trainsplits);
	//MakeDiscreteData(TestData,TestDiscrete,goodvars,testsplits);

        //we determine the splits in the training data
        MakeDiscreteDataTrain(TrainData,TrainDiscrete,goodvars,trainsplits, score, flags);
        //we discretize the test data using the splits from the training data
        MakeDiscreteDataTest(TestData,TestDiscrete,trainsplits,flags);

	//save the discrete training data and the corresponding splits
	TrainDiscrete.DataFile = model.mstrTrainOutputFile;
	TrainDiscrete.WriteData(goodvars);

	// Write the scores
	ofstream out (model.mstrScoresFile.c_str()); //FILE* out = fopen(model.mstrScoresFile.c_str(),"w");
	//if(NULL==out)
	if(!out.is_open())
	{
		printf("Cannot open data file %s.\n",model.mstrScoresFile.c_str());
		exit(1);
	}
	for(int i=0;i<(TrainData.NumberOfGenes-1);i++)
	{
		for(int j=0; j<scorecols; j++) 
			out << score[i][j] << "\t"; //fprintf(out,"%5.10lf\t",score[i][j]);
		out << "\n"; //fprintf(out, "\n");
	}
	out.close(); //fclose(out); 


	/*FILE* out = fopen(model.mstrTrainIndexFile.c_str(),"w");
	if(NULL==out)
	{
		printf("Cannot open data file %s.\n",model.mstrTrainIndexFile.c_str());
		exit(1);
	}
	for(int i=0;i<TrainDiscrete.NumberOfGenes;i++)
	{
		if(goodvars[i])
		{
			fprintf(out,"%d\t%.5lf\n",i+1,trainsplits[i]);
		}
	}
	fclose(out); */

	//save the discrete test data and the corresponding splits
	TestDiscrete.DataFile = model.mstrTestOutputFile;
	TestDiscrete.WriteData(goodvars);
	/*out = fopen(model.mstrTestIndexFile.c_str(),"w");
	if(NULL==out)
	{
		printf("Cannot open data file %s.\n",model.mstrTestIndexFile.c_str());
		exit(1);
	}
	for(int i=0;i<TestDiscrete.NumberOfGenes;i++)
	{
		if(goodvars[i])
		{
			fprintf(out,"%d\t%.5lf\n",i+1,testsplits[i]);
		}
	}
	fclose(out);*/

	//clean memory
	//mytable->Reset(); delete mytable; mytable = NULL;
	delete[] trainsplits; trainsplits = NULL;
	//delete[] testsplits; testsplits = NULL;
	delete[] goodvars; goodvars = NULL;

	for(int i = 0; i < TrainData.NumberOfGenes; i++) {
                delete[] score[i]; score[i] = NULL;
        }
	delete[] score; score = NULL;
	delete[] flags; flags = NULL;

	TrainData.Cleanup();
	TestData.Cleanup();
	TrainDiscrete.Cleanup();
	TestDiscrete.Cleanup();
	//vslDeleteStream(&mystream);
	gsl_rng_free(mystream);
	return(1);
}


extern "C" void rsplits(char **fname)
{
        pid_t pID = fork();

        if (pID == 0) {                 // child
                char *fname2 = fname[0];
                splits_main(fname2);
                exit(0);
        }
        else if (pID <0) {
                printf("Failed to fork to process the request.\n");
        }
        else {                          // parent
                int childstatus;
                waitpid(pID, &childstatus, 0);
        }
}

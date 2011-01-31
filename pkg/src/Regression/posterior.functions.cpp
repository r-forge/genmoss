#include "f2c.h"
#include "clapack.h"
#include "blaswrap.h"

void MakeMarginals(LPTable Y,LPTable dataTable)
{
	int i,j;
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
		
		dataTable->SetIndex(Y->Index);
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

void InitPriorTable(LPTable prior,LPTable data,double cPrior)
{
	int i;
	
	if(!prior->Alloc(data->Dimens,data->nDimens))
	{
		printf("Failed to allocate prior table.\n"); exit(1);
	}
	
	double s = cPrior;
	if(s<0)
	{
		s = 1.0/((double)data->Total);
	}
	for(i=0;i<prior->Total;i++)
	{
		prior->Data[i] = s;
	}

	return;
}

void InitThetaTable(int* amodel,LPTable Theta, gsl_rng * stream)
{
	int i;
	double r[Theta->Total];
	
	//vdRngGaussian(METHOD,stream,Theta->Total,r,0,M);
	for(i=0; i<Theta->Total; i++)
		r[i] = gsl_ran_gaussian(stream, M);

	Theta->Data[0] = 0;
	for(i=1;i<Theta->Total;i++)
	{
		if(1==amodel[i])
		{
			Theta->Data[i] = r[i];
		}
		else
		{
			Theta->Data[i] = 0;
		}
	}

	return;
}

double probYgivenG(LPTable Y,LPTable Theta,int* amodel)
{
	int i,j;
	double lprob = 0.0;
	double emptyterm;
	double s;
	
	for(i=1;i<Y->Total;i++)
	{
		if(amodel[i])
		{
			lprob += (Y->Data[i]*Theta->Data[i]);
		}
	}
	
	emptyterm = 1.0;
	
	Y->GetFirst();
	while(Y->GetNext())
	{
		int thereisone = 0;
		s = 0.0;
		j = 1;
		Theta->GetFirst();
		while(Theta->GetNext())
		{
			if(amodel[j])
			{
				if(subset(Theta->nDimens,Theta->Index,Y->Index))
				{
					s += Theta->Get();
					thereisone = 1;
				}
			}
			j++;
		}
		if(thereisone)
		{
			emptyterm += exp(s);
		}
	}	
	lprob = lprob-(Y->Data[0]*log(emptyterm));
	
	return(lprob);
}


int GetNextIndexCond(int n,int* index,LPNTable shape)
{
	int i,j;
		
	for(i=0;i<n;i++)
	{
		if(0==index[i])
		{
			index[i]=1;
			return(1);
		}
		else if(0==shape->Index[i])
		{
			index[i]=0;
		}
	}
	return(0);
}

void CreateShape(int lenC,int* C,LPTable S,LPNTable shape,
                 int** VarSets,int* lenVarSets,int nVarSets)
{
	int i,j,k;
	int indexC[S->nDimens];
	int lenF;
	int lenD;
	double s;
	int okay;
	
	for(i=0;i<S->nDimens;i++)
	{
		indexC[i]=0;
	}
	for(i=0;i<lenC;i++)
	{
		indexC[C[i]] = 1;
	}
	
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
		if(subset(S->nDimens,S->Index,indexC))
		{
			lenF = 0;
			for(k=0;k<S->nDimens;k++)
			{
				lenF += S->Index[k];
			}
			s += pow(-1.0,lenF)*S->Get();
		}
		okay = S->GetNext();
	}
	shape->Set(s);
	
	for(i=0;i<nVarSets;i++)
	{
		if(subset(S->nDimens,VarSets[i],indexC))
		{
			s = 0.0;
			S->GetFirst();
			okay = 1;
			while(okay)
			{
				if(subset(S->nDimens,S->Index,indexC))
				{
					if(subset(S->nDimens,VarSets[i],S->Index))
					{
						lenD = lenVarSets[i];
						lenF = 0;
						for(k=0;k<S->nDimens;k++)
						{
							lenF += S->Index[k];
						}
						s += pow(-1.0,lenF-lenD)*S->Get();
					}
				}
				okay = S->GetNext();
			}
			shape->SetIndex(VarSets[i]);
			shape->Set(s);
		}
	}
			
	return;
}


void PfromTheta(LPTable P,LPTable Theta)
{
	int i,j;
	double s;
	double total;
	
	total = 0.0;
	P->GetFirst();
	while(P->GetNext())
	{
		Theta->GetFirst();
		s = 0.0;
		while(Theta->GetNext())
		{
			if(subset(Theta->nDimens,Theta->Index,P->Index))
			{
				s += Theta->Get();
			}
		}
		P->Set(exp(s));
		total += P->Get();
	}
	P->Data[0] = 1.0/(1.0+total);
	for(i=1;i<P->Total;i++)
	{
		P->Data[i] *= P->Data[0];
	}
	
	return;
}

void ThetafromP(LPTable theta,LPTable p,int* aModel)
{
	int i,j;
	double s;
	
	theta->GetFirst();
	theta->Set(log(p->Data[0]));
	j = 1;
	while(theta->GetNext())
	{
		if(aModel[j])
		{
			int lenE = 0;
			for(i=0;i<theta->nDimens;i++) lenE+=theta->Index[i];
		
			p->GetFirst();
			s = pow(-1,lenE)*log(p->Data[0]);
			while(p->GetNext())
			{
				if(subset(p->nDimens,p->Index,theta->Index))
				{
					int lenF = 0;
					for(i=0;i<p->nDimens;i++) lenF+=p->Index[i];
				
					s+= pow(-1.0,lenE-lenF)*log(p->Get());
				}
			}
		}
		else
		{
			s = 0;
		}
		theta->Set(s);
		j++;
	}

	return;
}


void getThetaMarginal(LPNTable theta,LPNTable p)
{
	int i;
	double s;
	
	theta->GetFirst();
	theta->Set(log(p->Data[0]));
	while(theta->GetNext())
	{
		int lenE = 0;
		for(i=0;i<theta->nDimens;i++) lenE+=theta->Index[i];
		
		p->GetFirst();
		s = pow(-1,lenE)*log(p->Data[0]);
		while(p->GetNext())
		{
			if(subset(p->nDimens,p->Index,theta->Index))
			{
				int lenF = 0;
				for(i=0;i<p->nDimens;i++) lenF+=p->Index[i];
				
				s+= pow(-1.0,lenE-lenF)*log(p->Get());
			}
		}
		theta->Set(s);
	}

	return;
}

void findComplete(int n,int nVarSets,int** VarSets,int* lenVarSets,LPGraph G,int* isComplete)
{
	int i,j,k;
	
	for(i=0;i<nVarSets;i++)
	{
		isComplete[i] = 1;
		for(j=0;j<n-1;j++)
		{
			if(VarSets[i][j])
			{
				for(k=j+1;k<n;k++)
				{
					if(VarSets[i][k])
					{
						if(0==G->Edge[j][k])
						{
							isComplete[i] = 0;
							break;
						}
					}
				}
			}
			if(0==isComplete[i])
			{
				break;
			}
		}
	}

	//the empty set is not complete
	isComplete[0] = 0;
	//make sure all the main effects are in the model
	for(i=0;i<nVarSets;i++)
	{
		if(1==lenVarSets[i])
		{
			isComplete[i] = 1;
		}
	}

	return;
}

void ipf(int** VarSetsMarg,int* lenVarSetsMarg,int nVarSetsMarg,
         LPTable P,LPTable dataTable,int* aModelGenerators,int* aModel,
		 gsl_rng * stream)
{
	int i,j,k;
	double delta = 0.0000001;
	double m[P->Total];
	double mold[P->Total];
	double s;
	
	int nMssShape = 0;
	for(i=0;i<nVarSetsMarg;i++)
	{
		nMssShape += aModelGenerators[i];
	}

	int* lenC = new int[nMssShape];
	int** C = new int*[nMssShape];
	for(i=0;i<nMssShape;i++)
	{
		C[i] = new int[dataTable->nDimens];
	}
	
	LPNTable margin = new NTable[nMssShape];
	nMssShape = 0;
	for(i=0;i<nVarSetsMarg;i++)
	{
		if(aModelGenerators[i])
		{
			k = 0;
			for(j=0;j<dataTable->nDimens;j++)
			{
				if(VarSetsMarg[i][j])
				{
					C[nMssShape][k]=j;
					k++;
				}
			}
			lenC[nMssShape] = k;
		
			margin[nMssShape].Create(lenC[nMssShape],C[nMssShape],dataTable);
			nMssShape++;
		}
	}
	
	LPTable Theta = new Table;
	if(!Theta->Alloc(dataTable->Dimens,dataTable->nDimens))
	{
		printf("Failed to allocate the theta's table.\n"); exit(1);
	}
	
	LPTable oldTheta = new Table;
	if(!oldTheta->Alloc(dataTable->Dimens,dataTable->nDimens))
	{
		printf("Failed to allocate the theta's table.\n"); exit(1);
	}
	
	//generate a random theta
	InitThetaTable(aModel,Theta,stream);
	PfromTheta(P,Theta);
	
	int notdone = 1;
	while(notdone)
	{
		for(j=0;j<P->Total;j++)
		{
			mold[j] = P->Data[j];
			oldTheta->Data[j] = Theta->Data[j];
		}
		for(i=0;i<nMssShape;i++)
		{
			LPNTable p = new NTable;
			p->Create(lenC[i],C[i],P);
		
			P->GetFirst();
			j=0;
			m[j] = P->Data[j]*margin[i].GetI(P->Index)/p->GetI(P->Index);
			while(P->GetNext())
			{
				j++;
				m[j] = P->Data[j]*margin[i].GetI(P->Index)/p->GetI(P->Index);
			}
			
			for(j=0;j<P->Total;j++)
			{
				P->Data[j] = m[j]; 
			}
			
			p->Reset(); delete p ; p = NULL;

			s = 0.0;
			for(j=0;j<P->Total;j++)
			{
				s += P->Data[j];
			}
	
			for(j=0;j<P->Total;j++)
			{
				P->Data[j] /= s;
			}

			//transform in the space of thetas and impose the constraints
			ThetafromP(Theta,P,aModel);
			PfromTheta(P,Theta);
		}
	
		notdone = 0;
		for(i=0;i<P->Total;i++)
		{
			/*
			if(fabs(mold[i]-P->Data[i])>delta)
			{
				notdone = 1;
			}
			*/
			//printf("\t%.5lf",Theta->Data[i]);
			if(fabs(oldTheta->Data[i]-Theta->Data[i])>delta)
			{
				notdone = 1;
			}
		}
		//printf("\n");
	}
	
	s = 0.0;
	for(i=0;i<P->Total;i++)
	{
		s += P->Data[i];
	}
	
	for(i=0;i<P->Total;i++)
	{
		P->Data[i] /= s;
	}
	
	//clean memory
	oldTheta->Reset(); delete oldTheta; oldTheta = NULL;
	Theta->Reset(); delete Theta; Theta = NULL;
	for(i=0;i<nMssShape;i++)
	{
		margin[i].Reset();
	}
	delete[] margin; margin = NULL;
	for(i=0;i<nMssShape;i++)
	{
		delete[] C[i]; C[i] = NULL;
	}
	delete[] C; C = NULL;
	delete[] lenC; lenC = NULL;
	return;
}

double hessian(int** VarSetsMarg,int* lenVarSetsMarg,int nVarSetsMarg,
               int* isComplete,LPTable Y,LPTable Theta)
{
	double logdet = 0.0;
	int i,j,k;
	int d1, indexd1;
	int d2, indexd2;
	double s;
	
	double emptyterm = 1.0;
	Y->GetFirst();
	while(Y->GetNext())
	{
		int thereisone = 0;
		s = 0.0;
		j = 1;
		Theta->GetFirst();
		while(Theta->GetNext())
		{
			if(isComplete[j])
			{
				if(subset(Theta->nDimens,Theta->Index,Y->Index))
				{
					s += Theta->Get();
					thereisone = 1;
				}
			}
			j++;
		}
		if(thereisone)
		{
			emptyterm += exp(s);
		}
	}	
	
	int p = 0;
	for(i=0;i<Theta->Total;i++)
	{
		if(isComplete[i]) p++;
	}

	double single[p];
	double couple[p][p];
	double H[p][p];
	for(i=0;i<p;i++)
	{
		single[i]=0;
		for(j=0;j<p;j++)
		{
			couple[i][j] = 0;
			H[i][j] = 0;
		}
	}
	
	indexd1 = 0;
	for(d1=0;d1<Theta->Total;d1++)
	{
		if(isComplete[d1])
		{
			single[indexd1] = 0;
			Y->GetFirst();
			while(Y->GetNext())
			{
				if(subset(Theta->nDimens,VarSetsMarg[d1],Y->Index))
				{
					int thereisone = 0;
					s = 0.0;
					j=1;
					Theta->GetFirst();
					while(Theta->GetNext())
					{
						if(isComplete[j])
						{
							if(subset(Theta->nDimens,Theta->Index,Y->Index))
							{
								s += Theta->Get();
								thereisone = 1;
							}
						}
						j++;
					}
					if(thereisone)
					{
						single[indexd1] += exp(s);
					}
				}
			}
			indexd1++;
		}
	}
	//printf("single :: \n");
	//for(i=0;i<p;i++) printf("%.6lf\t",single[i]);
	//printf("\n");
	
	indexd1=0;
	for(d1=0;d1<Theta->Total;d1++)
	{
		if(isComplete[d1])
		{
			indexd2=0;
			for(d2=0;d2<nVarSetsMarg;d2++)
			{
				if(isComplete[d2])
				{
					if(d1<d2)
					{
						couple[indexd1][indexd2] = 0;
						Y->GetFirst();
						while(Y->GetNext())
						{
							if(subset(Theta->nDimens,VarSetsMarg[d1],Y->Index))
							{
								if(subset(Theta->nDimens,VarSetsMarg[d2],Y->Index))
								{
									int thereisone = 0;
									s = 0.0;
									j=1;
									Theta->GetFirst();
									while(Theta->GetNext())
									{
										if(isComplete[j])
										{
											if(subset(Theta->nDimens,Theta->Index,Y->Index))
											{
												s += Theta->Get();
												thereisone = 1;
											}
										}
										j++;
									}
									if(thereisone)
									{
										couple[indexd1][indexd2] += exp(s);
									}
								}
							}
						}
					}
					couple[indexd2][indexd1] = couple[indexd1][indexd2];
					indexd2++;
				}
			}
			indexd1++;
		}
	}
	
	for(d1=0;d1<p;d1++)
	{
		H[d1][d1] = -Y->Data[0]*pow(emptyterm,-2)*pow(single[d1],2)+Y->Data[0]*pow(emptyterm,-1)*single[d1];
		for(d2=d1+1;d2<p;d2++)
		{
			H[d1][d2] = -Y->Data[0]*pow(emptyterm,-2)*single[d1]*single[d2]+Y->Data[0]*pow(emptyterm,-1)*couple[d1][d2];
			H[d2][d1] = H[d1][d2];
		}
	}
	
	/*
	printf("H ::\n");
	for(i=0;i<p;i++)
	{
		for(j=0;j<p;j++)
		{
			printf("%.6lf\t",H[i][j]);
		}
		printf("\n");
	}
	*/

	char jobvl = 'N';
	char jobvr = 'N';
	integer lda = (integer)p;
	double wr[2*p];
	double wi[2*p];
	double vl[p][p];
	integer ldvl = (integer)(p*p);
	double vr[p][p];
	integer ldvr = (integer)(p*p);
	double work[p*p];
	integer lwork = (integer)(p*p);
	integer info;
	integer pp = (integer)p;

	dgeev_(&jobvl,&jobvr,&pp,(double*)H,&lda,(double*)wr,(double*)wi,(double*)vl, 
		   &ldvl,(double*)vr,&ldvr,(double*)work,&lwork,&info);
	
	for(i=0;i<p;i++)
	{
		logdet+=log(wr[i]);
		//printf("logdet = %.5lf :: wr[%d] = %.5lf\n",logdet,i+1,wr[i]);
	}
	
	return(logdet);
}

double laplace(int** VarSetsMarg,int* lenVarSetsMarg,int nVarSetsMarg,
               int* currentModel,int* currentModelGenerators,LPTable dataTable,LPTable Ydata,
			   gsl_rng * stream)
{
	int i;
	double dev = 0;
	double nSampleSize = dataTable->GetGrandTotal();
							
	int nParams = 0;
	for(i=0;i<nVarSetsMarg;i++)
	{
		nParams += currentModel[i];
	}
	//printf("nParams = %d :: samplesize = %.2lf\n",nParams,nSampleSize);
	
	LPTable Theta = new Table;
	if(!Theta->Alloc(dataTable->Dimens,dataTable->nDimens))
	{
		printf("Failed to allocate the theta's table.\n"); exit(1);
	}
	LPTable P = new Table;
	if(!P->Alloc(dataTable->Dimens,dataTable->nDimens))
	{
		printf("Failed to allocate the probs' table.\n"); exit(1);
	}
	
	LPTable junkTheta = new Table;
	if(!junkTheta->Alloc(dataTable->Dimens,dataTable->nDimens))
	{
		printf("Failed to allocate the theta's table.\n"); exit(1);
	}
	
	LPTable nextTheta = new Table;
	if(!nextTheta->Alloc(dataTable->Dimens,dataTable->nDimens))
	{
		printf("Failed to allocate the nextTheta's table.\n"); exit(1);
	}
	
	ipf(VarSetsMarg,lenVarSetsMarg,nVarSetsMarg,
		P,dataTable,currentModelGenerators,currentModel,
		stream);
	ThetafromP(Theta,P,currentModel);
		
	//printf("loglik = %.3lf\n",probYgivenG(Ydata,Theta,currentModel));
	double loglik = probYgivenG(Ydata,Theta,currentModel); 
	dev = loglik-0.5*hessian(VarSetsMarg,lenVarSetsMarg,nVarSetsMarg,
							 currentModel,dataTable,Theta);
	
	//printf("dev = %.5lf\n",dev);
			
	//clean memory
	P->Reset(); delete P; P = NULL;
	Theta->Reset(); delete Theta; Theta = NULL;
	junkTheta->Reset(); delete junkTheta; junkTheta = NULL;
	nextTheta->Reset(); delete nextTheta; nextTheta = NULL;
	return(dev);
}

double getLaplace(int* currentModel,int* currentModelGenerators,
                   LPTable priorTable,LPTable Yprior,
				   LPTable datapriorTable,LPTable Ydataprior,
				   int** VarSets,int* lenVarSets,int nVarSets,
				   gsl_rng * stream)
{
	double prior = laplace(VarSets,lenVarSets,nVarSets,
	                       currentModel,currentModelGenerators,priorTable,Yprior,
						   stream);
	double post = laplace(VarSets,lenVarSets,nVarSets,
						  currentModel,currentModelGenerators,datapriorTable,Ydataprior,
						  stream);
	return(post-prior);
}

void doRJMCMCstart(CModel* models,CModel* localmodels,int startpoint,
                   LPTable dataTable,LPTable priorTable,LPTable mytable,
				   int** VarSets,int* lenVarSets,int nVarSets,
				   int** DownLinks,int* nDownLinks,int** UpLinks,int* nUpLinks,
				   gsl_rng * stream)
{
	int i,j,k;
	int iteration = -1;
	double s;
				
	int n = dataTable->nDimens;			
	int* currentModel = new int[nVarSets];
	int* currentModelGenerators = new int[nVarSets];
	int* currentModelDualGenerators = new int[nVarSets];
	int* currentAddDel = new int[nVarSets];
	int nCurrentAddDel = -1;
	
	int* nextModel = new int[nVarSets];
	int* nextModelGenerators = new int[nVarSets];
	int* nextModelDualGenerators = new int[nVarSets];
	int* nextAddDel = new int[nVarSets];
	int nNextAddDel = -1;
	
	LPTable Ydata = new Table;
	MakeMarginals(Ydata,mytable);
	
	LPTable Yprior = new Table;
	MakeMarginals(Yprior,priorTable);
	
	//generate a random starting point
	RandomModel(currentModel,nVarSets,VarSets,lenVarSets,rndModelBern,stream);
	MakeModelHierarchical(nVarSets-1,currentModel[nVarSets-1],currentModel,DownLinks,nDownLinks);
	for(i=0;i<nVarSets;i++)
	{
		nextModel[i] = currentModel[i];
	}
	
	//find the generators of the current model
	memset(currentModelGenerators,0,nVarSets*sizeof(int));
	findGenerators(0,currentModel,nVarSets,UpLinks,nUpLinks,currentModelGenerators);
	
	//find the dual generators of the current model
	memset(currentModelDualGenerators,0,nVarSets*sizeof(int));
	findDualGenerators(nVarSets-1,currentModel,nVarSets,DownLinks,nDownLinks,currentModelDualGenerators);
					
	memset(currentAddDel,0,nVarSets*sizeof(int));
	nCurrentAddDel = 0;
	for(i=0;i<nVarSets;i++)
	{
		if(lenVarSets[i]>=2) //cannot add or delete the main effects
		{
			if(currentModelGenerators[i]||currentModelDualGenerators[i])
			{
				currentAddDel[nCurrentAddDel] = i;
				nCurrentAddDel++;
			}
		}
	}
	
	s = getLaplace(currentModel,currentModelGenerators,priorTable,Yprior,mytable,Ydata,
	               VarSets,lenVarSets,nVarSets,stream);	
	models->AddLaplace(nVarSets,currentModel,currentModelGenerators,s);
	localmodels->AddLaplace(nVarSets,currentModel,currentModelGenerators,s);		
						
	iteration = 1;
	while(1)
	{
		printf("startpoint[%d] iteration[%d]\n",startpoint,iteration);
		iteration++;
		
		//propose a model
		for(i=0;i<nVarSets;i++)
		{
			nextModel[i] = currentModel[i];
		}
		
		//add all the neighbors in the list of models
		for(k=0;k<nCurrentAddDel;k++)
		{
			int a = currentAddDel[k];
			nextModel[a] = 1-currentModel[a];
		
			//find the generators of the next model
			memset(nextModelGenerators,0,nVarSets*sizeof(int));
			findGenerators(0,nextModel,nVarSets,UpLinks,nUpLinks,nextModelGenerators);
			int ngenerators = 0;
			for(i=0;i<nVarSets;i++)
			{
				ngenerators += nextModelGenerators[i];
			}
			if(ngenerators>=2)
			{
				CModel* nextmodel = localmodels->Find(nVarSets,nextModel);
				if(NULL==nextmodel)
				{
					s = getLaplace(nextModel,nextModelGenerators,priorTable,Yprior,mytable,Ydata,
					               VarSets,lenVarSets,nVarSets,stream);
					localmodels->AddLaplace(nVarSets,nextModel,nextModelGenerators,s);
				}
			}
			nextModel[a] = 1-nextModel[a];
		}
				
		localmodels->DeleteNotRelevant();
		printf("%.5lf\n",localmodels->laplace);
		CModel* nextmodel = localmodels->Shotgun(stream,stdout);
		if(NULL==nextmodel)
		{
			break;//done with the search
		}
		for(i=0;i<nVarSets;i++)
		{
			currentModel[i] = nextmodel->vars[i];
			currentModelGenerators[i] = nextmodel->generators[i];
		}
		memset(currentModelDualGenerators,0,nVarSets*sizeof(int));
		findDualGenerators(nVarSets-1,currentModel,nVarSets,DownLinks,nDownLinks,currentModelDualGenerators);
		memset(currentAddDel,0,nVarSets*sizeof(int));
		nCurrentAddDel = 0;
		for(i=0;i<nVarSets;i++)
		{
			if(lenVarSets[i]>=2) //cannot add or delete the main effects
			{
				if(currentModelGenerators[i]||currentModelDualGenerators[i])
				{
					currentAddDel[nCurrentAddDel] = i;
					nCurrentAddDel++;
				}
			}
		}
		models->AddLaplace(nVarSets,currentModel,currentModelGenerators,nextmodel->laplace);
	}
	
	//clean memory
	Yprior->Reset(); delete Yprior; Yprior = NULL;
	Ydata->Reset(); delete Ydata; Ydata = NULL;
	delete[] currentModelDualGenerators; currentModelDualGenerators = NULL;
	delete[] currentModelGenerators; currentModelGenerators = NULL;
	delete[] currentModel; currentModel = NULL;
	delete[] currentAddDel; currentAddDel = NULL;
	delete[] nextModelDualGenerators; nextModelDualGenerators = NULL;
	delete[] nextModelGenerators; nextModelGenerators = NULL;
	delete[] nextModel; nextModel = NULL;
	delete[] nextAddDel; nextAddDel = NULL;
	return;
}

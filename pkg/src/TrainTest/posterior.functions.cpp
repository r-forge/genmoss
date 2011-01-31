void MakeMarginals(LPTable Y,LPTable dataTable)
{
	int i;
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
		LPNTable p = new NTable;
		p->Create(lenC,C,dataTable);
		Y->Set(p->Data[p->Total-1]);
		
		p->Reset(); delete p; p = NULL;
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

void InitThetaTable(int* amodel,LPTable Theta)
{
	int i;
	const int TT = Theta->Total;
	double r[TT];
	double M = 0.01;
	
	//vdRngGaussian(METHOD,stream,Theta->Total,r,0,M);
	for(i=0;i<TT;i++)
		r[i] = gsl_ran_gaussian(stream, M);

	Theta->Data[0] = 0;
	for(i=1;i<TT;i++)
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
	int i;
		
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

void CreateShape(int** VarSets,int* lenVarSets,int nVarSets,
				 int lenC,int* C,LPTable S,LPNTable shape)
{
	int i,k;
	const int SD = S->nDimens;
	int indexC[SD];
	int lenF;
	int lenD;
	double s;
	int okay;
	
	for(i=0;i<SD;i++)
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
		if(subset(SD,S->Index,indexC))
		{
			lenF = 0;
			for(k=0;k<SD;k++)
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
		if(subset(SD,VarSets[i],indexC))
		{
			s = 0.0;
			S->GetFirst();
			okay = 1;
			while(okay)
			{
				if(subset(SD,S->Index,indexC))
				{
					if(subset(SD,VarSets[i],S->Index))
					{
						lenD = lenVarSets[i];
						lenF = 0;
						for(k=0;k<SD;k++)
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
	int i;
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

void getNextTheta(int** VarSets,int* lenVarSets,int nVarSets,
                  int* amodel,int lenCL,int* CL,LPNTable smalltheta,LPTable Theta,LPTable nextTheta)
{
	int i,iF;
	int okay1;
	int iE, iL, iC;
	int* CLcomplement = new int[nextTheta->nDimens];
	int* CLfull = new int[nextTheta->nDimens];
	int* FunionL = new int[nextTheta->nDimens];
	const int TT = Theta->Total;
	double g[TT];
	int len[TT];
	
	for(i=0;i<nextTheta->nDimens;i++)
	{
		CLcomplement[i] = 1;
		CLfull[i] = 0;
	}
	for(i=0;i<lenCL;i++)
	{
		CLcomplement[CL[i]] = 0;
		CLfull[CL[i]] = 1;
	}
	
	for(i=0;i<nextTheta->Total;i++)
	{
		nextTheta->Data[i] = Theta->Data[i];
		g[i] = 0;
		len[i]=0;
	}
	
	////////////////////////////////////////////////
	Theta->GetFirst();
	okay1 = 1;
	iF = 0;
	while(okay1)
	{
		len[iF] = lenVarSets[iF];
										
		double s1 = 1.0;
		for(iL=1;iL<nVarSets;iL++)
		{	
			if(subset(Theta->nDimens,VarSets[iL],CLcomplement))
			{
				for(i=0;i<Theta->nDimens;i++)
				{
					FunionL[i] = 0;
					if(VarSets[iL][i]==1) FunionL[i] = 1;
					if(Theta->Index[i]==1) FunionL[i] = 1;
				}
				
				int thereisone = 0;
				double sC = 0.0;
				for(iC=1;iC<nVarSets;iC++)
				{
					if(0==subset(Theta->nDimens,VarSets[iC],Theta->Index))
					{
						if(1==subset(Theta->nDimens,VarSets[iC],FunionL))
						{
							sC += Theta->GetI(VarSets[iC]);
							thereisone = 1;
						}
					}
				}
				if(thereisone)
				{
					s1 += exp(sC);
				}
			}
		}
		g[iF]=log(s1);
		iF++;
		okay1 = Theta->GetNext();
	}
	////////////////////////////////////////////////
	
	nextTheta->GetFirst();
	nextTheta->Set(0);
	iE = 0;
	while(nextTheta->GetNext())
	{
	    iE++;
		if(amodel[iE])
		{
		    if(1==subset(Theta->nDimens,nextTheta->Index,CLfull))
			{
				double s0 = smalltheta->GetI(nextTheta->Index);
				int lenE = 0;
				for(i=0;i<nextTheta->nDimens;i++) lenE+=nextTheta->Index[i];
		
				Theta->GetFirst();
				int okay1 = 1;
				iF = 0;
				while(okay1)
				{
					if(len[iF]>=0)
					{
						if(subset(Theta->nDimens,Theta->Index,nextTheta->Index))
						{
							s0 += pow(-1,lenE-len[iF]-1)*g[iF];
						}
					}
					iF++;
					okay1 = Theta->GetNext();
				}
				nextTheta->Set(s0);
			}
			else
			{
				Theta->SetIndex(nextTheta->Index);
				nextTheta->Set(Theta->Get());
			}
		}
		else
		{
			nextTheta->Set(0);
		}
	}

	delete[] CLcomplement; CLcomplement = NULL;
	delete[] CLfull; CLfull = NULL;
	delete[] FunionL; FunionL = NULL;
	return;
}

void sampleNewThetas(int** VarSetsMarg,int* lenVarSetsMarg,int nVarSetsMarg,
                     int* amodel,int* aModelGenerators,LPTable mytable,LPTable Y,LPTable Theta,LPTable nextTheta)
{
	int i,j,k;
	double r,s;
		
	int nMssShape = 0;
	for(i=0;i<nVarSetsMarg;i++)
	{
		nMssShape += aModelGenerators[i];
	}

	int* lenC = new int[nMssShape];
	int** C = new int*[nMssShape];
	for(i=0;i<nMssShape;i++)
	{
		C[i] = new int[mytable->nDimens];
	}

	LPNTable PriorMssShape = new NTable[nMssShape];
	LPNTable tau = new NTable[nMssShape];
	LPNTable theta = new NTable[nMssShape];
	
	nMssShape = 0;
	for(i=0;i<nVarSetsMarg;i++)
	{
		if(aModelGenerators[i])
		{
			k = 0;
			for(j=0;j<mytable->nDimens;j++)
			{
				if(VarSetsMarg[i][j])
				{
					C[nMssShape][k]=j;
					k++;
				}
			}
			lenC[nMssShape] = k;
		
			tau[nMssShape].Create(lenC[nMssShape],C[nMssShape],mytable);
			theta[nMssShape].Create(lenC[nMssShape],C[nMssShape],mytable);
			PriorMssShape[nMssShape].Create(lenC[nMssShape],C[nMssShape],mytable);
			CreateShape(VarSetsMarg,lenVarSetsMarg,nVarSetsMarg,
						lenC[nMssShape],C[nMssShape],Y,&PriorMssShape[nMssShape]);
			nMssShape++;
		}
	}
	
	double beta = 1.0/Y->Data[0];
	for(k=0;k<nMssShape;k++)
	{
		s = 0;
		for(j=0;j<tau[k].Total;j++)
		{
			r = -1;
			while(r<=0)
			{
				//vdRngGamma(METHOD,stream,1,&r,PriorMssShape[k].Data[j],0,beta);
				r = gsl_ran_gamma(stream, PriorMssShape[k].Data[j], beta);
			}
			tau[k].Data[j] = r;
			s += r;
		}
		for(j=0;j<tau[k].Total;j++)
		{
			tau[k].Data[j] /= s;
		}
			
		getThetaMarginal(&theta[k],&tau[k]);
		getNextTheta(VarSetsMarg,lenVarSetsMarg,nVarSetsMarg,
					 amodel,lenC[k],C[k],&theta[k],Theta,nextTheta);
		for(i=0;i<Theta->Total;i++)
		{
			Theta->Data[i] = nextTheta->Data[i];
		}
	}	
	
	//clean memory
	for(i=0;i<nMssShape;i++)
	{
		PriorMssShape[i].Reset();
		tau[i].Reset();		
		theta[i].Reset();
	}
	delete[] PriorMssShape; PriorMssShape = NULL;
	delete[] tau; tau = NULL;
	delete[] theta; theta = NULL;
	for(i=0;i<nMssShape;i++)
	{
		delete[] C[i]; C[i] = NULL;
	}
	delete[] C; C = NULL;
	delete[] lenC; lenC = NULL;

	return;
}

void ipf(int** VarSetsMarg,int* lenVarSetsMarg,int nVarSetsMarg,
         LPTable P,LPTable dataTable,int* aModelGenerators,int* aModel)
{
	int i,j,k;
	double delta = 0.0000001;
	const int PT = P->Total;
	double m[PT];
	double mold[PT];
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
	InitThetaTable(aModel,Theta);
	PfromTheta(P,Theta);
	
	int notdone = 1;
	while(notdone)
	{
		for(j=0;j<PT;j++)
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
			
			for(j=0;j<PT;j++)
			{
				P->Data[j] = m[j]; 
			}
			
			p->Reset(); delete p ; p = NULL;

			s = 0.0;
			for(j=0;j<PT;j++)
			{
				s += P->Data[j];
			}
	
			for(j=0;j<PT;j++)
			{
				P->Data[j] /= s;
			}

			//transform in the space of thetas and impose the constraints
			ThetafromP(Theta,P,aModel);
			PfromTheta(P,Theta);
		}
	
		notdone = 0;
		for(i=0;i<PT;i++)
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
	for(i=0;i<PT;i++)
	{
		s += P->Data[i];
	}
	
	for(i=0;i<PT;i++)
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

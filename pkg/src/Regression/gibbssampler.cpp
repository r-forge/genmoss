/*
AUTHOR: ADRIAN DOBRA (ADOBRA@STAT.WASHINGTON.EDU)
*/
// gibbssampler.cpp: implementation of the CGibbsSampler class.
//
//////////////////////////////////////////////////////////////////////
#include "gibbssampler.h"

#include "covariance.cpp"

extern int subset(int n,int* set1,int* set2);
extern void MakeMarginals(LPTable Y,LPTable dataTable);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
CGibbsSampler::CGibbsSampler()
{
	//vslNewStream(&stream,BRNG,SEED);

	stream = gsl_rng_alloc(BRNG);
	gsl_rng_set(stream, SEED);

	tauPrior = 1;
	myPI = 3.1415926535897931159979635;
	nMaxRegressors = 3; //this gives the maximum number of accepted regressors
	//nMaxRegressorsInteraction = 5;
	//nDoInteractions = 0;
	//nDoSpaceRatio = 0;
	//nDoMetropolisH = 0;
	mdProbMax = 0.1;
	return;
}

CGibbsSampler::~CGibbsSampler()
{
	Cleanup();
	return;
}

void CGibbsSampler::Cleanup()
{
	//vslDeleteStream(&stream);
	gsl_rng_free(stream);
	return;
}

void CGibbsSampler::FirstIndex(int n,int* index,set<int>& VarInModel)
{
	int i;
	for(i=0;i<n;i++)
	{
		index[i]=0;
		VarInModel.erase(i);
	}
	return;
}

int CGibbsSampler::GetNextIndex(int n,int* index,set<int>& VarInModel)
{
	int i,j;
	
	for(i=0;i<n;i++)
	{
		if(0==index[i])
		{
			index[i]=1;
			VarInModel.insert(i);
			return(1);
		}
		else
		{
			VarInModel.erase(i);
			index[i]=0;
		}
	}
	return(0);
}

int CGibbsSampler::MakeInteractionIndex(int p,int i,int j)
{
	int k;
	
	//make sure i < j
	if(i>j)
	{
		k = i;
		i = j;
		j = k;
	}
	//the index of each variable starts at zero
	i = i+1;
	j = j+1;
	return((int)(p+(i-1)*(p-(((double)i)/2.0))+j));
}

void CGibbsSampler::GetInteractionIndex(int ind,int p,int& vi,int& vj)
{
	int i,j;
    if(ind<=p)
	{
		vi = -1;
		vj = ind;
		return;
	}
	ind = ind-p;
	for(i=1;i<=p;i++)
	{
		j = (int) (ind-(i-1)*(p-(((double)i)/2.0)));
		if((i<=j)&&(j<=p))
		{
			vi = i-1;
			vj = j-1;
			return;
		}
	}
	fprintf(stderr,"Fatal error :: GetInteractionIndex\n");
	exit(10);
	return;
} 

double* CGibbsSampler::NormalizeWeights(double* w,int nmax)
{
	int i;
	double* cumw = new double[nmax+1];

	//the weights are in increasing order
	//divide by the largest one, i.e. substract its log
	for(i=0;i<nmax;i++)
	{
		w[i] -= w[nmax-1];
	}
 
	for(i=0;i<nmax;i++)
	{
		w[i] = AdjustDouble(exp(w[i]));
	}

    cumw[0] = 0;
	for(i=1;i<=nmax;i++)
	{
		cumw[i] = cumw[i-1]+w[i-1];
    }
    for(i=1;i<=nmax;i++)
	{
		cumw[i] /= cumw[nmax];
    }
    return cumw;
}

int CGibbsSampler::WeightedSampling(int maxk,double* weights)
{
	double s;
	int i;
	int q1 = 0;
	int q2 = maxk;

	//vdRngUniform(METHOD,stream,1,&s,0,1);
	s = gsl_rng_uniform(stream);

	while(q1+1!=q2)
	{
		int q = (q1+q2)/2;
		if(weights[q]<s)
		{
			q1 = q;
		}
		else
		{
			q2 = q;
		}
	}
	return(q1);  
}

int CGibbsSampler::SelectEquation(int NumberOfGenes,CEquations* eq)
{
	int i,j;
	int equation;
	int chosenReg = -1;
	double delta;

	double* weights = new double[NumberOfGenes];
	int* ids = new int[NumberOfGenes];

	int nmax = eq->harvest(ids,weights);
	if(0==nmax)
	{
		return(-1);
	}
	
	//normalize the weights
	int chosenNeighbor = 0;
	double* cumw = NormalizeWeights(weights,nmax);
	chosenNeighbor = WeightedSampling(nmax,cumw);
	delete[] cumw; cumw = NULL;

	chosenReg = ids[chosenNeighbor];

	//clean memory
	delete[] weights; weights = NULL;
	delete[] ids; ids = NULL;
	return(chosenReg);
}

double** CGibbsSampler::gmatrix(double** beta,int p,CData& Data)
{
	int i,j,k;
	int n = Data.SampleSize;
	double s;
	
	double theta[n];
	for(i=0;i<n;i++)
	{
		s = 0;
		for(j=0;j<p;j++)
		{
			s += Data.X[i][j]*beta[j][0];
		}
		theta[i] = 1/(1+exp(-s));
	}
	
	double** G = allocmatrix(p,p);
	for(j=0;j<p;j++)
	{
		for(k=j;k<p;k++)
		{
			s = 0;
			for(i=0;i<n;i++)
			{
				s -= (1-theta[i])*theta[i]*Data.X[i][j]*Data.X[i][k];
			}
			if(k==j)
			{
				G[j][j] = s-(1.0/tauPrior); 
			}
			else
			{
				G[j][k] = s;
				G[k][j] = s;
			}
		}
	}
	return(G);
}

double** CGibbsSampler::gvector(double** beta,int p,CData& Data)
{
	int i,j;
	int n = Data.SampleSize;
	double s;
	
	double theta[n];
	for(i=0;i<n;i++)
	{
		s = 0;
		for(j=0;j<p;j++)
		{
			s += Data.X[i][j]*beta[j][0];
		}
		theta[i] = 1/(1+exp(-s));
	}
	
	double** g = allocmatrix(p,1);
	for(j=0;j<p;j++)
	{
		g[j][0] = -beta[j][0]/tauPrior;
		for(i=0;i<n;i++)
		{
			g[j][0] += (Data.Y[i]-theta[i])*Data.X[i][j];
		}
	}
	
	return(g);
}

double CGibbsSampler::loglikprior(double** beta,int p,CData& Data)
{
	int i,j;
	int n = Data.SampleSize;
	double s;
	double h = 0;
	double epsilon = 0.00001;
	
	for(i=0;i<p;i++)
	{
		h += pow(beta[i][0],2);
	}
	h = -p*log(2*myPI*tauPrior)/2-h/(2*tauPrior);
	
	double theta[n];
	for(i=0;i<n;i++)
	{
		s = 0;
		for(j=0;j<p;j++)
		{
			s += Data.X[i][j]*beta[j][0];
		}
		theta[i] = 1/(1+exp(-s));
		if(theta[i]<epsilon)
		{
			theta[i] = epsilon;
		}
		if(theta[i]>1-epsilon)
		{
			theta[i] = 1-epsilon;
		}
	}
	
	for(i=0;i<n;i++)
	{
		h += Data.Y[i]*log(theta[i])+(1-Data.Y[i])*log(1-theta[i]);
	}
	
	return(h);
}

double CGibbsSampler::logitnr(int p, CData& Data)
{
	int i,j;
	double s;
	double** betaSigma = NULL;
	double epsilon = 0.0001;
	double** A = allocmatrix(p,1);
	double** betaMu = allocmatrix(p,1);
	double** newbetaMu = allocmatrix(p,1);
	double oldvalue = loglikprior(newbetaMu,p,Data);
	double newvalue = 0;
	double* M = new double[p*p];
	
	while(1)
	{
		if(NULL!=betaSigma)
		{
			freematrix(p,betaSigma);
		}
		betaSigma = gmatrix(betaMu,p,Data);
		//inverse(p,betaSigma);
		memset(M,0,p*p*sizeof(double));
		for(i=0;i<p;i++)
		{
			for(j=0;j<p;j++)
			{
				M[i*p+j] = betaSigma[i][j];
			}
		}
		if(!Solve(p,M))
		{
			return(DBL_MAX);
		}
		for(i=0;i<p;i++)
		{
			for(j=0;j<p;j++)
			{
				betaSigma[i][j] = M[i*p+j];
			}
		} 
		
		double** g = gvector(betaMu,p,Data);
		matrixproduct(p,p,1,betaSigma,g,A);
		freematrix(p,g);
		for(i=0;i<p;i++)
		{
			newbetaMu[i][0] = betaMu[i][0]-A[i][0];
			//printf("newbetaMu[%d][0] = %.5lf\n",i+1,newbetaMu[i][0]);
		}
		newvalue = loglikprior(newbetaMu,p,Data);
		
		//printf("new value = %.5lf :: oldvalue = %.5lf\n",newvalue,oldvalue);
		
		if(fabs(oldvalue-newvalue)<=epsilon)
		{
			s = 0;
			for(i=0;i<p;i++)
			{
				s += fabs(newbetaMu[i][0]-betaMu[i][0]);
			}
			if(s<=epsilon)
			{
				break;
			}
		}
		for(i=0;i<p;i++)
		{
			betaMu[i][0] = newbetaMu[i][0];
		}
		oldvalue = newvalue;
	}
	for(i=0;i<p;i++)
	{
		betaMu[i][0] = newbetaMu[i][0];
	}
	oldvalue = newvalue;
	
	memset(M,0,p*p*sizeof(double));
	for(i=0;i<p;i++)
	{
		for(j=0;j<p;j++)
		{
			M[i*p+j] = betaSigma[i][j];
		}
	}
	double logDetM = getlogdet(p,M);
	
	double marglik = (p*log(2*myPI)/2.0)+(logDetM/2.0)+oldvalue;
	//clean memory
	delete[] M; M = NULL;
	freematrix(p,A);
	freematrix(p,betaSigma);
	freematrix(p,betaMu);
	freematrix(p,newbetaMu);
	return(marglik);
}

void CreateShape(LPTable S,LPTable shape)
{
	int i,j,k;
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

LPTable getMarginal(set<int>& aclique,CData& Discrete)
{
   int i,j;

   LPTable mytable = new Table;
   int nvars = aclique.size();
   int* index = new int[nvars];
   for(j=0;j<nvars;j++)
   {
      index[j] = 2;
   }
   if(!mytable->Alloc(index,nvars))
   {
      printf("Error allocating memory!\n");
      exit(1);
   }
   int* vars = new int[nvars];
   i = 0;
   for(set<int>::iterator it=aclique.begin();it!=aclique.end();it++)
   {
      vars[i] = *it;
	  i++;
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

LPTable getPrior(int nvars)
{
	int i,j;
	
	LPTable mytable = new Table;
	int* index = new int[nvars];
    for(j=0;j<nvars;j++)
    {
	   index[j] = 2;
    }
    if(!mytable->Alloc(index,nvars))
    {
	   printf("Error allocating memory!\n");
       exit(1);
    }
	
	double s = 1.0/((double)mytable->Total);	
	for(i=0;i<mytable->Total;i++)
	{
		mytable->Data[i] = s;
	}
	delete[] index; index = NULL;
	return(mytable);
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

double getScore(set<int>& aclique,CData& Discrete)
{
	double mypost = 0;
	
	if(0==aclique.size())
	{
		return(lgamma((double)(1+Discrete.SampleSize)));
	}
	
	LPTable dataTable = getMarginal(aclique,Discrete);
	LPTable priorTable = getPrior(aclique.size());
	LPTable datapriorTable = new Table;
	if(!datapriorTable->Alloc(dataTable->Dimens,dataTable->nDimens))
	{
		printf("Failed to allocate data+prior table.\n"); exit(1);
	}
	for(int i=0;i<datapriorTable->Total;i++)
	{
		datapriorTable->Data[i] = dataTable->Data[i]+priorTable->Data[i];
	}
	
	mypost = ExactPost(datapriorTable)-ExactPost(priorTable);
	
	datapriorTable->Reset(); delete datapriorTable; datapriorTable = NULL;
	dataTable->Reset(); delete dataTable; dataTable = NULL;
	priorTable->Reset(); delete priorTable; priorTable = NULL;
	return(mypost);
}

double CGibbsSampler::calculateLogPost(int TargetGene, set<int>& ExplanatoryGenes, CData& Data)
{
	int i,j;
	double logpost = 0.0;

	set<int> cliqueY;
	cliqueY.insert(TargetGene);

	set<int> cliqueX = ExplanatoryGenes;

        // add the confounding variables into X and YX cliques.
        j = Data.NumOfConfoundingVars;
        for(i=0; i<j; i++)
        {
                cliqueX.insert(TargetGene-j+i);
        }

	set<int> cliqueYX = cliqueX;
	cliqueYX.insert(TargetGene);
	
	logpost = getScore(cliqueYX,Data)-getScore(cliqueY,Data)-getScore(cliqueX,Data);
		
	return(logpost);
}

int CGibbsSampler::ChooseNextGene(vector<int>& PossiblePredictors)
{
	int index;
	if(0==PossiblePredictors.size())
	{
		return(-1);
	}
	//viRngUniform(METHOD,stream,1,&index,0,PossiblePredictors.size());
	index = gsl_rng_uniform_int(stream, PossiblePredictors.size());

	int selectedGene = PossiblePredictors[index];
	PossiblePredictors.erase(PossiblePredictors.begin() + index);
	return selectedGene;
}

int CGibbsSampler::ChooseNextGeneSet(set<int>& PossiblePredictors)
{
	int i;
	int index;
	if(0==PossiblePredictors.size())
	{
		return(-1);
	}
	//viRngUniform(METHOD,stream,1,&index,0,PossiblePredictors.size());
	index = (int)gsl_rng_uniform_int(stream, PossiblePredictors.size());

	set<int>::iterator it;
	i = 0;
	for(it=PossiblePredictors.begin();it!=PossiblePredictors.end();it++)
	{
		if(i==index)
		{
			return(*it);
		}
		i++;
	}
	return(-1);
}

set<int> CGibbsSampler::RandomRegression(set<int> Predictors)
{
	int i,j;
	int index;
	set<int> reg;
	int reglen;
	int maxlen = nMaxRegressors;
	
	if(maxlen>Predictors.size())
	{
		maxlen = Predictors.size();
	}
	
	if(0==maxlen) return(reg);
	//viRngUniform(METHOD,stream,1,&reglen,0,maxlen);
	reglen = (int) gsl_rng_uniform_int(stream, maxlen);

	for(i=0;i<reglen;i++)
	{
		//viRngUniform(METHOD,stream,1,&index,0,Predictors.size());
		index = (int) gsl_rng_uniform_int(stream, Predictors.size());
		set<int>::iterator it;
		j = 0;
		for(it=Predictors.begin();it!=Predictors.end();it++)
		{
			if(j==index)
			{
				reg.insert(*it);
				Predictors.erase(*it);
				break;
			}
			j++;
		}
	}
	
	return(reg);
}

double CGibbsSampler::metropolisReg(int targetGene,CData& Data,CRegression* hybridreglist)
{
		int i;
	double r;
	set<int> VarInModel;
	double currentMarginalLik;
	double nextMargLik;
	double ratio;
	
	set<int> allPredictors;
    for(i=0;i<Data.NumberOfGenes;i++)
	{
		if(i!=targetGene)
		{
			allPredictors.insert(i);
		}
	}
	
	//generate a random starting model
	set<int> Predictors = allPredictors;
	VarInModel = RandomRegression(Predictors);
	currentMarginalLik = calculateLogPost(targetGene, VarInModel, Data);
	
	//add the fake variables
	int fakevar = Data.NumberOfGenes+1;
	for(i=VarInModel.size();i<nMaxRegressors;i++)
	{
		VarInModel.insert(fakevar);
		fakevar++;
	}
	
	set<int> VarOutsideModel = allPredictors;
	for(set<int>::iterator it=VarInModel.begin();it!=VarInModel.end();it++)
	{
		int anothergene = *it;
		VarOutsideModel.erase(anothergene);
		VarOutsideModel.insert(fakevar);
		fakevar++;
	}

	double spaceratio = 0;
	double totaliterations = 0;
	for(int iteration=0;iteration<nChainIterations;iteration++)
	{
		//pick a variable outside the model
		set<int> shotgunPredictors = VarOutsideModel;
		int geneToAdd = ChooseNextGeneSet(shotgunPredictors);
		if(geneToAdd<0)
		{
			fprintf(stderr,"There are no variables outside the model. ERROR!\n");
			exit(1);
		}
		
		//pick a variable currently in the model
		shotgunPredictors = VarInModel;
		int geneToDelete = ChooseNextGeneSet(shotgunPredictors);
		if(geneToDelete<0)
		{
			fprintf(stderr,"There are no variables in the model. ERROR!\n");
			exit(1);
		}
		
		//perform the exchange
		set<int> newVarOutsideModel = VarOutsideModel;
		newVarOutsideModel.erase(geneToAdd);
		newVarOutsideModel.insert(geneToDelete);
		
		set<int> newVarInModel = VarInModel;
		newVarInModel.erase(geneToDelete);
		newVarInModel.insert(geneToAdd);

		set<int> cVarInModel = newVarInModel;
		for(set<int>::iterator it=newVarInModel.begin();it!=newVarInModel.end();it++)
		{
			int anothergene = *it;
			if(anothergene>Data.NumberOfGenes)
			{
				cVarInModel.erase(anothergene);
			}
		}
	
		nextMargLik = calculateLogPost(targetGene, cVarInModel, Data);
		ratio = exp(nextMargLik-currentMarginalLik);
		if(ratio>1)
		{
			ratio = 1;
		}

		double u;
		//vdRngUniform(METHOD,stream,1,&u,0,1);
		u = gsl_rng_uniform(stream);
 
		if(u<ratio)
		{	
			//move is accepted
			currentMarginalLik = nextMargLik;
			VarInModel = newVarInModel;
			VarOutsideModel = newVarOutsideModel;
		}
		
		cVarInModel = VarInModel;
		for(set<int>::iterator it=VarInModel.begin();it!=VarInModel.end();it++)
		{
			int anothergene = *it;
			if(anothergene>Data.NumberOfGenes)
			{
				cVarInModel.erase(anothergene);
			}
		}
		if(hybridreglist->existReg(currentMarginalLik,cVarInModel))
		{
			spaceratio += 1;
		}
		totaliterations += 1;
		
		printf("spaceratio = %.5lf\n",spaceratio);
	}
	
	spaceratio = spaceratio/totaliterations;
    return(spaceratio);
}

void CGibbsSampler::exhaustiveReg(int targetGene,CData& Data,double* probs,FILE* modelsfile)
{
	int i;
	set<int> VarInModel;
	double currentMarginalLik;
	int* index = new int[Data.NumberOfGenes];
	CRegression* reglist = new CRegression;
	
	int model = 0;
	FirstIndex(Data.NumberOfGenes,index,VarInModel);
	do
	{
		if(VarInModel.find(targetGene) == VarInModel.end())
		{
			model++;
			if(model%10000==0)
			{
				printf("Processing regression [%d] :: [%d]\n",model,VarInModel.size());
				reglist->Delete();
			}
			currentMarginalLik = calculateLogPost(targetGene, VarInModel, Data);
			reglist->AddReg(currentMarginalLik,VarInModel);
		}		
	}while(GetNextIndex(Data.NumberOfGenes,index,VarInModel));
	
	reglist->NormalizeWeights();
    reglist->MakeCovProbs(Data.NumberOfGenes,probs);
	reglist->SaveAll(nMaxRegressors,modelsfile);
	reglist->DeleteAll();
		
	//clean memory
	delete[] index; index = NULL;
	delete reglist; reglist = NULL;
	return;
}

CRegression* CGibbsSampler::modelSelection(int targetGene,CData& Data,double* probs,FILE* modelsfile,FILE* logfile,int* countmodels)
{
	int i,j;
	double s, r;
	double maxProb[nChainReplicates];
	CRegression* reglists = new CRegression[nChainReplicates];
	CRegression* allregs = new CRegression;
	
	for(i=0;i<nChainReplicates;i++)
	{
		reglists[i].SetCutoffs(mdCutoffMax,mdCutoffMin);
		countmodels[i] = shotgunReg(targetGene,Data,&reglists[i],i+1,logfile);
	}
	allregs->SetCutoffs(mdCutoffMax,mdCutoffMax);
	for(i=0;i<nChainReplicates;i++)
	{
		allregs->AddList(&reglists[i]);
		reglists[i].DeleteAll();
	}
	delete[] reglists; reglists = NULL;

	//entireNeighborhood(targetGene,Data,allregs,0,logfile);
	
	//if(nDoSpaceRatio)
	//{
	//	double spr = metropolisReg(targetGene,Data,allregs);
	//	fprintf(spaceratiofile,"%d\t%d\t%.5lf\t%.5lf\n",Data.NumberOfGenes-1,allregs->NumberRegressions(),mdCutoffMax,spr);
	//}
				
	allregs->NormalizeWeights();
    allregs->MakeCovProbs(Data.NumberOfGenes,probs);
	allregs->SaveAll(nMaxRegressors,modelsfile);
	//allregs->DeleteAll();
	//delete allregs; allregs = NULL;
	return(allregs);
}

void CGibbsSampler::modelSelectionInteraction(int targetGene,CData& Data,double* probs,FILE* modelsfile,FILE* logfile,CRegression* allregs)
{
	int i,j;
	double s, r;
	CRegression* interactregs = new CRegression;
	interactregs->SetCutoffs(mdCutoffMax,mdCutoffMax);
	
	//printf("targetGene = %d :: Data.NumberOfGenes = %d\n",targetGene,Data.NumberOfGenes);
	
	for(i=0;i<nChainReplicates;i++)
	{
		j = 0;
		CRegression* p = allregs->Next;
		while(NULL!=p)
		{
			/*
			printf("processing regression with prob [%.5lf] :: %.5lf\n",p->Weight,calculateLogPost(targetGene,p->Vars,Data));
			for(set<int>::iterator it=(p->Vars).begin();it!=(p->Vars).end();it++)
			{
				j = *it;
				printf("\t%d",j+1);
			}
			printf("\n");
			*/
			s = interactReg(targetGene,p->Vars,Data,interactregs);
			fprintf(logfile,"%d\t%d\t%.5lf\n",i+1,j+1,s);
			printf("interactions %d\t%d\t%.5lf\n",i+1,j+1,s);
			p = p->Next;
			j++;
		}
	}
				
	interactregs->NormalizeWeights();
    interactregs->MakeCovProbs(Data.NumberOfGenes,probs);
	//interactregs->SaveAll(nMaxRegressorsInteraction,modelsfile);
	interactregs->DeleteAll();
	delete interactregs; interactregs = NULL;
	return;
}

CRegression* CGibbsSampler::modelSelectionMH(int targetGene,CData& Data,double* probs,FILE* modelsfile,FILE* logfile,int* countmodels)
{
	int i,j;
	double s, r;
	double maxProb[nChainReplicates];
	CRegression* reglists = new CRegression[nChainReplicates];
	CRegression* allregs = new CRegression;
	
	for(i=0;i<nChainReplicates;i++)
	{
		reglists[i].SetCutoffs(mdCutoffMax,mdCutoffMax);
		mhReg(targetGene,Data,&reglists[i],i+1,logfile,countmodels[i]);
	}
	allregs->SetCutoffs(mdCutoffMax,mdCutoffMax);
	for(i=0;i<nChainReplicates;i++)
	{
		allregs->AddList(&reglists[i]);
		reglists[i].DeleteAll();
	}
	delete[] reglists; reglists = NULL;

	allregs->NormalizeWeights();
    allregs->MakeCovProbs(Data.NumberOfGenes,probs);
	allregs->SaveAll(nMaxRegressors,modelsfile);
	return(allregs);
}

//given a set of main effects, identify good regressions with interactions
//returns the maximum posterior probability of a model identified
double CGibbsSampler::interactReg(int targetGene,set<int>& MainEffects,CData& Data,CRegression* reglist)
{
	double r;
	double currentMarginalLik;
	CRegression* localreglist = new CRegression;
	localreglist->SetCutoffs(mdCutoffMax,mdCutoffMin);
	
	currentMarginalLik = calculateLogPost(targetGene,MainEffects,Data);		
	localreglist->AddReg(currentMarginalLik,MainEffects);
	reglist->AddReg(currentMarginalLik,MainEffects);
	//printf("first time added [%.5lf]\n",currentMarginalLik);
	
	while(1)
	{
		int nmodels = 0;
		int ntotalmodels = 0;
		CRegression* nextmodel = localreglist->Shotgun(stream,nmodels,ntotalmodels);
		if(NULL==nextmodel)
		{
			break;
		}		
		set<int> VarInModel = nextmodel->Vars;
		
		//interaction terms
		set<int>::iterator itI;
		set<int>::iterator itJ;
		for(itI=MainEffects.begin();itI!=MainEffects.end();itI++)
		{
			int geneI = *itI;
			for(itJ=itI;itJ!=MainEffects.end();itJ++)
			{
				int geneJ = *itJ;
				int gene = MakeInteractionIndex(Data.NumberOfGenes-1,geneI,geneJ);
				//printf("geneI = %d :: geneJ = %d :: gene = %d\n",geneI,geneJ,gene);
								
				if(VarInModel.find(gene) == VarInModel.end())
				{
					//add the variable to the model
					VarInModel.insert(gene);
					//if(VarInModel.size()<=nMaxRegressorsInteraction)
					//{
					//	currentMarginalLik = calculateLogPost(targetGene, VarInModel, Data);
					//	localreglist->AddReg(currentMarginalLik,VarInModel);
					//	reglist->AddReg(currentMarginalLik,VarInModel);
					//}
					//delete the variable from the model
					VarInModel.erase(gene);
				}  
				else if(VarInModel.find(gene) != VarInModel.end())
				{
					//take the variable out of the model
					VarInModel.erase(gene);
					if(VarInModel.empty())
					{
						printf("Delete an interaction term and ended up with no repressors.\n");
						exit(1);
					}
					else
					{
						currentMarginalLik = calculateLogPost(targetGene, VarInModel, Data);
					}
					localreglist->AddReg(currentMarginalLik,VarInModel);
					reglist->AddReg(currentMarginalLik,VarInModel);
					//put the variable back in the model
					VarInModel.insert(gene);
				}
			}
		}
		
		//vdRngUniform(METHOD,stream,1,&r,0,1);
		r = gsl_rng_uniform(stream);

		if(r<=mdProbMax)
		{
			localreglist->Delete();
		}
	}
	
	currentMarginalLik = localreglist->getMaxProb();
	localreglist->DeleteAll();
	delete localreglist; localreglist = NULL;
	return(currentMarginalLik);
}

int CGibbsSampler::shotgunReg(int targetGene,CData& Data,CRegression* reglist,int astartpoint,FILE* logfile)
{
	int i, k;
	double r;
	set<int> VarInModel;
	double currentMarginalLik;
	
	set<int> allPredictors;
	// k is the number of SNPs in the Data:
        // = total number of cols - confounding vars - target variable
        k = Data.NumberOfGenes - Data.NumOfConfoundingVars - 1;

    	for(i=0;i<k;i++)
	{
		//if(i!=targetGene)
		//{
			allPredictors.insert(i);
		//}
	}
	
	int nevalmodels = 0;
	while(1)
	{
		reglist->nIterationNumber++;
		int nmodels = 0;
		int ntotalmodels = 0;
		CRegression* nextmodel = reglist->Shotgun(stream,nmodels,ntotalmodels);
			
		if(NULL==nextmodel)
		{
			if(0==ntotalmodels)
			{
				//generate a random starting model
				set<int> Predictors = allPredictors;
				VarInModel = RandomRegression(Predictors);
				currentMarginalLik = calculateLogPost(targetGene, VarInModel, Data);
				reglist->AddReg(currentMarginalLik,VarInModel);
				nextmodel = reglist->Shotgun(stream,nmodels,ntotalmodels);
			}
			else
			{
				return(nevalmodels);
			}
		}

		fprintf(logfile,"%d\t%d\t%.5lf\t%.5lf\t%d\t%d\n",astartpoint,reglist->nIterationNumber,reglist->mdCutoffMax,reglist->Next->Weight,ntotalmodels,nmodels);
		printf("%d\t%d\t%.5lf\t%.5lf\t%d\t%d\n",astartpoint,reglist->nIterationNumber,reglist->mdCutoffMax,reglist->Next->Weight,ntotalmodels,nmodels);
	
		VarInModel = nextmodel->Vars;			
		//consider all the neighbors of the current model
		//and add them to the list of models
		set<int> shotgunPredictors = allPredictors;
		while(!shotgunPredictors.empty())
		{
			int gene = ChooseNextGeneSet(shotgunPredictors);
			if(gene<0) break; //just in case
			shotgunPredictors.erase(gene);
				
			//test if the variable is currently in the model
			if(VarInModel.find(gene) == VarInModel.end())
			{
				//add the variable to the model
				VarInModel.insert(gene);
				if(VarInModel.size()<=nMaxRegressors)
				{
					nevalmodels++;
					currentMarginalLik = calculateLogPost(targetGene, VarInModel, Data);
					reglist->AddReg(currentMarginalLik,VarInModel);
				}
							
				//substitute each variable currently in the model with this variable
				for(set<int>::iterator it=VarInModel.begin();it!=VarInModel.end();it++)
				{
					int anothergene = *it;
					if(gene!=anothergene)
					{
						set<int> cVarInModel = VarInModel;
						cVarInModel.erase(anothergene);
						nevalmodels++;
						currentMarginalLik = calculateLogPost(targetGene, cVarInModel, Data);
						if(cVarInModel.size()>nMaxRegressors)
						{
							printf("Found regression of size [%d] greater than [%d]\n",cVarInModel.size(),nMaxRegressors);
							exit(1);
						}
						reglist->AddReg(currentMarginalLik,cVarInModel);
					}
				}
				//delete the variable from the model
				VarInModel.erase(gene); 
			}
			else if(VarInModel.find(gene) != VarInModel.end())
			{
				//take the variable out of the model
				VarInModel.erase(gene);
				nevalmodels++;
				currentMarginalLik = calculateLogPost(targetGene, VarInModel, Data);
				reglist->AddReg(currentMarginalLik,VarInModel);
				//put the variable back in the model
				VarInModel.insert(gene);
			}
		}
		
		//vdRngUniform(METHOD,stream,1,&r,0,1);
		r = gsl_rng_uniform(stream);

		if(r<=mdProbMax)
		{
			reglist->Delete();
		}
	}
	return(nevalmodels);
}

void CGibbsSampler::mhReg(int targetGene,CData& Data,CRegression* reglist,int astartpoint,FILE* logfile,int nMhIterations)
{
	int i;
	double r;
	set<int> VarInModel;
	double currentMarginalLik;
	double nextMargLik;
	double ratio;
	
	set<int> allPredictors;
    for(i=0;i<Data.NumberOfGenes;i++)
	{
		if(i!=targetGene)
		{
			allPredictors.insert(i);
		}
	}
	
	//generate a random starting model
	set<int> Predictors = allPredictors;
	VarInModel = RandomRegression(Predictors);
	currentMarginalLik = calculateLogPost(targetGene, VarInModel, Data);
	reglist->AddReg(currentMarginalLik,VarInModel);
	
	//add the fake variables
	int fakevar = Data.NumberOfGenes+1;
	for(i=VarInModel.size();i<nMaxRegressors;i++)
	{
		VarInModel.insert(fakevar);
		fakevar++;
	}
	
	set<int> VarOutsideModel = allPredictors;
	for(set<int>::iterator it=VarInModel.begin();it!=VarInModel.end();it++)
	{
		int anothergene = *it;
		VarOutsideModel.erase(anothergene);
		VarOutsideModel.insert(fakevar);
		fakevar++;
	}
	
	for(int iteration=0;iteration<nMhIterations;iteration++)
	{
		int ntotalmodels = reglist->NumberRegressions();
		if(iteration%1000==0)
		{
			fprintf(logfile,"%d\t%d\t%.5lf\t%.5lf\t%d\t%d\n",astartpoint,iteration+1,reglist->mdCutoffMax,reglist->Next->Weight,ntotalmodels,0);
			printf("%d\t%d\t%.5lf\t%.5lf\t%d\t%d\n",astartpoint,iteration+1,reglist->mdCutoffMax,reglist->Next->Weight,ntotalmodels,0);
		}

		//pick a variable outside the model
		set<int> shotgunPredictors = VarOutsideModel;
		int geneToAdd = ChooseNextGeneSet(shotgunPredictors);
		if(geneToAdd<0)
		{
			fprintf(stderr,"There are no variables outside the model. ERROR!\n");
			exit(1);
		}
		
		//pick a variable currently in the model
		shotgunPredictors = VarInModel;
		int geneToDelete = ChooseNextGeneSet(shotgunPredictors);
		if(geneToDelete<0)
		{
			fprintf(stderr,"There are no variables in the model. ERROR!\n");
			exit(1);
		}
		
		//perform the exchange
		set<int> newVarOutsideModel = VarOutsideModel;
		newVarOutsideModel.erase(geneToAdd);
		newVarOutsideModel.insert(geneToDelete);
		
		set<int> newVarInModel = VarInModel;
		newVarInModel.erase(geneToDelete);
		newVarInModel.insert(geneToAdd);

		set<int> cVarInModel = newVarInModel;
		for(set<int>::iterator it=newVarInModel.begin();it!=newVarInModel.end();it++)
		{
			int anothergene = *it;
			if(anothergene>Data.NumberOfGenes)
			{
				cVarInModel.erase(anothergene);
			}
		}
		nextMargLik = calculateLogPost(targetGene, cVarInModel, Data);
		reglist->AddReg(nextMargLik,cVarInModel);
		
		ratio = exp(nextMargLik-currentMarginalLik);
		if(ratio>1)
		{
			ratio = 1;
		}

		double u;
		//vdRngUniform(METHOD,stream,1,&u,0,1);
		u = gsl_rng_uniform(stream);
 
		if(u<ratio)
		{	
			//move is accepted
			currentMarginalLik = nextMargLik;
			VarInModel = newVarInModel;
			VarOutsideModel = newVarOutsideModel;
		}
	}

	return;
}

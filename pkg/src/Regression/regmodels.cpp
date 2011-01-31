/*
AUTHOR: ADRIAN DOBRA (ADOBRA@STAT.WASHINGTON.EDU)
*/

#include "regmodels.h"

CRegression::CRegression()
{
	notStudied = 1;
	Next = NULL;
	Weight = -999999999999999.99;
	normWeight = 0;
	nmodels = 0;
	mdCutoffMax = 1.0;
	mdCutoffMin = 1.0;
	nIterationNumber = 0;
}

CRegression::~CRegression()
{
	return;
}

void CRegression::SetCutoffs(double vmax,double vmin)
{
	mdCutoffMax = vmax;
	mdCutoffMin = vmin;
	return;
}

void CRegression::ResetStudied()
{
	CRegression* p = this->Next;
	while(NULL!=p)
	{
		p->notStudied = 1;
		p = p->Next;
	}
	return;	
}

void CRegression::DeleteAll()
{
	CRegression* p = Next;
	CRegression* pnext;
	while(NULL!=p) {
		pnext = p->Next;
		p->Weight = 0;
		p->Next = NULL;
		delete p;
		p = pnext;
	}
	Next = NULL;
	return;
}

void CRegression::Delete()
{
	CRegression* p = Next;
	if(NULL==Next)
	{
		return;
	}
	double maxp = p->Weight;
	DeleteNotRelevant(maxp);
}

void CRegression::DeleteNotRelevant(double maxweight)
{
	CRegression* p = this;
	CRegression* pnext;
	
	while(NULL!=p->Next)
	{
		if(exp(p->Next->Weight-maxweight)<=mdCutoffMax)
		{
			pnext = p->Next;
			p->Next = pnext->Next;
			delete pnext;
		}
		else
		{
			p = p->Next;
		}
	}
	return;
}

int CRegression::existReg(double weight,set<int> vars)
{
	CRegression* p = this;
	
	if(NULL==p->Next)
	{
		return(0);
	}
	while(NULL!=p->Next)
	{
		p = p->Next;
		if((p->Weight==weight)&&(p->Vars==vars))
		{
			return(1);
		}
	}
	return(0);
}

double CRegression::getMaxProb()
{
	if(NULL==this->Next)
	{
		return(0);
	}
	return(this->Next->Weight);
}

void CRegression::AddList(CRegression* reglist)
{
	CRegression* p = reglist->Next;
	while(NULL!=p)
	{
		AddReg(p->Weight,p->Vars);
		p = p->Next;
	}
	//Delete();

	return;
}

int CRegression::NumberRegressions()
{
	CRegression* p = this->Next;
	int n = 0;
	while(NULL!=p)
	{
		n++;
		p = p->Next;
	}
	return(n);
}

int CRegression::AddReg(double weight,set<int> vars)
{
	CRegression* p = this;
	CRegression* pnext;
	int newModelAdded = 0;
	
	//printf("weight = %.5lf\n",weight);
	if(NULL==p->Next)
	{
		CRegression* newp = new CRegression;
		newp->Vars = vars;
		newp->Weight = weight;
		newp->Next = p->Next;
		p->Next = newp;
		nmodels++;
		newModelAdded = 1;
	}
	else
	{
		double maxweight = p->Next->Weight;
		if(weight>maxweight) //new best model is found
		{
			/*
			printf("new model %.5lf\n",weight);
			for(set<int>::iterator it=vars.begin();it!=vars.end();it++)
			{
				printf("\t%d",(*it)+1);
			}
			printf("\n");
			*/
			
			CRegression* newp = new CRegression;
			newp->Vars = vars;
			newp->Weight = weight;
			newp->Next = p->Next;
			p->Next = newp;
			nmodels++;
			newModelAdded = 1;
		}
		else if(exp(weight-maxweight)>=mdCutoffMin)
		{
			int found = 0;
			while(NULL!=p->Next)
			{
				if(p->Next->Weight>weight)
				{
					p = p->Next;
				}
				else if((p->Next->Weight==weight)&&(p->Next->Vars==vars))
				{
					found = 1;
					break;
				} else if(p->Next->Weight==weight)
				{
					p = p->Next;
				}
				else
				{
					break;
				}
			}
	
			if(!found)
			{
				CRegression* newp = new CRegression;
				newp->Vars = vars;
				newp->Weight = weight;
				newp->Next = p->Next;
				p->Next = newp;
				nmodels++;
				newModelAdded = 1;
			}
		}
	}
	return(newModelAdded);
}

double CRegression::NormalizeWeights()
{
	CRegression* p = this;
	if(NULL==p->Next)
	{
		return(0);
	}
	p = p->Next;
	double maxp = p->Weight;
	p->normWeight = 1;
	double sump = 1;
	while(NULL!=p->Next)
	{
		p = p->Next;
		p->normWeight = exp(p->Weight-maxp);
		sump += p->normWeight;
	}
	p = this;
	while(NULL!=p->Next)
	{
		p = p->Next;
		p->normWeight = p->normWeight/sump;
	}
	return(sump);
}

void CRegression::MakeCovProbs(int ncov,double* probs)
{
	int i;
	double s;
	memset(probs,0,ncov*sizeof(double));
	CRegression* p = this;
	while(NULL!=p->Next)
	{
		p = p->Next;
		set<int>::iterator it;
		for(it=(p->Vars).begin();it!=(p->Vars).end();it++)
		{
			if(*it<ncov)
			{
				probs[*it] += p->normWeight;
			}
		}
	}
	return;
}

void CRegression::LoadAll(int nMaxRegressors,FILE* file)
{
	int i,j;
	double weight,normweight;
	CRegression* p = this;
	while(2==fscanf(file,"%lf %lf",&weight,&normweight))
	{
		set<int> vars;
		for(i=0;i<nMaxRegressors;i++)
		{
			if(1!=fscanf(file,"%d",&j))
			{
				printf("This file does not seem to be a correct regression file.\n");
				exit(1);
			}
			if(j>=1)
			{
				vars.insert(j-1);
			}
		}
	
		CRegression* newp = new CRegression;
		newp->Vars = vars;
		newp->Weight = weight;
		newp->normWeight = normweight;
		newp->Next = p->Next;
		p->Next = newp;
		p = newp;
	}
	return;
}

void CRegression::SaveAll(int nMaxRegressors,FILE* file)
{
	int i,j;
	CRegression* p = this;
	while(NULL!=p->Next)
	{
		p = p->Next;
		fprintf(file,"%.5lf\t%.5lf",p->Weight,p->normWeight);
		set<int>::iterator it;
		i = 0;
		for(it=(p->Vars).begin();it!=(p->Vars).end();it++)
		{
			j = *it;
			fprintf(file,"\t%d",j+1);
			i++;
		}
		for(j=i;j<nMaxRegressors;j++)
		{
			fprintf(file,"\t%d",-1);
		}
		fprintf(file,"\n");
	}
	return;
}

double* CRegression::NormalizeWeightsShotgun(double* w,int nmax)
{
	int i;
	double* cumw = new double[nmax+1];

	double maxw = w[0];
	for(i=1;i<nmax;i++)
	{
		if(maxw<w[i])
		{
			maxw = w[i];
		}
	}
	
	for(i=0;i<nmax;i++)
	{
		w[i] -= maxw;
	}
 
	for(i=0;i<nmax;i++)
	{
		w[i] = exp(w[i]);
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

int CRegression::WeightedSamplingShotgun(int maxk,double* weights, gsl_rng * stream)
{
	double s;
	int i;
	int q1 = 0;
	int q2 = maxk;

	//vdRngUniform(0,stream,1,&s,0,1);
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


CRegression* CRegression::Shotgun(gsl_rng * stream,int& nPmodels,int& nPtotalmodels)
{
	int nmodels = 0;
	int ntotalmodels = 0;
	CRegression* p = this;
	
	if(NULL==p->Next)
	{
		return(NULL);
	}
	while(NULL!=p->Next)
	{
		if(p->Next->notStudied)
		{
			nmodels++;
		}
		ntotalmodels++;
		p = p->Next;
	}
	
	nPmodels = nmodels;
	nPtotalmodels = ntotalmodels;

	if(0==nmodels)
	{
		return(NULL);
	}
	double* w = new double[nmodels];
	p = this;
	nmodels = 0;
	while(NULL!=p->Next)
	{
		if(p->Next->notStudied)
		{
			w[nmodels] = p->Next->Weight;
			nmodels++;
		}
		p = p->Next;
	}
	
	double* cumw = NormalizeWeightsShotgun(w,nmodels);
	int chosenmodel = WeightedSamplingShotgun(nmodels,cumw,stream);
	delete[] cumw; cumw = NULL;
	delete[] w; w = NULL;
	
	p = this;
	nmodels = 0;
	while(NULL!=p->Next)
	{
		if(p->Next->notStudied)
		{
			if(nmodels==chosenmodel)
			{
				p->Next->notStudied = 0;
				return(p->Next);
			}
			nmodels++;
		}
		p = p->Next;
	}
	return(NULL);
}

CRegression* CRegression::FirstModel(gsl_rng * stream,int& nPmodels,int& nPtotalmodels)
{
	int nmodels = 0;
	int ntotalmodels = 0;
	CRegression* p = this;
	
	if(NULL==p->Next)
	{
		return(NULL);
	}
	while(NULL!=p->Next)
	{
		nmodels += p->Next->notStudied;
		ntotalmodels++;
		p = p->Next;
	}
	
	nPmodels = nmodels;
	nPtotalmodels = ntotalmodels;
	
	if(0==nmodels)
	{
		return(NULL);
	}

	p = this;
	while(NULL!=p->Next)
	{
		if(p->Next->notStudied)
		{
			p->Next->notStudied = 0;
			break;
		}
		p = p->Next;
	}
	return(p->Next);
}

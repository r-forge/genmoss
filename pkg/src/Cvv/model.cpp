/*
AUTHOR: ADRIAN DOBRA (ADOBRA@STAT.WASHINGTON.EDU)
*/

#include "model.h"

CModel::CModel()
{
	mdCutoffMax = 1.0;
	mdCutoffMin = 1.0;
	notStudied = 1;
	Counter = 0;
	Length = 0;
	vars = NULL;
	laplace = 0.0;
	Next = NULL;
}

CModel::~CModel()
{
	if(NULL!=vars)
	{
		delete[] vars; vars = NULL;
	}
	Length = 0;
	laplace = 0;
	Next = NULL;
}

void CModel::SetCutoffs(double vmax,double vmin)
{
	mdCutoffMax = vmax;
	mdCutoffMin = vmin;
	return;
}

CModel* CModel::Find(int len,int* amodel)
{
	int i;
	CModel* p = this;
	
	while(NULL!=p->Next)
	{
		int r = LexCompare(len,amodel,p->Next->vars);
		if(0==r)//found this model
		{
			return(p->Next);
		}
		else if(-1==r)
		{
			p = p->Next;
		}
		else //1==r
		{
			break;
		}
	}
	return(NULL);
}

void CModel::NormalizeWeights()
{
	CModel* p = this;
	
	if(NULL==p->Next)
	{
		return;
	}
	double maxlaplace = p->Next->laplace;
	p = p->Next;
	while(NULL!=p->Next)
	{
		if(maxlaplace<p->Next->laplace)
		{
			maxlaplace = p->Next->laplace;
		}
		p = p->Next;
	}
	
	double sumlaplace = 0.0;
	p = this;
	while(NULL!=p->Next)
	{
		p->Next->laplace = exp(p->Next->laplace-maxlaplace);
		sumlaplace += p->Next->laplace;
		p = p->Next;
	}
	p = this;
	while(NULL!=p->Next)
	{
		p->Next->laplace = p->Next->laplace/sumlaplace;
		p = p->Next;
	}
	
	return;
}

double* CModel::NormalizeWeightsShotgun(double* w,int nmax)
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

int CModel::WeightedSamplingShotgun(int maxk,double* weights, gsl_rng * stream)
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

CModel* CModel::Shotgun(gsl_rng * stream,FILE* logfile)
{
	int nmodels = 0;
	int ntotalmodels = 0;
	CModel* p = this;
	
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
	printf("Models in list [%d] not studied [%d]\n",ntotalmodels,nmodels);
	fprintf(logfile,"\t%d\t%d\n",ntotalmodels,nmodels);
	
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
			w[nmodels] = p->Next->laplace;
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

void CModel::Diagnostics(double& maxprob,int& nrelmodels)
{
	maxprob = 0;
	nrelmodels = 0;
	CModel* p = this;
	if(NULL==p->Next)
	{
		return;
	}
	maxprob = p->Next->laplace;
	p = p->Next;
	while(NULL!=p->Next)
	{
		if(maxprob<p->Next->laplace)
		{
			maxprob = p->Next->laplace;
		}
		p = p->Next;
	}
	p = this;
	while(NULL!=p->Next)
	{
		if(exp(p->Next->laplace-maxprob)>=mdCutoffMax)
		{
			nrelmodels++;
		}
		p = p->Next;
	}
	this->laplace = maxprob;

	return;
}


void CModel::MakeCovProbs(int ncovprob,double* covprob)
{
	int i;
	CModel* p = this;
	
	if(NULL==p->Next)
	{
		return;
	}
	memset(covprob,0,ncovprob*sizeof(double));
	while(NULL!=p->Next)
	{
		for(i=0;i<ncovprob;i++)
		{
			if(p->Next->vars[i])
			{
				covprob[i] += p->Next->laplace;
			}
		}
		p = p->Next;
	}

	return;
}

void CModel::AddLaplace(int len,int* amodel,string modelPath,double alaplace)
{
	int i;
	CModel* p = this;
	CModel* pnext;
	
	//do not add a model that is already known
	//not to be relevant
	if(0!=this->laplace)
	{
		if(alaplace>this->laplace)
		{
			this->laplace = alaplace;
		}
		else if(exp(alaplace-this->laplace)<mdCutoffMin)
		{
			return;
		}
	}
	else
	{
		this->laplace = alaplace;
	}
	
	//keep the total number of models in the counter
	//of the first cell
	p->Counter++;
	
	while(NULL!=p->Next)
	{
		int r = LexCompare(len,amodel,p->Next->vars);
		if(0==r)//found this model
		{
			return;
		}
		else if(-1==r)
		{
			p = p->Next;
		}
		else //1==r
		{
			break;
		}
	}
		
	CModel* newp = new CModel;
	newp->Length = len;
	newp->vars = new int[len];
	newp->modelPath = modelPath;
	newp->laplace = alaplace;
	for(i=0;i<len;i++)
	{
		newp->vars[i] = amodel[i];
	}
	newp->Next = p->Next;
	p->Next = newp;
		
	return;
}

void CModel::DeleteAll()
{
	CModel* p = Next;
	CModel* pnext;
	
	while(NULL!=p)
	{
		pnext = p->Next;
		delete p;
		p = pnext;
	}
	Next = NULL;
	return;
}

void CModel::SaveAll(FILE* file)
{
	int i;
	CModel* p = Next;
	
	while(NULL!=p)
	{
		fprintf(file,"%d",p->vars[0]);
		for(i=1;i<p->Length;i++)
		{
			fprintf(file," %d",p->vars[i]);
		}
		fprintf(file," %.10lf\n",p->laplace);//((double)p->Counter)/((double)Counter));
		p = p->Next;
	}
	return;
}

void CModel::SaveBestModel(FILE* file)
{
	int i;
	CModel* p = getBestModel();
	
	if(NULL!=p)
	{
		fprintf(file,"%d",p->vars[0]);
		for(int i=1;i<p->Length;i++)
		{
			fprintf(file,"\t%d",p->vars[i]);
		}
		fprintf(file,"\t%.10lf\n",p->laplace);
	}
	return;
}

void CModel::SaveRelevant(FILE* file)
{
	double maxprob;
	int nrelmodels;
	Diagnostics(maxprob,nrelmodels);
	
	fprintf(file,"%d\n",nrelmodels);
	CModel* p = Next;
	while(NULL!=p)
	{
		if(exp(p->laplace-maxprob)>=mdCutoffMax)
		{
			fprintf(file,"%d",p->vars[0]);
			for(int i=1;i<p->Length;i++)
			{
				fprintf(file,"\t%d",p->vars[i]);
			}
			fprintf(file,"\t%.10lf\n",p->laplace);
		}
		p = p->Next;
	}
	
	return;
}

void CModel::DeleteNotRelevant()
{
	double maxprob;
	int nrelmodels;
	
	if(0.0!=this->laplace)
	{
		maxprob = this->laplace;
	}
	else
	{
		Diagnostics(maxprob,nrelmodels);
	}

	CModel* p = this;
	CModel* pnext;
	
	while(NULL!=p->Next)
	{
		double s = exp(p->Next->laplace-maxprob);
		if(s<mdCutoffMax)
		{
			pnext = p->Next;
			p->Next = p->Next->Next;
			delete pnext;
		}
		else
		{
			p = p->Next;
		}
	}
	
	return;
}


CModel* CModel::getBestModel()
{
	CModel* p = this->Next;
	if(NULL==p)
	{
		return(NULL);
	}
	
	double maxw = p->laplace;
	CModel* bestp = p;
	while(NULL!=p->Next)
	{
		if(p->Next->laplace>maxw)
		{
			maxw = p->Next->laplace;
			bestp = p->Next;
		}
		p = p->Next;
	}
	return(bestp);
}

int LexCompare(int len,int* m1,int* m2)
{
	int i;
	
	for(i=0;i<len;i++)
	{
		if(m1[i]>m2[i])
		{
			return(-1);
		}
		if(m1[i]<m2[i])
		{
			return(1);
		}
	}
	return(0);
}

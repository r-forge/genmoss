/*
 Functions that are need to specify a hierarchical model.
*/

void InitLattice(const int n,
				 int** VarSets,int* lenVarSets,int nVarSets,
				 int**& DownLinks,int*& nDownLinks, 
				 int**& UpLinks,int*& nUpLinks)
{
	int i,j,k;
	int index[n];
	
	DownLinks = new int*[nVarSets];
	nDownLinks = new int[nVarSets];
	memset(nDownLinks,0,nVarSets*sizeof(int));
	
	UpLinks = new int*[nVarSets];
	nUpLinks = new int[nVarSets];
	memset(nUpLinks,0,nVarSets*sizeof(int));
	
	//need a dichotomous table
	for(i=0;i<n;i++) index[i]=2;
	LPTable p = new Table;
	p->Alloc(index,n);
	
	//down links
	for(i=0;i<nVarSets;i++)
	{
		k = 0;
		for(j=0;j<n;j++)
		{
			index[j] = VarSets[i][j];
			if(1==index[j])
			{
				k++;
			}
		}
		
		nDownLinks[i] = k;
		if(0==k) //nothing to do
		{
			DownLinks[i] = NULL;
		}
		else
		{
			DownLinks[i] = new int[k];
			memset(DownLinks[i],0,k*sizeof(int));
			
			k = 0;
			for(j=0;j<n;j++)
			{
				if(1==index[j])
				{
					index[j] = 0;
					DownLinks[i][k] = p->GetIndex(index);
					k++;
					index[j] = 1;
				}
			}
		}
	}
	
	//uplinks
	for(i=0;i<nVarSets;i++)
	{
		k = 0;
		for(j=0;j<n;j++)
		{
			index[j] = VarSets[i][j];
			if(0==index[j])
			{
				k++;
			}
		}
		
		nUpLinks[i] = k;
		if(0==k) //nothing to do
		{
			UpLinks[i] = NULL;
		}
		else
		{
			UpLinks[i] = new int[k];
			memset(UpLinks[i],0,k*sizeof(int));
			
			k = 0;
			for(j=0;j<n;j++)
			{
				if(0==index[j])
				{
					index[j] = 1;
					UpLinks[i][k] = p->GetIndex(index);
					k++;
					index[j] = 0;
				}
			}
		}		
	}
	
	p->Reset(); delete p; p = NULL;
	return;
}

void DeleteLattice(int** VarSets,int* lenVarSets,int nVarSets,
				   int**& DownLinks,int*& nDownLinks, 
				   int**& UpLinks,int*& nUpLinks)
{
	int i;
	
	for(i=0;i<nVarSets;i++)
	{
		if(nDownLinks[i]>=1)
		{
			delete[] DownLinks[i];
			DownLinks[i] = NULL;
		}
		if(nUpLinks[i]>=1)
		{
			delete[] UpLinks[i];
			UpLinks[i] = NULL;
		}
	}
	delete[] DownLinks; DownLinks = NULL;
	delete[] nDownLinks; nDownLinks = NULL;
	delete[] UpLinks; UpLinks = NULL;
	delete[] nUpLinks; nUpLinks = NULL;
	return;
}

//generates a random model
void RandomModel(int* amodel,const int lenmodel,int** VarSets,int* lenVarSets,double prob, gsl_rng * stream)
{
	int i;

	int r[lenmodel];
	
	//viRngBernoulli(METHOD, stream, lenmodel, r, prob);
	for(i=0;i<lenmodel;i++)
	{
		r[i] = gsl_ran_bernoulli(stream, prob);
		amodel[i] = 0;
	}
	for(i=1;i<lenmodel;i++)
	{
		amodel[i] = r[i];
	}
	
	//make sure all the main effects are in the model
	for(i=0;i<lenmodel;i++)
	{
		if(1==lenVarSets[i])
		{
			amodel[i] = 1;
		}
	}
	return;
}

void MakeModelHierarchical(int node,int mark,int* amodel,int** DownLinks,int* nDownLinks)
{
	int i;
	int newmark = 0;
	
	if(node>=1)
	{
		if(1==mark)
		{
			amodel[node] = 1;
		}
		newmark = amodel[node];
		for(i=0;i<nDownLinks[node];i++)
		{
			MakeModelHierarchical(DownLinks[node][i],newmark,
								  amodel,DownLinks,nDownLinks);
		}
	}
	return;
}

void PrettyWriteModel(int* amodel,int lenmodel,int** VarSets,int n)
{
	int i,j;
	
	for(i=0;i<lenmodel;i++)
	{
		if(1==amodel[i])
		{
			printf("%d",i+1);
			for(j=0;j<n;j++)
			{
				if(1==VarSets[i][j])
				{
					printf("\t%d",j+1);
				}
			}
			printf("\n");
		}
	}
	return;
}

int sprintfmodel(char* buffer,int bufsize,int* amodel,int lenmodel,int** VarSets,int n)
{
	int i,j,k;
	
	k = 0;
	for(i=0;i<lenmodel;i++)
	{
		if(1==amodel[i])
		{
			k += sprintf(buffer+k,"[");
			int first = 1;
			for(j=0;j<n;j++)
			{
				if(1==VarSets[i][j])
				{
					if(first)
					{
						first = 0;
						k += sprintf(buffer+k,"%d",j+1);
					}
					else
					{
						k += sprintf(buffer+k,",%d",j+1);
					}
				}
			}
			k += sprintf(buffer+k,"]");
		}
	}
	
	return(k<=bufsize);
}

void findGenerators(int node,int* amodel,int lenmodel,int** UpLinks,int* nUpLinks,int* generators)
{
	int i;
	
	int ismaximal = 1;
	for(i=0;i<nUpLinks[node];i++)
	{
		if(1==amodel[UpLinks[node][i]])
		{
			if(1==ismaximal)
			{
				ismaximal = 0;
			}
			findGenerators(UpLinks[node][i],amodel,lenmodel,UpLinks,nUpLinks,generators);
		}
	}
	generators[node] = ismaximal;
	return;
}


void findDualGenerators(int node,int* amodel,int lenmodel,int** DownLinks,int* nDownLinks,int* generators)
{
	int i;
	
	int isminimal = 1;
	for(i=0;i<nDownLinks[node];i++)
	{
		if(0==amodel[DownLinks[node][i]])
		{
			if(1==isminimal)
			{
				isminimal = 0;
			}
			findDualGenerators(DownLinks[node][i],amodel,lenmodel,DownLinks,nDownLinks,generators);
		}
	}
	generators[node] = isminimal;
	return;
}

void WriteModel(FILE* fid,int startpoint,int* amodel,int lenmodel)
{
	int i;
	
	fprintf(fid,"%d",amodel[0]);
	for(i=1;i<lenmodel;i++)
	{
		fprintf(fid," %d",amodel[i]);
	}
	fprintf(fid," %d\n",startpoint);
	return;	
}

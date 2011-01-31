/*
 *  varsets.functions.cpp
 *  
 *
 *  Created by Adrian Dobra on 10/8/06.
 *  Copyright 2006 __University of Washington__. All rights reserved.
 *
 */

void InitVarSets(int n,int**& VarSets,int*& lenVarSets,int& nVarSets)
{
    int i,j;
	int index[n];

	for(i=0;i<n;i++) index[i]=2;
	LPTable p = new Table;
	p->Alloc(index,n);

	nVarSets = (int)(pow(2.0,n));
	VarSets = new int*[nVarSets];
	for(i=0;i<nVarSets;i++)
	{
		VarSets[i] = new int[n];
	}
	lenVarSets = new int[nVarSets];
	
	p->GetFirst();
	j = 0;
	while(1)
	{
		lenVarSets[j] =0;
		for(i=0;i<n;i++)
		{
			VarSets[j][i]=p->Index[i];
			lenVarSets[j] += p->Index[i];
		}
		j++;
		if(!p->GetNext())
		{
			break;
		}
	}

	p->Reset(); delete p; p = NULL;
	return;
}


void DeleteVarSets(int**& VarSets,int*& lenVarSets,int& nVarSets)
{
	int i;
	
	for(i=0;i<nVarSets;i++)
	{
		delete[] VarSets[i];
		VarSets[i] = NULL;
	}
	delete[] VarSets;
	VarSets = NULL;
	delete[] lenVarSets;
	lenVarSets = NULL;
	return;
}

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

int intersect(int n,int* set1,int* set2,int* set12)
{
	int i;
	int r = 0;
	
	for(i=0;i<n;i++)
	{
		set12[i] = set1[i]*set2[i];
		r += set12[i];
	}
	return(r);
}



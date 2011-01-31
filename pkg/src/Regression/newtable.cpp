#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "newtable.h"

#define MAXPATHLENGTH 1024

/*
void CheckPointer(void* pointer)
{
   if(pointer==NULL)
   {
      printf("Memory allocation error.\n");
      exit(1);
   }
   return;
}

int choose(int m, int n)
{
   if(m > n) return 0;
   if(m == 0) return 1;
   if(m == 1) return n;
   if(m == n) return 1;
   if(m == (n-1)) return n;
   return(choose(m,n-1)+choose(m-1, n-1));
}
*/

///////////////////////////////////////////////////////////

//class Table :: BEGIN
Table::Table()
{
	nDimens = 0;
	Total = 0;
	Dimens = NULL;
	nAr = NULL;
	Data = NULL;
	Index = NULL;
	GrandTotal = 0.0;
	return;
}

Table::~Table()
{
   if((nDimens != 0)||(Total != 0)||(Dimens != NULL)||
      (Index != NULL) || (nAr != NULL) ||(Data != NULL) ||
      (GrandTotal != 0.0))
   {
      printf("Table::You didn't free the memory!!\n");
   }
   return;
}

//returns the total size of the table
//i.e. the product of all the dimensions
int Table::Alloc(int* dim, int ndim)
{
	int i,s;
	
	Reset(); //cleans the memory first
	if(ndim < 1)
	{
		printf("Invalid dimension!!\n");
		return 0;
	}
	nDimens = ndim;
	Dimens  = new int[ndim]; CheckPointer(Dimens);
	nAr     = new int[ndim]; CheckPointer(nAr);
	Index   = new int[ndim]; CheckPointer(Index);
	for(i=0; i<nDimens; i++)
	{
		Dimens[i] = dim[i];
	}
	s = 1;
    //compute nAr
	for(i = nDimens - 1; i>=0; i--)
	{
		nAr[i] = s;
		s *= Dimens[i];
	}
	Total = s;
        //allocate storage
	Data = new double[Total]; CheckPointer(Data);
        memset(Data,0,Total*sizeof(double));
	return Total;
}

void Table::Reset()
{
	nDimens = 0;
	Total = 0;
	delete[] Dimens; Dimens = NULL;
	delete[] nAr;    nAr = NULL;
	delete[] Data;   Data = NULL;
	delete[] Index;  Index = NULL;
	GrandTotal = 0.0;
	return;
}

int Table::ReadTable(const char* sFileName)
{
	FILE* in;
	int ndim;
	int dim[MAXDIMENS];
	int s;
	double rs;
	int i;
	
	if(NULL == (in = fopen(sFileName, "r")))
	{
		printf("Could not open file %s!!\n", sFileName);
		return 0;
	}
	if(!fscanf(in, "%d", &s))
	{
		printf("File %s has an unknown format!\n", sFileName);
		fclose(in);
		return 0;
	}
	if(s>MAXDIMENS)
	{
		printf("Increase MAXDIMENS!\n");
		return 0;
	}
	ndim = s;
     //read the dimensions
	for(i=0; i<ndim; i++)
	{
		if(!fscanf(in, "%d", &s))
		{
			printf("File %s has an unknown format!\n", sFileName);
			fclose(in);
			return 0;
		}
		dim[i] = s;
	}
	if(!Alloc(dim, ndim))
	{
		fclose(in);
		return 0;
	}
	for(i=0; i<Total; i++)
	{
		if(!fscanf(in, "%lf", &rs))
		{
			printf("File %s has an unknown format!\n", sFileName);
			fclose(in);
			return 0;
		}
		Data[i] = rs;
	}
	fclose(in);
	return 1;
}

int Table::ReadTableConsole(const char* sFileName)
{
	int ndim;
	int dim[MAXDIMENS];
	int s;
	double rs;
	int i, j;
	
	printf("Number of Dimensions :: ");
	if(!scanf("%d", &s))
	{
		return 0;
	}
	if(s > MAXDIMENS)
	{
		printf("Increase MAXDIMENS!\n");
		return 0;
	}
	ndim = s;
    //read the dimensions
	for(i=0; i<ndim; i++)
	{
		printf("Index %d :: ", i+1);
		if(!scanf("%d", &s))
		{
			return 0;
		}
		dim[i] = s;
	}
	if(!Alloc(dim, ndim))
	{
		return 0;
	}
	
	int NotFinished = 1;
	GetFirst();
	for(i=0; i<Total; i++)
	{
		printf("(");
		for(j=0; j<nDimens-1; j++)
		{
			printf("%d,", Index[j]+1);
		}
		printf("%d) :: ", Index[nDimens-1]+1);
		if(!scanf("%lf", &rs))
		{
			return 0;
		}
		Data[i] = rs;
		if(NotFinished)
		{
			NotFinished = GetNext();
		}
	}
	return WriteTable(sFileName);
}

int Table::WriteTable(const char* sFileName)
{
	FILE* out;
	int i;
	
	if(NULL == (out = fopen(sFileName, "w")))
	{
		printf("Could not open file %s!!\n", sFileName);
		return 0;
	}
	fprintf(out, "%d\n", nDimens);
     //save dimensions
	for(i=0; i<nDimens; i++)
	{
		fprintf(out, "%d  ", Dimens[i]);
	}
	fprintf(out, "\n");
     //save data
	for(i=0; i<Total; )
	{
		fprintf(out, "%.10lf  ", Data[i]);
		i++;
		if(i%5 == 0)
		{
			fprintf(out, "\n");
		}
	}
	fclose(out);
	return 1;
}

int Table::PrettyWriteTable(const char* sFileName)
{
	FILE* out;
	int i,j;
	
	if(NULL == (out = fopen(sFileName, "w")))
	{
		printf("Could not open file %s!!\n", sFileName);
		return 0;
	}
	
	GetFirst();
    for(i=0;i<Total;i++)
    {
      fprintf(out,"[%d,",Index[0]+1);
      for(j=1;j<nDimens-1;j++) fprintf(out,"%d,",Index[j]+1);
      fprintf(out,"%d] :: ",Index[nDimens-1]+1);
      fprintf(out,"%.10lf\n",Data[i]);
	  GetNext();
	}
	fclose(out);
	return 1;
}


double Table::GetI(int* index)
{
	int i;
	int s = 0;
	for(i=0; i<nDimens; i++)
	{
		s += nAr[i] * index[i];
	}
	return(*(Data + s));
}

int Table::GetIndex(int* index)
{
	int i;
	int s = 0;
	for(i=0; i<nDimens; i++)
	{
		s += nAr[i] * index[i];
	}
	return(s);
}

double Table::Get()
{
	int i;
	int s = 0;
	for(i=0; i<nDimens; i++)
	{
		s += nAr[i] * Index[i];
	}
	return(*(Data + s));
}

void Table::Set(double val)
{
	int i;
	int s = 0;
	for(i=0; i<nDimens; i++)
	{
		s += nAr[i] * Index[i];
	}
	*(Data + s) = val;
	return;
}

void Table::SetIndex(int* index)
{
	int i;
	for(i=0; i<nDimens; i++)
	{
		Index[i] = index[i];
	}
	return;
}

void Table::GetFirst()
{
	for(int i=0; i<nDimens; i++)
	{
		Index[i] = 0;
	}
	return;
}

void Table::GetFirstLoglin()
{
   for(int i=0;i<nDimens;i++)
   {
      Index[i] = 1;
   }
   return;
}

void Table::GetFirstWithIndex(int howmany,int* positions,int* values)
{
   int i;

   for(i=0;i<nDimens;i++)
   {
      Index[i]=0;
   }
   for(i=0;i<howmany;i++)
   {
      Index[positions[i]]=values[i];
   }

   return;
}

//return 0 if it's over
int Table::GetNext()
{
   int j = nDimens - 1;
   while(j >= 0)
   {
      if(Index[j] < Dimens[j] - 1)
      {
         Index[j]++;
	 return 1;
      }
      else
      {
         Index[j] = 0;
	 j--;
      }
   }
   return 0;
}

int Table::GetNextLoglin()
{
   int j = nDimens - 1;
   while(j >= 0)
   {
      if(Index[j] < Dimens[j] - 1)
      {
         Index[j]++;
	 return 1;
      }
      else
      {
         Index[j] = 1;
	 j--;
      }
   }
   return 0;
}

int  Table::GetNextWithIndex(int howmany,int* positions)
{
   int j=nDimens-1;
   int k=howmany-1;
   while(j>=0)
   {
	  if(k>=0)
	  {
	     if(j==positions[k])
		 {
			j--;
			k--;
			continue;
		 }
	  }
      if(Index[j]<Dimens[j]-1)
      {
		 Index[j]++;
         return 1;
      }
      else
      {
		 Index[j]=0;
         j--;
      }
   }
   return 0;
}

double Table::GetGrandTotal()
{
  //if(GrandTotal > 0)
  //{
  //	return GrandTotal;
  //}
	int i;
	double s = 0.0;
	for(i=0; i<Total; i++)
	{
		s += Data[i];
	}
	GrandTotal = s;
	return GrandTotal;
}

int Table::CheckCell(double minp,int acell)
{
	int i;
	double s = log(GetGrandTotal())+log(minp);
	
	if(log(Data[acell])<s)
	{
		return(0);
	}
	return(1);
}

int Table::CheckTable(double minp)
{
	int i;
	double s = log(GetGrandTotal())+log(minp);

	for(i=0;i<Total;i++)
	{
		if(log(Data[i])<s)
		{
			return(0);
		}
	}
	return(1);
}

/*
   Creates a new table from the initial table.
   The new table has a dimension less than
   the original table i.e. we have a '+'
   in the place of the dimension ind
*/

LPTable Table::ReduceOne(int ind)
{
	int i, j;
	int dim[MAXDIMENS]; //dimensions of the new table
	int ndim; //number of dimensions for the new table
	int NotFinished = 1;
	int* myIndex = new int[nDimens];CheckPointer(myIndex);
	memset(myIndex,0,nDimens*sizeof(int));
	double s;
	
	LPTable ctab = new Table;
	if(NULL == ctab)
	{
		printf("Error creating new table :: ReduceOne.\n");
		return(NULL);
	}
	ndim = nDimens - 1;
	j = 0;
	for(i=0; i<nDimens; i++)
	{
		if(i != ind)
		{
			dim[j] = Dimens[i];
			j++;
		}
	}
     //allocate storage for the new table
	ctab->Alloc(dim, ndim);
	ctab->GetFirst();
	while(NotFinished)
	{
		j = 0;
		for(i=0; i<nDimens; i++)
		{
			if(i != ind)
			{
				myIndex[i] = ctab->Index[j];
				j++;
			}
		}
		
		myIndex[ind] = 0;
		s = 0.0;
		
		for(i=0; i<Dimens[ind]; i++)
		{
			s += GetI(myIndex);
			myIndex[ind]++;
		}
		ctab->Set(s);
		NotFinished = ctab->GetNext();
	}
	delete[] myIndex;
	return(ctab);
}

//pind is a vector of length m of dimensions
//the dimensions in pind are in an increasing order
LPTable Table::Reduce(int m, int* pind)
{
	int i;
	LPTable newtab;
	LPTable oldtab;
	
	if(m == 0)
	{
		return NULL;
	}
	newtab = ReduceOne(pind[m-1]);
	if(m == 1)
	{
		return newtab;
	}
	for(i=m-2; i>=0; i--)
	{
		oldtab = newtab;
		newtab = oldtab->ReduceOne(pind[i]);
		oldtab->Reset();
		delete oldtab;
	}
	return newtab;
}
//class Table :: END

//class NTable :: BEGIN
NTable::NTable()
{
	nParentDimens = 0;
	ConvertIndex = NULL;
	return;
}

NTable::~NTable()
{
	if((nParentDimens != 0) || (ConvertIndex != NULL))
	{
		printf("NTable::You didn't free the memory!!\n");
	}
	return;
}

int NTable::Alloc(int* dim, int ndim)
{
	if(NULL == (ConvertIndex = new int[ndim]))
	{
		printf("Could not allocate memory :: NTable :: Alloc.\n");
		return 0;
	}
	return Table::Alloc(dim, ndim);
}

void NTable::Reset()
{
	nParentDimens = 0;
	delete[] ConvertIndex; ConvertIndex = NULL;
	Table::Reset();
}

void NTable::SetConvertIndex(int* pind)
{
	int i;
	for(i=0; i<nDimens; i++)
	{
		ConvertIndex[i] = pind[i];
	}
	return;
}

double NTable::GetI(int* index)
{
	int i;
	int s = 0;
	for(i=0; i<nDimens; i++)
	{
		s += nAr[i] * index[ConvertIndex[i]];
	}
	return(*(Data + s));
}

//index has length nParentDimens
void NTable::SetIndex(int* index)
{
	int i;
	for(i=0; i<nDimens; i++)
	{
		Index[i] = index[ConvertIndex[i]];
	}
	return;
}

//reduce table by the dimensions NOT contained in c
//the vector c has length m
int NTable::Create(int m, int* c, LPTable table)
{
	int i, j, k;
	int index[MAXDIMENS];
	
	k = 0;
	for(i=0; i<c[0]; i++)
	{
		index[k] = i;
		k++;
	}
	j = 1;
	while(j < m)
	{
		for(i=c[j-1]+1; i<c[j]; i++)
		{
			index[k] = i;
			k++;
		}
		j++;
	}
	for(i=c[m-1]+1; i<table->nDimens; i++)
	{
		index[k] = i;
		k++;
	}
	
	if(k>=1)
	{
		LPTable newtab = table->Reduce(k, index);
		if(newtab == NULL)
		{
			printf("NTable :: Error in Create.\n");
			return 0;
		}
		if(!Alloc(newtab->Dimens, newtab->nDimens))
		{
			printf("NTable :: Could not allocate memory.\n");
			return 0;
		}
		for(i=0; i<Total; i++)
		{
			Data[i] = newtab->Data[i];
		}
		newtab->Reset(); delete newtab; newtab = NULL;
	}
	else
	{
		if(!Alloc(table->Dimens,table->nDimens))
		{
			printf("NTable :: Could not allocate memory.\n");
			return 0;
		}
		for(i=0;i<Total;i++)
		{
			Data[i] = table->Data[i];
		}
	}
	nParentDimens = table->nDimens;
	SetConvertIndex(c);
	return 1;
}
//class NTable :: END

//class STable :: BEGIN
STable::STable()
{
	nTables = 0;
	pTables = NULL;
}

STable::~STable()
{
	if((nTables != 0) || (pTables != NULL))
	{
		printf("STable::You didn't free the memory!!\n");
	}
}

int STable::Create(int m, LPTable table)
{
	//if m == table->nDimens we don't need a STable
	if((m < 1) || (table == NULL) || (m >= table->nDimens))
	{
		printf("STable :: Error in Create.\n");
		return 0;
	}
	
	int n = table->nDimens;
	int c[MAXDIMENS];
	int i, j;
	int ind;
	
	nTables = choose(m, n);
	
	pTables = new NTable[nTables];
	if(pTables == NULL)
	{
		printf("STable :: Error in Create.\n");
		return 0;
	}
	
   //generate all m-subsets of {0,...,n-1}
	for(i=0; i<m; i++)
	{
		c[i] = i;
	}
	
	ind = 0;
	while(1)
	{
      //create a new NTable based on table
      //and indexes given by the vector c
		if(!(pTables[ind].Create(m, c, table)))
		{
			printf("STable :: Failed to create NTable.\n");
			return 0;
		}
		
		for(j = m - 1; j>=0; )
		{
			if(c[j] == (n - m + j))
			{
				j--;
			}
			else
			{
				break;
			}
		}
		
		if(j < 0) break; //we've already generated all the subsets
		
		c[j]++;
		for(i=j+1; i<m; i++)
		{
			c[i] = c[i-1] + 1;
		}
		
		ind++;
	}
	
	return 1;
}

void STable::Reset()
{
	int i;
	
	for(i=0; i<nTables; i++)
	{
		pTables[i].Reset();
	}
	delete[] pTables;
	
	pTables = NULL;
	nTables = 0;
	return;
}

double STable::Get(int* index)
{
	double s = 0.0;
	int i;
	
	for(i=0; i<nTables; i++)
	{
		s += pTables[i].GetI(index);
	}
	return s;
}
double STable::GetMin(int* index)
{
	double s;
	double s1;
	int i;
	
	s = pTables[0].GetI(index);
	for(i=1; i<nTables; i++)
	{
		s1 = pTables[i].GetI(index);
		if(s1 < s)
		{
			s = s1;
		}
	}
	
	return s;
}
//class STable :: END

//BEGIN
LPTable CreateS(int m, LPTable parent)
{
	STable stab;
	
	if(!stab.Create(m, parent))
	{
		printf("Error creating STable :: CreateS.\n");
		return(NULL);
	}
	
	LPTable newtab = new Table;
	if(newtab == NULL)
	{
		printf("Error creating STable :: CreateS.\n");
		return(NULL);
	}
	
	if(!newtab->Alloc(parent->Dimens, parent->nDimens))
	{
		printf("Error creating STable :: CreateS.\n");
		return(NULL);
	}
	
	int NotFinished = 1;
	newtab->GetFirst();
	while(NotFinished)
	{
		newtab->Set(stab.Get(newtab->Index));
		NotFinished = newtab->GetNext();
	}
	
	stab.Reset();
	
	return newtab;
}

LPTable CreateMin(int m, LPTable parent)
{
	STable stab;
	
	if(!stab.Create(m, parent))
	{
		printf("Error creating STable :: CreateS.\n");
		return(NULL);
	}
	
	LPTable newtab = new Table;
	if(newtab == NULL)
	{
		printf("Error creating STable :: CreateS.\n");
		return(NULL);
	}
	
	if(!newtab->Alloc(parent->Dimens, parent->nDimens))
	{
		printf("Error creating STable :: CreateS.\n");
		return(NULL);
	}
	
	int NotFinished = 1;
	newtab->GetFirst();
	while(NotFinished)
	{
		newtab->Set(stab.GetMin(newtab->Index));
		NotFinished = newtab->GetNext();
	}
	
	stab.Reset();
	
	return newtab;
}

LPTable CreateNBar(LPTable parent)
{
	if(parent == NULL)
	{
		printf("Error :: CreateNBar.\n");
		return NULL;
	}
	
	LPTable newtab = new Table;
	if(newtab == NULL)
	{
		printf("Error :: CreateNBar\n");
		return(NULL);
	}
	
	if(!newtab->Alloc(parent->Dimens, parent->nDimens))
	{
		printf("Error creating STable :: CreateS.\n");
		return(NULL);
	}
	
	double n;
	n = parent->GetGrandTotal();
	int k;
	k = parent->nDimens;
	int m;
	int plus = -1; //we begin with '-'
	int first = 1; //first iteration?
	int i;
	
	for(m=1; m<k; m++)
	{
		LPTable smtab = CreateS(m, parent);
		if(smtab == NULL)
		{
			printf("Error :: CreateNBar.\n");
			return NULL;
		}
		
		if(first)
		{
			for(i=0; i<newtab->Total; i++)
			{
				newtab->Data[i] = 1 - smtab->Data[i]; //n
			}
			
			first = 0;
		}
		else //first == 0
		{
			if(plus == -1)
			{
				for(i=0; i<newtab->Total; i++)
				{
					newtab->Data[i] -= smtab->Data[i];
				}
			}
			else //plus == 1
			{
				for(i=0; i<newtab->Total; i++)
				{
					newtab->Data[i] += smtab->Data[i];
				}
			}
		}
		
		smtab->Reset();
		delete smtab;
		
		plus = - plus;
	}
	
	if(plus == -1)
	{
		for(i=0; i<newtab->Total; i++)
		{
			newtab->Data[i] -= parent->Data[i];
		}
	}
	else
	{
		for(i=0; i<newtab->Total; i++)
		{
			newtab->Data[i] += parent->Data[i];
		}
	}
	
	return newtab;
}

int FrechetBounds1(LPTable table, const char* sFileName)
{
	if(table == NULL)
	{
		printf("Error :: FrechetBounds1.\n");
		return 0;
	}
	
	int k = table->nDimens;
	LPTable UpperBound = CreateMin(1, table);
	if(UpperBound == NULL)
	{
		printf("Error :: FrechetBounds1.\n");
		return 0;
	}
	
	LPTable LowerBound = CreateS(1, table);
	if(LowerBound == NULL)
	{
		printf("Error :: FrechetBounds1.\n");
		UpperBound->Reset();
		delete UpperBound;
		return 0;
	}
	
	double n = table->GetGrandTotal();
	int i;
	
	for(i=0; i<table->Total; i++)
	{
		LowerBound->Data[i] -= n*(k-1);
		if(LowerBound->Data[i] < 0.0)
		{
			LowerBound->Data[i] = 0.0;
		}
	}
	
	int rez = WriteBounds(UpperBound, LowerBound, sFileName);
	
	UpperBound->Reset();
	LowerBound->Reset();
	delete UpperBound;
	delete LowerBound;
	
	return rez;
}

int FrechetBounds(LPTable table, const char* sFileName)
{
	if(table == NULL)
	{
		printf("Error :: FrechetBounds.\n");
		return 0;
	}
	
	int k = table->nDimens;
	LPTable UpperBound = CreateMin(k-1, table);
	if(UpperBound == NULL)
	{
		printf("Error :: FrechetBounds.\n");
		return 0;
	}
	
	LPTable LowerBound = new Table;
	if(LowerBound == NULL)
	{
		printf("Error :: FrechetBounds.\n");
		UpperBound->Reset();
		delete UpperBound;
		return 0;
	}
	
	if(!LowerBound->Alloc(table->Dimens, table->nDimens))
	{
		printf("Error :: FrechetBounds.\n");
		UpperBound->Reset();
		delete UpperBound;
		return 0;
	}
	
	double n = table->GetGrandTotal();
	int m;
	int plus;
	int first = 1; //first iteration?
	int i;
	
	if(k%2 == 0) //k even
	{
		plus = 1;
		
		for(m=1; m<k; m++)
		{
			LPTable smtab = CreateS(m, table);
			if(smtab == NULL)
			{
				printf("Error :: FrechetBounds.\n");
				UpperBound->Reset();
				LowerBound->Reset();
				delete UpperBound;
				delete LowerBound;
				
				return 0;
			}
			
			if(first)
			{
				for(i=0; i<table->Total; i++)
				{
					LowerBound->Data[i] = smtab->Data[i] - n;
				}
				
				first = 0;
			}
			else
			{
				if(plus == -1)
				{
					for(i=0; i<table->Total; i++)
					{
						LowerBound->Data[i] -= smtab->Data[i];
					}
				}
				else //plus == 1
				{
					for(i=0; i<table->Total; i++)
					{
						LowerBound->Data[i] += smtab->Data[i];
					}
				}
			}
			
			plus = - plus;
			smtab->Reset();
			delete smtab;
		}
	}
	else // k odd
	{
		LPTable nBar = CreateNBar(table);
		if(nBar == NULL)
		{
			printf("Error :: FrechetBounds.\n");
			UpperBound->Reset();
			LowerBound->Reset();
			delete UpperBound;
			delete LowerBound;
			return 0;
		}
		
		LPTable nMinBar = CreateMin(k-1, nBar);
		nBar->Reset();
		delete nBar;
		if(nMinBar == NULL)
		{
			printf("Error :: FrechetBounds.\n");
			UpperBound->Reset();
			LowerBound->Reset();
			delete UpperBound;
			delete LowerBound;
			return 0;
		}
		
		plus = -1;
		
		for(m=1; m<k; m++)
		{
			LPTable smtab = CreateS(m, table);
			if(smtab == NULL)
			{
				printf("Error :: FrechetBounds.\n");
				UpperBound->Reset();
				LowerBound->Reset();
				delete UpperBound;
				delete LowerBound;
				return 0;
			}
			
			if(first)
			{
				for(i=0; i<table->Total; i++)
				{
					LowerBound->Data[i] = n - smtab->Data[i];
				}
				
				first = 0;
			}
			else
			{
				if(plus == -1)
				{
					for(i=0; i<table->Total; i++)
					{
						LowerBound->Data[i] -= smtab->Data[i];
					}
				}
				else //plus == 1
				{
					for(i=0; i<table->Total; i++)
					{
						LowerBound->Data[i] += smtab->Data[i];
					}
				}
			}
			
			plus = - plus;
			smtab->Reset();
			delete smtab;
		}
		
		for(i=0; i<table->Total; i++)
		{
			LowerBound->Data[i] -= nMinBar->Data[i];
		}
		
		nMinBar->Reset();
		delete nMinBar;
	}
	
	for(i=0; i<table->Total; i++)
	{
		if(LowerBound->Data[i] < 0)
		{
			LowerBound->Data[i] = 0;
		}
	}
	
	int rez = WriteBounds(UpperBound, LowerBound, sFileName);
	
	UpperBound->Reset();
	LowerBound->Reset();
	delete UpperBound;
	delete LowerBound;
	
	return rez;
}

int Bonferroni(int m, LPTable table, const char* sFileName)
{
	if(table == NULL)
	{
		printf("Error :: BonferroniBounds.\n");
		return 0;
	}
	
	int i;
	double nGrandTotal;
	nGrandTotal = table->GetGrandTotal();
	
	for(i=0; i<table->Total; i++)
	{
		table->Data[i] = nGrandTotal;
	}  
	
	table->WriteTable("tabinit.dat");
	
	LPTable UpperBound;
	LPTable LowerBound = new Table;
	if(LowerBound == NULL)
	{
		printf("Error :: BonferroniBounds.\n");
		return 0;
	}
	
	if(!LowerBound->Alloc(table->Dimens, table->nDimens))
	{
		printf("Error :: BonferroniBounds.\n");
		return 0;
	}
	
	LPTable nBar = CreateNBar(table);
	if(nBar == NULL)
	{
		printf("Error :: Bonferroni.\n");
		return 0;
	}
	
	nBar->WriteTable("nbar.dat");
	
	double nBarGrandTotal;
	nBarGrandTotal = nBar->GetGrandTotal();
	int m1;
	int plus;
	
	UpperBound = CreateS(1, nBar);
	if(UpperBound == NULL)
	{
		printf("Error :: Bonferroni.\n");
		nBar->Reset(); delete nBar;
		LowerBound->Reset(); delete LowerBound;
		return 0;
	}
	
	UpperBound->WriteTable("s1bar.dat");
	
	for(i=0; i<table->Total; i++)
	{
		UpperBound->Data[i] = 1 - UpperBound->Data[i]; //nBarGrandTotal
	}
	
	plus = 1;
	
	for(m1=2; m1<=m-1; m1++)
	{
		LPTable sm1 = CreateS(m1, nBar);
		if(sm1 == NULL)
		{
			printf("Error :: Bonferroni.\n");
			nBar->Reset(); delete nBar;
			UpperBound->Reset(); delete UpperBound;
			LowerBound->Reset(); delete LowerBound;
			return 0;
		}
		
		if(plus == 1)
		{
			for(i=0; i<table->Total; i++)
			{
				UpperBound->Data[i] += sm1->Data[i];
			}
		}
		else //plus == -1
		{
			for(i=0; i<table->Total; i++)
			{
				UpperBound->Data[i] -= sm1->Data[i];
			}
		}
		
		sm1->Reset(); delete sm1;
		plus = - plus;
	}
	
	for(i=0; i<table->Total; i++)
	{
		LowerBound->Data[i] = UpperBound->Data[i];
	}
	
	if(m < table->nDimens)
	{
		LPTable sm = CreateS(m, nBar);
		if(sm == NULL)
		{
			printf("Error :: Bonferroni.\n");
			nBar->Reset(); delete nBar;
			UpperBound->Reset(); delete UpperBound;
			LowerBound->Reset(); delete LowerBound;
			return 0;
		}
		
		if(plus == 1) //m even
		{
			for(i=0; i<table->Total; i++)
			{
				UpperBound->Data[i] += sm->Data[i];
			}
		}
		else //m odd
		{
			for(i=0; i<table->Total; i++)
			{
				LowerBound->Data[i] -= sm->Data[i];
			}
		}
		sm->Reset();
		delete sm;
	}
	else
	{
		if(plus == 1) //m even
		{
			for(i=0; i<table->Total; i++)
			{
				UpperBound->Data[i] += nBar->Data[i];
			}
		}
		else //m odd
		{
			for(i=0; i<table->Total; i++)
			{
				LowerBound->Data[i] -= nBar->Data[i];
			}
		}
	}
	
	WriteBounds(UpperBound, LowerBound, "extra.dat");
	
	double ratio = nGrandTotal; // /nBarGrandTotal;
	printf("ratio = %lf\n", ratio);
	for(i=0; i<table->Total; i++)
	{
		UpperBound->Data[i] = ceil(UpperBound->Data[i]*ratio);
		LowerBound->Data[i] = floor(LowerBound->Data[i]*ratio);
		if(LowerBound->Data[i] < 0)
		{
			LowerBound->Data[i] = 0;
		}
	}
	
	int rez = WriteBounds(UpperBound, LowerBound, sFileName);
	
	nBar->Reset();
	UpperBound->Reset();
	LowerBound->Reset();
	delete nBar;
	delete UpperBound;
	delete LowerBound;
	
	return rez;
}

LPTable ReduceOneShuttle(LPTable tab, int ind, LPTable tabS)
{
	LPTable tabR = tab->ReduceOne(ind);
	if(NULL == tabR)
	{
		printf("Error creating new table :: ReduceOneShuttle.\n");
		return(NULL);
	}
	
	LPTable rez = new Table;
	if(NULL == rez)
	{
		printf("Error creating new table :: ReduceOneShuttle.\n");
		tabR->Reset(); delete tabR;
		return(NULL);
	}
	
	if(!rez->Alloc(tab->Dimens, tab->nDimens))
	{
		printf("Error creating new table :: ReduceOneShuttle.\n");
		tabR->Reset(); delete tabR;
		delete rez;
		return(NULL);
	}
	
	rez->GetFirst();
	int i, j;
	int NotFinished = 1;
	while(NotFinished)
	{
		j = 0;
		for(i=0; i<rez->nDimens; i++)
		{
			if(i != ind)
			{
				tabR->Index[j] = rez->Index[i];
				j++;
				tabS->Index[i] = rez->Index[i];
			}
		}
		
		double s = 0.0;
		for(i=0; i<rez->Dimens[ind]; i++)
		{
			if(i != rez->Index[ind])
			{
				tabS->Index[ind] = i;
				s += tabS->Get();
			}
		}
		
		rez->Set(tabR->Get() - s);
		
		if(NotFinished)
		{
			NotFinished = rez->GetNext();
		}
	}
	
	tabR->Reset(); delete tabR;
	
	return rez;
}

int ShuttleBounds(LPTable table, int nIterations, const char* sFileName)
{
	int k = table->nDimens;
	LPTable UpperBound = CreateMin(k-1, table);
	if(UpperBound == NULL)
	{
		printf("Error :: ShuttleBounds.\n");
		return 0;
	}
	LPTable LowerBound = new Table;
	if(LowerBound == NULL)
	{
		printf("Error :: ShuttleBounds.\n");
		UpperBound->Reset(); delete UpperBound;
		return 0;
	}
	if(!LowerBound->Alloc(table->Dimens, table->nDimens))
	{
		printf("Error :: ShuttleBounds.\n");
		UpperBound->Reset(); delete UpperBound;
		return 0;
	}
	int i;
	for(i=0; i<table->Total; i++)
	{
		LowerBound->Data[i] = 0;
	}
	
	int ind;
	int iteration;
	for(iteration=0; iteration<nIterations; iteration++)
	{
		LPTable OldUpperBound = UpperBound;
		LPTable OldLowerBound = LowerBound;
		
		UpperBound = ReduceOneShuttle(table, 0, OldLowerBound);
		if(UpperBound == NULL)
		{
			OldUpperBound->Reset(); delete OldUpperBound;
			OldLowerBound->Reset(); delete OldLowerBound;
			return 0;
		}
		for(ind=1; ind<k; ind++)
		{
			LPTable newtab = ReduceOneShuttle(table, ind, OldLowerBound);
			if(newtab == NULL)
			{
				UpperBound->Reset();    delete UpperBound;
				OldUpperBound->Reset(); delete OldUpperBound;
				OldLowerBound->Reset(); delete OldLowerBound;
				return 0;
			}
			for(i=0; i<table->Total; i++)
			{
				if(newtab->Data[i] < UpperBound->Data[i])
				{
					UpperBound->Data[i] = newtab->Data[i];
				}
			}
			newtab->Reset(); delete newtab;
		}
		
		LowerBound = ReduceOneShuttle(table, 0, OldUpperBound);
		if(LowerBound == NULL)
		{
			UpperBound->Reset(); delete UpperBound;
			OldUpperBound->Reset(); delete OldUpperBound;
			OldLowerBound->Reset(); delete OldLowerBound;
			return 0;
		}
		
		for(ind=1; ind<k; ind++)
		{
			LPTable newtab = ReduceOneShuttle(table, ind, OldUpperBound);
			if(newtab == NULL)
			{
				UpperBound->Reset(); delete UpperBound;
				LowerBound->Reset(); delete LowerBound;
				OldUpperBound->Reset(); delete OldUpperBound;
				OldLowerBound->Reset(); delete OldLowerBound;
				return 0;
			}
			for(i=0; i<table->Total; i++)
			{
				if(newtab->Data[i] > LowerBound->Data[i])
				{
					LowerBound->Data[i] = newtab->Data[i];
				}
			}
			newtab->Reset(); delete newtab;
		}
		OldUpperBound->Reset(); delete OldUpperBound;
		OldLowerBound->Reset(); delete OldLowerBound;
	}
	int rez = WriteBounds(UpperBound, LowerBound, sFileName);
	UpperBound->Reset();
	LowerBound->Reset();
	delete UpperBound;
	delete LowerBound;
	return rez;
}

int WriteBounds(LPTable UpperBound, LPTable LowerBound, const char* sFileName)
{
	if((UpperBound == NULL) || (LowerBound == NULL))
	{
		printf("Error :: WriteBounds.\n");
		return 0;
	}
	
	FILE* out;
	if(NULL == (out = fopen(sFileName, "w")))
	{
		printf("Could not open file %s!!\n", sFileName);
		return 0;
	}
	
	int i, j;
	int k = UpperBound->nDimens;
	int NotFinished = 1;
	
	UpperBound->GetFirst();
	for(i=0; i<UpperBound->Total; i++)
	{
		fprintf(out, "(");
		for(j=0; j<k-1; j++)
		{
			fprintf(out, "%d,", UpperBound->Index[j]+1);
		}
		fprintf(out, "%d) :: (%.2lf, %.2lf)\n",
			UpperBound->Index[k-1]+1,
			UpperBound->Data[i], LowerBound->Data[i]);
		
		if(NotFinished)
		{
			NotFinished = UpperBound->GetNext();
		}
	}
	
	fclose(out);
	
	return 1;
}
//END

/*
int main()
{
   int i,j;

   LPTable table = new Table;
   CheckPointer(table);
   table->ReadTable("tab10.dat");

   int howmany = 1;
   int positions[6];
   int values[6];

   positions[0] = 0;
   values[0] = 0;

   table->GetFirstWithIndex(howmany,positions,values);
   for(j=0;j<table->nDimens;j++)
   {
     printf("%d ",table->Index[j]);
   }  
   printf("\n");

   int notdone = 1;
   int count = 1;
   while(1)
   {
     notdone = table->GetNextWithIndex(howmany,positions);
     if(!notdone) break;
     for(j=0;j<table->nDimens;j++)
     {
        printf("%d ",table->Index[j]);
     }
     printf("\n");
     printf("count = %d\n",++count);
   } 

   return 1;
}
*/




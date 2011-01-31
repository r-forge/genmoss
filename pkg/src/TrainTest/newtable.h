#ifndef _TABLE_H 
#define _TABLE_H

#define MAXDIMENS 100 //maximum number of dimensions

#include <cstdio>
#include <cstring>
#include <cstdlib>

int choose(int m, int n);
void CheckPointer(void* pointer);

class Table;
class NTable;
class STable;

typedef Table*  LPTable;
typedef NTable* LPNTable;
typedef STable* LPSTable;

//base class representing a table
class Table
{
public:
   //data members
	int nDimens; //number of dimensions
	int Total; //the total number of cells
	int* Dimens; //storage for dimensions
	int* nAr;
	double* Data; //storage for the useful data
	double GrandTotal;
	
	int* Index;//used to access the cells	

     //methods
	Table(); //contructor
	~Table(); //destructor
     //the destructor does not free the memory
	
	int  ReadTable(const char* sFileName); //reads a table from a file
	int  ReadTableConsole(char* sFileName); //reads a table from a file
	int  WriteTable(const char* sFileName); //writes the table in a file
	int  PrettyWriteTable(const char* sFileName); //writes the table with full info about cell indices
	int  Alloc(int* dim, int ndim); //allocates memory
	void Reset(); //frees memory
	
	double GetI(int* index); //returns the cell specified by index (NOT Index)
	int GetIndex(int* index); //returns the number of the cell specified by index
	double Get(); //gets the cell specified by Index
	void   Set(double val); //sets the value of the cell
                         //specified by Index
	void   SetIndex(int* index);
	void   GetFirst(); //use these two functions to iterate through
	int    GetNext(); //all the cells of the table
	double GetGrandTotal();
	int CheckTable(double minp);
	int CheckCell(double minp,int acell);

        //for LOG-LINEAR MODELS only
        void   GetFirstLoglin();
	int    GetNextLoglin();

        //methods for accessing only some elements	
        void GetFirstWithIndex(int howmany,int* positions,int* values);
        int  GetNextWithIndex(int howmany,int* positions); 

        //collapse a table across some variables
	LPTable ReduceOne(int ind);
	LPTable Reduce(int m, int* pind);
};

class NTable : public Table
{
public:
     //data members
	int nParentDimens;
	int* ConvertIndex;
	
     //methods
	NTable(); //constructor
	~NTable(); //destructor
	
	int Alloc(int* dim, int ndim);
	void Reset();
	
	void SetConvertIndex(int* pind); //initializes ConvertIndex
     //to be called after Alloc
     //pind has length nDimens and contains the indexes
     //of the original dimensions
	void SetIndex(int* index);
     //index is a vector of length nParentDimens
     //initializes Index by taking only the values from index
     //that correspond to the dimensions of the reduced table NTable
	double GetI(int* index);
	
	int Create(int m, int* c, LPTable table);
     //creates a NTable with parent table
     //and vector of dimensions c
};

class STable
{
public:
     //data members
	int nTables; //the number of NTables STable contains
	LPNTable pTables; //array of pointers to NTables
	
     //methods
	STable(); //contructor
	~STable(); //destructor
	
	int Create(int m, Table* table);
     //creates STable from the Table table
     //m-dimensional marginal bounds
	void Reset(); //frees the memory
	
	double Get(int* index);
     //returns the element specified by in index	
	double GetMin(int* index);
};

LPTable CreateS(int m, LPTable parent); //generates Sm from table parent
//generates min over all the components of Sm
LPTable CreateNBar(Table* parent);
int FrechetBounds1(LPTable table, char* sFileName);
//1-dimensional Frechet Marginal Bounds for k-way tables
int FrechetBounds(LPTable table, char* sFileName);
//(k-1)-dimensional Frechet Marginal Bounds for k-way tables
int Bonferroni(int m, LPTable table, char* sFileName);
//m-dimensional Bonferroni Marginal Bounds
LPTable ReduceOneShuttle(LPTable tab, int ind, LPTable tabS);
int ShuttleBounds(LPTable table, int nIterations, char* sFileName);
int WriteBounds(LPTable UpperBound, LPTable LowerBound, const char* sFileName);

/*
	The cells are stored in the file in the lexicographic order.
   Before the cells, the file contains the number of dimensions
   and the dimensions themselves.

*/

#endif

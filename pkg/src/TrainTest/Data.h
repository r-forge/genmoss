/*
AUTHOR: ADRIAN DOBRA (ADOBRA@STAT.DUKE.EDU)
		QUANLI WANG (QUANLI@STAT.DUKE.EDU)
*/

// Data.h: interface for the CData class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(DATA_H)
#define DATA_H

#include <cstdio>
#include <cstring>
#include <string>
#include <cstdlib>

using namespace std;
class CData  
{
public:
	CData();
	virtual ~CData();

	double* Y;
	double** X;
	double** data;

	int NumberOfGenes; //columns in the data
	int SampleSize; //rows in the data
	string DataFile;
	int ncolsX;

	void FreeData();
	void ReadData();
	void WriteData();
	void Cleanup();
	void Allocate(int ncolsx);
};

#endif

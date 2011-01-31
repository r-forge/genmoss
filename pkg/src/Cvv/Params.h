/*
AUTHOR: ADRIAN DOBRA (ADOBRA@STAT.WASHINGTON.EDU)
*/
// Params.h: interface for the Params class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(PARAMS_H)
#define PARAMS_H

#include <map>
#include <cstring>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

class Params  
{
public:
	Params();
	string& GetErrorMessage(){return mstrErrorMessage;};
	bool Load(string FileName);

public:
	int mnNumberOfGenes;

	string mstrDataFile;
	int mnSampleSize;
	
	int mnMaxRegressors;
	string mstrModelsFile;
	int mnChainIterations;
	int mnCvvFold;
	string mstrErrorMessage;

	static double ToDouble(string str);
	static int ToInt(string str);
	static string ToLower(string str);
	static string ToString(int value);
	static string ToString(double value);
	static string trim(string s);
private:
	int mnStartingGraph;
};

#endif

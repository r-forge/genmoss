/*
AUTHOR: ADRIAN DOBRA (ADOBRA@STAT.WASHINGTON.EDU)
*/
// Model.h: interface for the Model class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(MODEL_H)
#define MODEL_H

#include <map>
#include <cstring>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

class Model  
{
public:
	Model();
	string& GetErrorMessage(){return mstrErrorMessage;};
	bool Load(string FileName);

public:
	string mstrDataFile;
	int mnNumberOfGenes;
	int mnSampleSize;
	int mnMaxRegressors;
	//int mnMaxRegressorsInteraction;
	//int mnDoInteractions;
	//int mnDoSpaceRatio;
	//int mnDoMetropolisH;
	int mnChainIterations;
	//int mnChainIterationsBurnin;
	
	//model search for regressions
	int mnChainReplicates;
	double mdCutoffMax;
	double mdCutoffMin;
	double mdProbMax;

	int mnNumOfConfoundingVars;
	
	//model search for tables
	int mnShotgunChainReplicates;
	double mdShotgunCutoffMax;
	double mdShotgunCutoffMin;
	double mdShotgunProbMax;
	
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

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
	int mnNumberOfVariables;

	string mstrTrainDataFile;
	int mnTrainSampleSize;
	string mstrTrainOutputFile;
	string mstrScoresFile;
	//string mstrTrainIndexFile;

	string mstrTestDataFile;
	int mnTestSampleSize;
	string mstrTestOutputFile;
	//string mstrTestIndexFile;	

	double mdSplitsMin;
	double mdSplitsInc;
	double mdSplitsMax;

	//double mdCutoffMax;
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

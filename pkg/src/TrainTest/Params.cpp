/*
	ADRIAN DOBRA (ADOBRA@STAT.WASHINGTON.EDU)
*/

// Params.cpp: implementation of the Params class.
//
//////////////////////////////////////////////////////////////////////
#pragma warning(disable:4786)
#include "Params.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Params::Params()
{
	mnNumberOfGenes = 0;

	mstrTrainDataFile = "";
	mnTrainSampleSize = 0;
	
	mstrTestDataFile = "";
	mnTestSampleSize = 0;
	
	mnMaxRegressors = 0;
	mstrModelsFile = "";
	mnChainIterations = 0;
	mstrErrorMessage = "";
}

string Params::ToLower(string str)
{
	//new implementation for GNU
	char *newstr = strdup(str.c_str());
	int i = 0;
	while (newstr[i] != '\0') {
		newstr[i] = tolower(newstr[i]);
		i++;
	}
	return newstr;
	//return strlwr(strdup(str.c_str())); 
}

int Params::ToInt(string str)
{
	return atoi(str.c_str());
}

double Params::ToDouble(string str)
{
	return atof(str.c_str());
}

string Params::ToString(int value) {
	char  buffer[10]; 
	sprintf(buffer, "%d", value);
	string result(buffer);
	return result;
	cout << result.c_str() << endl;
	
}

string Params::ToString(double value) {
	char  buffer[10]; 
	sprintf(buffer, "%4.1f", value);
	string result(buffer);
	return result;
	
}

string Params::trim(string s)
{
	if (s.empty()) return s;
	string ret;
	for (int i = 0; i < s.length();  i++) {
		if (!isblank(s.at(i)) && !iscntrl(s.at(i)))
			ret.append(s.substr(i,1));
	}
	return ret;
}


bool Params::Load(string FileName){
	int BufferSize = 4096;
	char* theLine = new char[BufferSize];
	ifstream theFile(FileName.c_str());
	if (theFile.fail()) {
		mstrErrorMessage = "Failed to open the description file!";
		return false;
	}

	int nLineCount = 0;
	while (!theFile.eof()) {
		theFile.getline(theLine, BufferSize);
		nLineCount++;
		string theline(theLine);
		string Name(""), Value("");
		theline = trim(theline);
		if (theline.length() && (theline.c_str()[0] != '#'))
		{
			int pos = 0;
			if ((pos = theline.find("=")) != -1) {
				Name = theline.substr(0, pos);
				Value = theline.substr(pos + 1);
			}
			if (Name == "" && Value == "") {
			} else if (Name == "" || Value == "") {
				mstrErrorMessage = "Invalid <Name = Value> pair format";
				cout << theLine << endl;
				return false;
			} else {
				string name = ToLower(Name);
				string value = ToLower(Value);
				if (name == "numberofvariables") {
					mnNumberOfGenes = ToInt(value);
				} else if (name == "traindatafile") {
					mstrTrainDataFile = Value;
				} else if (name == "trainsamplesize") {
					mnTrainSampleSize = ToInt(value);
				} else if (name == "testdatafile") {
					mstrTestDataFile = Value;
				} else if (name == "testsamplesize") {
					mnTestSampleSize = ToInt(value);	
				} else if (name == "maxregressors") {
					mnMaxRegressors = ToInt(value);
				} else if (name == "modelsfile") {
					mstrModelsFile = Value;
				} else if (name == "chainiterations") {
					mnChainIterations = ToInt(value);										
				} else {
					mstrErrorMessage = "Unknown <Name = Value> pair!"; //to be refined later
					cout << theLine << endl;
					return false;
				}
				//cout << theLine << endl;	
			}
		}
	}
	delete[] theLine;
	mstrErrorMessage = "";
	return true;
}


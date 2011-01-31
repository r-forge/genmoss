/*
	ADRIAN DOBRA (ADOBRA@STAT.WASHINGTON.EDU)
*/

// Model.cpp: implementation of the Model class.
//
//////////////////////////////////////////////////////////////////////
#pragma warning(disable:4786)
#include "Model.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Model::Model()
{
	mnNumberOfVariables = 0;

	mstrTrainDataFile = "";	
	mnTrainSampleSize = 0;
	mstrTrainOutputFile = "";
	mstrScoresFile = "";
	//mstrTrainIndexFile = "";
	
	mstrTestDataFile = "";	
	mnTestSampleSize = 0;
	mstrTestOutputFile = "";
	//mstrTestIndexFile = "";
	
	mdSplitsMin = 0;
	mdSplitsInc = 0;
	mdSplitsMax = 0;
	//mdCutoffMax = 0;
	mstrErrorMessage = "";
}

string Model::ToLower(string str)
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

int Model::ToInt(string str)
{
	return atoi(str.c_str());
}

double Model::ToDouble(string str)
{
	return atof(str.c_str());
}

string Model::ToString(int value) {
	char  buffer[10]; 
	sprintf(buffer, "%d", value);
	string result(buffer);
	return result;
	cout << result.c_str() << endl;
	
}

string Model::ToString(double value) {
	char  buffer[10]; 
	sprintf(buffer, "%4.1f", value);
	string result(buffer);
	return result;
	
}

string Model::trim(string s)
{
	if (s.empty()) return s;
	string ret;
	for (int i = 0; i < s.length();  i++) {
		if (!isblank(s.at(i)) && !iscntrl(s.at(i)))
			ret.append(s.substr(i,1));
	}
	return ret;
}


bool Model::Load(string FileName){
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
					mnNumberOfVariables = ToInt(value);				
				} else if (name == "traindatafile") {
					mstrTrainDataFile = Value;
				} else if (name == "trainsamplesize") {
					mnTrainSampleSize = ToInt(value);
				} else if (name == "trainoutputfile") {
					mstrTrainOutputFile = Value;
				} else if (name == "scoresfile") {
					mstrScoresFile = Value;
				//} else if (name == "trainindexfile") {
				//	mstrTrainIndexFile = Value;
				} else if (name == "testdatafile") {
					mstrTestDataFile = Value;
				} else if (name == "testsamplesize") {
					mnTestSampleSize = ToInt(value);
				} else if (name == "testoutputfile") {
					mstrTestOutputFile = Value;
				//} else if (name == "testindexfile") {
				//	mstrTestIndexFile = Value;							
				} else if (name == "splitsmin") {
					mdSplitsMin = ToDouble(value);
				} else if (name == "splitsinc") {
					mdSplitsInc = ToDouble(value);
				} else if (name == "splitsmax") {
					mdSplitsMax = ToDouble(value);
				//} else if (name == "cutoffmax") {
				//	mdCutoffMax = ToDouble(value);					
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


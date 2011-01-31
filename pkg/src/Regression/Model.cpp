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
	mstrDataFile = "";
	mnNumberOfGenes = 0;
	mnSampleSize = 0;
	mnMaxRegressors = 0;
	//mnMaxRegressorsInteraction = 0;
	mnChainIterations = 0;
	//mnChainIterationsBurnin = 0;
	mnChainReplicates = 0;
	mdCutoffMax = 0.0;
	mdCutoffMin = 0.0;
	mdProbMax = 0.0;
	//mnDoInteractions = 0;
	//mnDoSpaceRatio = 0;
	//mnDoMetropolisH = 0;

	mnNumOfConfoundingVars = 0;
	
	//model search for tables
	mnShotgunChainReplicates = 0;
	mdShotgunCutoffMax = 0.0;
	mdShotgunCutoffMin = 0.0;
	mdShotgunProbMax = 0.0;
	
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
				if (name == "datafile") {
					mstrDataFile = Value;
				} else if (name == "numberofvariables") {
					mnNumberOfGenes = ToInt(value);
				} else if (name == "samplesize") {
					mnSampleSize = ToInt(value);
				} else if (name == "maxregressors") {
					mnMaxRegressors = ToInt(value);
				//} else if (name == "maxregressorsinteraction") {
				//	mnMaxRegressorsInteraction = ToInt(value);
				//} else if (name == "dointeractions") {
				//	mnDoInteractions = ToInt(value);	
				//} else if (name == "dospaceratio") {
				//	mnDoSpaceRatio = ToInt(value);
				//} else if (name == "dometropolish") {
				//	mnDoMetropolisH = ToInt(value);															
 				} else if (name == "numofconfoundingvars") {
                                        mnNumOfConfoundingVars = ToInt(value);
		 		} else if (name == "chainiterations") { 
                                        mnChainIterations = ToInt(value);
				//} else if (name == "chainiterationsburnin") {
				//	mnChainIterationsBurnin = ToInt(value);	
				} else if (name == "chainreplicates") {
					mnChainReplicates = ToInt(value);
				} else if (name == "cutoffmax") {
					mdCutoffMax = ToDouble(value);
				} else if (name == "cutoffmin") {
					mdCutoffMin = ToDouble(value);	
				} else if (name == "probmax") {
					mdProbMax = ToDouble(value);
				} else if (name == "shotgunchainreplicates") {
					mnShotgunChainReplicates = ToInt(value);
				} else if (name == "shotguncutoffmax") {
					mdShotgunCutoffMax = ToDouble(value);
				} else if (name == "shotguncutoffmin") {
					mdShotgunCutoffMin = ToDouble(value);	
				} else if (name == "shotgunprobmax") {
					mdShotgunProbMax = ToDouble(value);							
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


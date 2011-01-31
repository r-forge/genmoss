/*
AUTHOR: ADRIAN DOBRA (ADOBRA@STAT.WASHINGTON.EDU)
*/

// equations.h: interface for the CEquations class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(EQUATIONS_H)
#define EQUATIONS_H
class CEquations  
{
public:
	CEquations();
	virtual ~CEquations();

	short eqNumber;
	double Weight;
	CEquations* Next;

	void Delete(int id);
	void DeleteAll();
	void Add(int eqnumber,double weight);
	int harvest(int* eqids, double* weights);
	void UpdateScore(int chosenEquation,double* dagPost, double* CandidatePostA,double* CandidatePostB);
};

#endif

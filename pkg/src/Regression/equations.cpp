/*
AUTHOR: ADRIAN DOBRA (ADOBRA@STAT.WASHINGTON.EDU)
*/

// equations.cpp: implementation of the CEquations class.
//
//////////////////////////////////////////////////////////////////////
#include <string>
#include "equations.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CEquations::CEquations()
{
	Next = NULL;
}

CEquations::~CEquations()
{
	return;
}

void CEquations::Delete(int id)
{
	CEquations* p = this;
	CEquations* pnext;

	//find the equation with that id
	while(NULL!=p->Next) {
		if(p->Next->eqNumber!=id) {
			p = p->Next;
		} else {
			break;
		}
	}
	if(NULL!=p->Next) {
		pnext = p->Next;
		p->Next = pnext->Next;
		pnext->eqNumber = 0;
		pnext->Next = NULL;
		delete pnext;
	}
}

//deletes the entire list excluding the head
void CEquations::DeleteAll()
{
	CEquations* p = Next;
	CEquations* pnext;
	while(NULL!=p) {
		pnext = p->Next;
		p->eqNumber = 0;
		p->Weight = 0;
		p->Next = NULL;
		delete p;
		p = pnext;
	}
	Next = NULL;
}

void CEquations::Add(int eqnumber,double weight)
{
	CEquations* p = this;
	while(NULL!=p->Next) {
		if(p->Next->Weight<weight) {
			p = p->Next;
		} else { 
			break;
		}
	}

	CEquations* newp = new CEquations();
	newp->eqNumber = eqnumber;
	newp->Weight = weight;
	newp->Next = p->Next;
	p->Next = newp;
	return;
}

int CEquations::harvest(int* eqids, double* weights)
{
	int n = 0;
	CEquations* p  = Next;
	while(NULL!=p) {
		eqids[n] = p->eqNumber;
		weights[n] = p->Weight;
		n++;
		p = p->Next;
	}   
	return n;
}

void CEquations::UpdateScore(int chosenEquation,
			     double* dagPost, double* CandidatePostA,double* CandidatePostB) {
	Delete(chosenEquation);
	double delta = CandidatePostA[chosenEquation]+CandidatePostB[chosenEquation];
	delta -= (dagPost[chosenEquation]+dagPost[chosenEquation+1]);
	Add(chosenEquation,delta);
}

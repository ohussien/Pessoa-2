/***********************************************************************

	PESSOA Version 2.0  
	------------------

Cyber-Physical Systems Laboratory
http://www.cyphylab.ee.ucla.edu
Author(s):	Omar Hussien - ohussien@ucla.edu

Dependencies:	CUDD library, 
University of California, Los Angeles.
November 2015. 

************************************************************************/



#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "cudd.h"
#include "dddmp.h"
//#include <boost/numeric/odeint.hpp>

#include <vector>
#include <algorithm>
using namespace std;

//using namespace boost::numeric::odeint;


class System {
private:


	int n;							//  number of states
	int m;							//  number of inputs
	double tau;						// time discretization parameter
	double eta;						// state discretization parameter
	double mu;						// input discretization parameter
	double epsilon;					// model resoultion
	int* min;  						// xu minimum values wrt eta and mu 
	int* max;						// xu maximum values wrt eta and mu 
	int* num;						// xu max-min
	int* nume;						// xux' max-min
	int* nxbits;					// number of bits allocated to x
	int* nubits;					// number of bits allocated to u
	int nbitsx;						// total number of bits for x
	int nbitsloop;					// total number of bits for xu
	int* nbits;						// number of bits allocated to xux'
	int totbits;					// total number of bits for xux'
	int type;						// type of system, 1 for incrementally stable otherwise 2
	DdNode* trans;					// transition system, initially everything is mapped to 1/2
	//DdNode* targetset;				// target set, for now it is computed using pessoa 1.4 
	//DdNode* Controller;				// controller ADD
	DdNode* v;			// domain of controller ADD
	int refine_count;				// number of iterations for the refine block
public:
	System(int n_in, int m_in, double tau_in, double eta_in, double mu_in,double epsilon_in,double* dxmin, double* dumin,double* dxmax, double* dumax, int ref_cnt=1);
	~System();
//	void UpdateBlock(DdManager* ddman);
//	void RefineBlock(DdManager* ddman);
	void SynSafety1(DdManager* ddman,DdNode* targetset);
	void SynSafety2(DdManager* ddman,DdNode* targetset);
	void SynReach1(DdManager* ddman,DdNode* targetset);
	int SaveController(DdManager* ddman, string filename);
	int SaveParams(string filename);
	int GetTotbits();
	int CmpControllers(DdManager* ddman,DdNode* cont_old);
	DdNode* ComputeReach(DdManager* ddman,int* x_ref,int* u_ref);
};

#endif

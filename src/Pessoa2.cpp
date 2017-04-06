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
#include "System.h"

using namespace std;


static DdManager* ddman;     // global ddmanager pointer variable
int numroots_abstract=1;    // number of roots/BDDs/functions

int numroots_design=1; 

int main(int argc, char *argv[]) {

int n,m,i,j;
FILE *fp1,*fp2;
double tau,eta,mu,epsilon;
double* xdmin;
double* xdmax;
double* udmin;
double* udmax;
DdNode* targetset;
DdNode* control;
clock_t begin, end;
double time_spent;

begin = clock();
ifstream Input("input.txt");
	if (!Input.is_open()) {
		cerr << "Error opening file: input.txt\n";
		return 0;
	}
Input>>n>>m;
Input>>tau>>eta>>mu>>epsilon;
 xdmin=new double[n];
 xdmax=new double[n];
 udmin=new double[m];
 udmax=new double[m];


// printf("\n n=%d m=%d\n",n,m);
// printf("\n tau=%f eta=%f mu=%f epsilon=%f\n",tau,eta,mu,epsilon);

for(i=0;i<n;i++)
	Input>>xdmin[i]>>xdmax[i];
for(i=0;i<m;i++)
	Input>>udmin[i]>>udmax[i];
// for(i=0;i<n;i++){
// 	printf("xdmin[%d]=%f \t",i,xdmin[i]);
// 	printf("xdmax[%d]=%f \n",i,xdmax[i]);
// }
// for(i=0;i<m;i++){
// 	printf("udmin[%d]=%f \t",i,udmin[i]);
// 	printf("udmax[%d]=%f \n",i,udmax[i]);
// }
Input.close();
System* sys1=new System(n,m,tau,eta,mu,epsilon,xdmin,udmin,xdmax,udmax,5000);
delete xdmin;
delete xdmax;
delete udmin;
delete udmax;

//BDD Variables	
	short numVars;  // num of vars; if unknown set to 0
	short numVarsZ; // num of vars for ZBDDs; if unknown set to 0
	int numSlots; // default for CUDD package
	int cacheSize; // default for CUDD package
	//int maxCacheSize;   // default for CUDD package
numVars=sys1->GetTotbits();  // num of vars; if unknown set to 0
	numVarsZ=0; // num of vars for ZBDDs; if unknown set to 0
	numSlots=CUDD_UNIQUE_SLOTS; // default for CUDD package
	cacheSize=CUDD_CACHE_SLOTS; // default for CUDD package
	long maxCacheSize=0; //10485760*2;   // default for CUDD package

	ddman = Cudd_Init(numVars, numVarsZ, numSlots, cacheSize, maxCacheSize); //maxCacheSize);




sys1->SaveParams("out");

DdNode* vx = Cudd_ReadZero(ddman);
	Cudd_Ref(vx);
// fp1=fopen("v_test.dot", "w");
// 			Cudd_DumpDot(ddman, numroots_design, &vx,NULL,NULL,fp1); 
// 			fclose(fp1);
// Dddmp_cuddBddStore(ddman, NULL, vx, NULL,
// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
// 					"v_test.bdd", NULL);

Cudd_RecursiveDeref(ddman,vx);



//targetset=Dddmp_cuddBddLoad(ddman, DDDMP_VAR_MATCHIDS, NULL, NULL, NULL,DDDMP_MODE_DEFAULT, "target.bdd", NULL);
targetset=Dddmp_cuddBddLoad(ddman, DDDMP_VAR_MATCHIDS, NULL, NULL, NULL,DDDMP_MODE_DEFAULT, "Constraint.bdd", NULL);

printf("read target set done\n");

sys1->SynSafety1(ddman,targetset);
printf("safety  done\n");



//control=Dddmp_cuddBddLoad(ddman, DDDMP_VAR_MATCHIDS, NULL, NULL, NULL,DDDMP_MODE_DEFAULT, "IPController_dom.bdd", NULL);
control=Dddmp_cuddBddLoad(ddman, DDDMP_VAR_MATCHIDS, NULL, NULL, NULL,DDDMP_MODE_DEFAULT, "unicyle_control_dom.bdd", NULL);

Cudd_Ref(control);
printf("now comparing controllers\n");
int test=sys1->CmpControllers(ddman,control);

//Cudd_RecursiveDeref(ddman,control);
printf("now calling function to save controller \n");
sys1->SaveController(ddman,"cont1");
printf("saving controller file done\n");
delete sys1;


Cudd_Quit(ddman);

end = clock();
time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

printf("time to compute safety controller %f\n",time_spent);
return 0;


}

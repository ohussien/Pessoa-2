/***********************************************************************

	PESSOA Version 2.0  
	------------------

Cyber-Physical Systems Laboratory
http://www.cyphylab.ee.ucla.edu
Author(s):	Omar Hussien - ohussien@ucla.edu

Dependencies:	CUDD lib, System.h, 
University of California, Los Angeles.
November 2015. 

************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include "util.h"
#include "cudd.h"
#include "dddmp.h"
#include "System.h"


using namespace std;

System::System(int n_in, int m_in, double tau_in, double eta_in, double mu_in,double epsilon_in,double* dxmin, double* dumin,double* dxmax, double* dumax,int ref_cnt) {
	
	int i,j,k,temp;
	n=n_in;
	m=m_in;
	tau=tau_in;
	eta=eta_in;
	mu=mu_in;
	epsilon=epsilon_in;
	min= new int[n+m];
	max= new int[n+m];
	num= new int[n+m];
	nume= new int[2*n+m];
	nxbits= new int[n];
	nubits=new int[m];
	nbits=new int[2*n+m];
	for(i=0;i<n;i++)
	{
		min[i]=floor(dxmin[i]/eta);
		max[i]=ceil(dxmax[i]/eta);
		num[i]=max[i]-min[i];
		nume[i]=max[i]-min[i];
		nume[i+m+n]=max[i]-min[i];//nume[i];
		k=0;
		temp=num[i];
		if(temp !=0)
		while(temp>0)
		{
			k++;
			temp=temp>>1;
		}
		else
			k=1;
		nxbits[i]=k;
		nbits[i]=k;
		nbits[i+m+n]=k;

	}
	for(i=0;i<n;i++){
//	printf("min[%d]=%f \t",i,dxmin[i]);
//	printf("max[%d]=%f \n",i,dxmax[i]);
}
	j=0;
	for(i=n;i<n+m;i++)
	{
		min[i]=ceil(dumin[j]/mu);
		max[i]=floor(dumax[j]/mu);
		num[i]=max[i]-min[i];
		nume[i]=max[i]-min[i];
		k=0;
		temp=num[i];
		if(temp !=0)
		while(temp>0)
		{
			k++;
			temp=temp>>1;
		}
		else
			k=1;
		nubits[j]=k;
		nbits[j+n]=k;
		j++;
	}

	nbitsx=0;
	for(i=0;i<n;i++)
		nbitsx+=nxbits[i];
	nbitsloop=nbitsx;
	for(i=0;i<m;i++)
		nbitsloop+=nubits[i];

totbits=nbitsloop+nbitsx;


refine_count=ref_cnt;


// for(i=0;i<n;i++)
// printf("---------  nxbits[%d]=%d\n", i,nxbits[i]);

}

System::~System() {


	delete min; 
	delete max;
	delete num;
	delete nume;
	delete nxbits;
	delete nubits;
	delete nbits;
}




void System::SynSafety1(DdManager* ddman,DdNode* targetset){

	//// vars for debugging

	int cubes;
	int* testu;
	DdGen *gen11;
	int suc11,cubes11;
	double* test_x;
	/////  end of vars for debugging	

	int fixed=0;
	int break_refine=0;
	int refine_count_current;
	int numroots_design=1; 
	FILE *fp1,*fp2;
	double* u_oob;
	int *permutation, *existential, *existential2, *u_check, *boundu;
	int i,j,suc,k,jj;
	long count;
	DdNode *tmp1, *tmp2, *tmp3, *tmp4, *tmp5, *tmp6, *cons, *cons2, *cons2m, *trans_prev, *v_prev, *tcube_u, *tcube_xp;
	DdGen *gen;
	CUDD_VALUE_TYPE value;
	int *x_ref, *u_ref;
	unsigned short num_dc;
	trans=Cudd_ReadOne(ddman);
	Cudd_Ref(trans);

	u_oob=new double[m];

	test_x=new double[n];



	///////// debug code ////////////


	//printf("just saved v\n");

	// printf("$$$$$$ printing states inside targetset  at the beginning  $$$$ \n");
	// gen11=Cudd_FirstCube(ddman,targetset,&testu,&value);
	// if(!Cudd_IsGenEmpty(gen11))
	//     {
	//     	if(value>0.5)
	//     	{
	// for(k=0;k<n;k++){
	// 	test_x[k]=min[k];
	// 	for(jj=0;jj<nbits[k];jj++)
	// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
	// 	test_x[k]=test_x[k]*eta;
	// 	printf("x[%d]=%f\n",k,test_x[k]);
	// 	}
	// 	if(count==1)
	// 		for(jj=0;jj<nbitsx;jj++)
	// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
	// printf("value=%f\n",value);
	// }
	// }


	// for(cubes=0;cubes<30;cubes++)
	// if(!(Cudd_IsGenEmpty(gen11)))
	// {
	// suc11=Cudd_NextCube(gen11,&testu,&value);
	// if(value>0.5){
	// for(k=0;k<n;k++){
	// 	test_x[k]=min[k];
	// 	for(jj=0;jj<nbits[k];jj++)
	// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
	// 	test_x[k]=test_x[k]*eta;
	// 	printf("x[%d]=%f\n",k,test_x[k]);
	// 	}
	// 		if(count==1)
	// 			for(jj=0;jj<nbitsx;jj++)
	// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
	// printf("value=%f\n",value);
	// }
	// }

	////// end of debugging code ///////////


	////////////////////////////////  preprocessing the target set from {0,1} to {-1,0,1}  ///////////////




	tmp1=Cudd_BddToAdd(ddman,targetset);
	Cudd_Ref(tmp1);

	v = tmp1;

	trans_prev=Cudd_ReadOne(ddman);
	Cudd_Ref(trans_prev);

	v_prev=Cudd_ReadZero(ddman);
	Cudd_Ref(v_prev);

	////////////////////////////////////////////////////////////////////////////////////////
	u_check= new int[totbits];

	for (i=totbits-nbitsloop; i<totbits; i++)
		u_check[i] = 2;	

	//Build for permutation of x,u,x' to x',u,x
	permutation =  new int[totbits];
	
	// nbitsloop=bits(x)+bits(u)

	//initialize u
	for (i=totbits-nbitsloop; i<nbitsloop; i++)
		permutation[i] = i;

	//reorder initial x (u remains the same)
	for (i=0; i<totbits-nbitsloop; i++)
		permutation[i] = nbitsloop+i;

	//reorder final x
	for (i=nbitsloop; i<totbits; i++)
		permutation[i] = i-nbitsloop;	

	////////////////////////////////////////////////////////////////////////////////////////

	// Build cube for existential check of u and x'
	existential = new int[totbits];

	//initialize everything to 2
	for (i=0; i<totbits-nbitsloop; i++)
		existential[i] = 2;

	//final x and u set to 1 -- rest is ignored
	for (i=totbits-nbitsloop; i<totbits; i++)
		existential[i] = 1;	

	// another existential for x' only
	existential2 = new int[totbits];

	//initialize everything to 2 
	for (i=0; i<totbits-nbitsx; i++)
		existential2[i] = 2;

	//final x to 1 -- rest is ignored
	for (i=totbits-nbitsx; i<totbits; i++)
		existential2[i] = 1;

	tmp1 = Cudd_CubeArrayToBdd(ddman, existential);
	Cudd_Ref(tmp1);

	tcube_u =Cudd_BddToAdd(ddman,tmp1);
	Cudd_Ref(tcube_u);

	Cudd_RecursiveDeref(ddman,tmp1);

	tmp1 = Cudd_CubeArrayToBdd(ddman, existential2);
	Cudd_Ref(tmp1);

	tcube_xp =Cudd_BddToAdd(ddman,tmp1);
	Cudd_Ref(tcube_xp);

	Cudd_RecursiveDeref(ddman,tmp1);

	////////////////////////////////////////////////////////////////////////////////////////


	/////////  Cons is used to get the negation of T(x,u,x') //////////////////////////

	cons=Cudd_addConst(ddman,1);
	Cudd_Ref(cons);

	cons2=Cudd_addConst(ddman,0.5);
	Cudd_Ref(cons2);

	cons2m=Cudd_addConst(ddman,-0.5);
	Cudd_Ref(cons2m);


	tmp1=Cudd_addApply(ddman,Cudd_addTimes,trans,cons2);
	Cudd_Ref(tmp1);

	Cudd_RecursiveDeref(ddman,trans);

	trans=tmp1;


 	suc=nbitsx;

 	boundu = new int[totbits];
 	for(i=0;i<m;i++)
 	{

 		for(j=num[n+i]+1;j<(1<<nbits[n+i]);j++)
 		{

 			for (k=0; k<totbits; k++)
 				boundu[k]=2;

			// construct the out of bound u
			for(k=0;k<nbits[n+i];k++)
				boundu[suc+k]=(j>>(nbits[n+i]-k-1))& 1;

			// printf("for input %d boundu[%d]= ",i,j);
			// for(k=0;k<totbits;k++)
			// 	printf("%d",boundu[k]);
			// printf("\n");


			// convert to bdd using cubearraytobdd function
			tmp1 = Cudd_CubeArrayToBdd(ddman, boundu);
			Cudd_Ref(tmp1);
			// invert then convert to ADD	
			tmp2=Cudd_BddToAdd(ddman,tmp1);
			Cudd_Ref(tmp2);

			Cudd_RecursiveDeref(ddman,tmp1);

			tmp1=Cudd_addApply(ddman,Cudd_addTimes,cons2,tmp2);
			Cudd_Ref(tmp1);

			Cudd_RecursiveDeref(ddman,tmp2);

			// multiply with trans to get rid of this input

			tmp2=Cudd_addApply(ddman,Cudd_addPlus,tmp1,trans);
			Cudd_Ref(tmp2);

			Cudd_RecursiveDeref(ddman,tmp1);
			Cudd_RecursiveDeref(ddman,trans);

			trans=tmp2;

 		}
 		suc+=nbits[n+i];	
 	}

	delete boundu;

	count=0;






	do
	{
		//printf("Update started\n");
		count++;
		refine_count_current=refine_count;

	

		Cudd_RecursiveDeref(ddman, trans_prev);
		trans_prev=trans;

		Cudd_RecursiveDeref(ddman, v_prev);

		tmp1=Cudd_addBddInterval(ddman,v, 0.5,1);
		Cudd_Ref(tmp1);

		Cudd_RecursiveDeref(ddman,v);



		v_prev=Cudd_BddToAdd(ddman,tmp1);
		Cudd_Ref(v_prev);

		Cudd_RecursiveDeref(ddman,tmp1);

		//printf("start the loop\n");

		////////// debug code



		// printf("----------------------------------------- printing states inside v_prev at the begining of iteration %d !!!!!!!!!!! \n",count);
		// gen11=Cudd_FirstCube(ddman,v_prev,&testu,&value);
		// if(!Cudd_IsGenEmpty(gen11))
		//     {
		// for(jj=0;jj<totbits;jj++)
		//     printf("testu[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }


		// for(cubes=0;cubes<6;cubes++)
		// if(!(Cudd_IsGenEmpty(gen11)))
		// {
		// suc11=Cudd_NextCube(gen11,&testu,&value);
		// for(jj=0;jj<totbits;jj++)
		//     printf("testu[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		//}

		////// end of debugging code ///////////


		///////////////////  Update V ///////////////////


		//printf("just started update V\n");
		tmp1=Cudd_addPermute(ddman, v_prev, permutation);
		Cudd_Ref(tmp1);

		tmp2=Cudd_addApply(ddman,Cudd_addMinimum,trans_prev,tmp1);
		Cudd_Ref(tmp2);

		Cudd_RecursiveDeref(ddman,tmp1);


		tmp3=Cudd_addApply(ddman,Cudd_addMinus,cons,trans_prev);
		Cudd_Ref(tmp3);


		tmp4=Cudd_addApply(ddman,Cudd_addMaximum,tmp2,tmp3);
		Cudd_Ref(tmp4);

		Cudd_RecursiveDeref(ddman,tmp2);
		Cudd_RecursiveDeref(ddman,tmp3);

		tmp5=Cudd_addUnivAbstract(ddman,tmp4,tcube_xp);
		Cudd_Ref(tmp5);

		Cudd_RecursiveDeref(ddman,tmp4);

		tmp6=Cudd_addExistAbstract(ddman,tmp5,tcube_u);
		Cudd_Ref(tmp6);

		Cudd_RecursiveDeref(ddman,tmp5);

		v=Cudd_addApply(ddman,Cudd_addMinimum,v_prev,tmp6);
		Cudd_Ref(v);

		Cudd_RecursiveDeref(ddman,tmp6);

		//printf("Update done\n");


		/////////////////  End of Update V ///////////////////



		///////// debug code ////////////



		// printf("$$$$$$ printing states inside V_prev  at iteration %d $$$$ \n",count);
		// gen11=Cudd_FirstCube(ddman,v_prev,&testu,&value);
		// if(!Cudd_IsGenEmpty(gen11))
		//     {
		//     	if(value>0.2)
		//     	{

		// 	for(k=0;k<n;k++){
		// 	test_x[k]=min[k];
		// 	for(jj=0;jj<nbits[k];jj++)
		// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
		// 	test_x[k]=test_x[k]*eta;
		// 	printf("x[%d]=%f\n",k,test_x[k]);
		// 	}
		    
		// printf("value=%f\n",value);
		// }
		// for(jj=0;jj<nbitsx;jj++)
		// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }


		// for(cubes=0;cubes<30;cubes++)
		// if(!(Cudd_IsGenEmpty(gen11)))
		// {
		// suc11=Cudd_NextCube(gen11,&testu,&value);
		// if(value>0.2){
		// for(k=0;k<n;k++){
		// 	test_x[k]=min[k];
		// 	for(jj=0;jj<nbits[k];jj++)
		// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
		// 	test_x[k]=test_x[k]*eta;
		// 	printf("x[%d]=%f\n",k,test_x[k]);
		// 	}
		// printf("value=%f\n",value);
		// }
		// for(jj=0;jj<nbitsx;jj++)
		// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }









		// printf("$$$$$$ printing states inside V  at iteration %d $$$$ \n",count);
		// gen11=Cudd_FirstCube(ddman,v,&testu,&value);
		// if(!Cudd_IsGenEmpty(gen11))
		//     {
		//     	if(value>0.2)
		//     	{

		// 	for(k=0;k<n;k++){
		// 	test_x[k]=min[k];
		// 	for(jj=0;jj<nbits[k];jj++)
		// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
		// 	test_x[k]=test_x[k]*eta;
		// 	printf("x[%d]=%f\n",k,test_x[k]);
		// 	}
		    
		// printf("value=%f\n",value);
		// }
		// for(jj=0;jj<nbitsx;jj++)
		// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }


		// for(cubes=0;cubes<30;cubes++)
		// if(!(Cudd_IsGenEmpty(gen11)))
		// {
		// suc11=Cudd_NextCube(gen11,&testu,&value);
		// if(value>0.2){
		// for(k=0;k<n;k++){
		// 	test_x[k]=min[k];
		// 	for(jj=0;jj<nbits[k];jj++)
		// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
		// 	test_x[k]=test_x[k]*eta;
		// 	printf("x[%d]=%f\n",k,test_x[k]);
		// 	}
		// printf("value=%f\n",value);
		// }
		// for(jj=0;jj<nbitsx;jj++)
		// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }




			////// end of debugging code ///////////




		if(Cudd_addLeq(ddman,v,v_prev) && Cudd_addLeq(ddman,v_prev,v))
		{
			printf("now v and v_prev are the same\n");
			fixed=1;
			break;
		}


		// printf("just finished update v in iteration number %d\n",count);


		// 	///// debugging code to be removed /////
		// if(count==1)
		// {
		// fp1=fopen("v1_test.dot", "w");
		// 			Cudd_DumpDot(ddman, numroots_design, &v,NULL,NULL,fp1); 
		// 			fclose(fp1);
		// Dddmp_cuddAddStore(ddman, NULL, v, NULL,
		// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
		// 					"v1_test.bdd", NULL);
		// }
		// if(count==2)
		// {
		// fp1=fopen("v2_test.dot", "w");
		// 			Cudd_DumpDot(ddman, numroots_design, &v,NULL,NULL,fp1); 
		// 			fclose(fp1);
		// Dddmp_cuddAddStore(ddman, NULL, v, NULL,
		// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
		// 					"v2_test.bdd", NULL);
		// }

		// printf("####################################### printing states inside V  before changing from -1 0 1 for iteration number %d \n",count);
		// gen11=Cudd_FirstCube(ddman,v,&testu,&value);
		// if(!Cudd_IsGenEmpty(gen11))
		//     {
		// for(jj=0;jj<totbits;jj++)
		//     printf("testu[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }


		// for(cubes=0;cubes<6;cubes++)
		// if(!(Cudd_IsGenEmpty(gen11)))
		// {
		// suc11=Cudd_NextCube(gen11,&testu,&value);
		// for(jj=0;jj<totbits;jj++)
		//     printf("testu[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }

			////// end of debugging code ///////////

		/// old code due to using -1 0 1

		tmp1=Cudd_addApply(ddman,Cudd_addMinus,v,cons);
		Cudd_Ref(tmp1);

		tmp2=Cudd_addApply(ddman,Cudd_addTimes,tmp1,cons2);
		Cudd_Ref(tmp2);

		Cudd_RecursiveDeref(ddman,tmp1);

		// end of old code

		//printf("about to get the set of states need to be discovered\n");


		////////// debug code



		// printf("####  printing states 1 1/2 0 inside V  at iteration %d ######## \n",count);
		// gen11=Cudd_FirstCube(ddman,v,&testu,&value);
		// if(!Cudd_IsGenEmpty(gen11))
		//     {
		//     	if(value>0.2)
		//     	{

		// 	for(k=0;k<n;k++){
		// 	test_x[k]=min[k];
		// 	for(jj=0;jj<nbits[k];jj++)
		// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
		// 	test_x[k]=test_x[k]*eta;
		// 	printf("x[%d]=%f\n",k,test_x[k]);
		// 	}
		    
		// printf("value=%f\n",value);
		// }
		// for(jj=0;jj<nbitsx;jj++)
		// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }


		// for(cubes=0;cubes<30;cubes++)
		// if(!(Cudd_IsGenEmpty(gen11)))
		// {
		// suc11=Cudd_NextCube(gen11,&testu,&value);
		// if(value>0.2){
		// for(k=0;k<n;k++){
		// 	test_x[k]=min[k];
		// 	for(jj=0;jj<nbits[k];jj++)
		// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
		// 	test_x[k]=test_x[k]*eta;
		// 	printf("x[%d]=%f\n",k,test_x[k]);
		// 	}
		// printf("value=%f\n",value);
		// }
		// for(jj=0;jj<nbitsx;jj++)
		// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }





		////// end of debugging code ///////////

		int num_cubes=0;
		vector<vector<int> > x_ref_new;
		int tmp_num_dc;
		int test_endof_cubes;
		num_dc=0;

		gen=Cudd_FirstCube(ddman,v,&x_ref,&value);
		if(Cudd_IsGenEmpty(gen))
		{
			//	printf("breaking because no more states to verify is empty 1st break\n");
			refine_count_current=0;
		}
		while(1)    //value != 0.5
		{
			//	printf("######getting states that needs to be explored, so far we have %d states after discovering %d cubes####\n",num_dc,num_cubes);
			if(num_dc>=refine_count)
			{
			//printf("breaking because we got the states needed to be discovered\n");
			break;
			}
			
			if(value != 0.5)
			{
				if(Cudd_IsGenEmpty(gen))
				{
					//printf("breaking because no more states to verify is empty 1st break\n");
		    		refine_count_current=0;
					break;
				}
				suc=Cudd_NextCube(gen,&x_ref,&value);
				//printf("we have the next cube now!!!!\n");
			}
			if(value==0.5)
			{
				test_endof_cubes=1;
				for(i=0;i<nbitsx;i++)
					if(x_ref[i]!=2)
						test_endof_cubes=0;
				if(test_endof_cubes==0)
				{
					// 	printf("currentcube");
					// 		for(i=0;i<nbitsx;i++)
					// 			{
					// printf("%d",x_ref[i]);//x_ref_new[i][j];
					// 			}
					// 		printf("\n");
					tmp_num_dc=1;
		       		x_ref_new.push_back(vector<int>());  
		        	for (i=0;i<nbitsx;i++)
					{	
		     			x_ref_new[num_cubes].push_back(x_ref[i]);
		     			if(x_ref[i]==2)
		         			tmp_num_dc=tmp_num_dc<<1;
		   			}
		   			num_dc+= tmp_num_dc;
		   			num_cubes++;
		   		}
		   	}
		  value=0;

		}
       	
       
   
		if(num_cubes==0)
	    {
	   			// 	printf("breaking because no more states to verify is empty 2nd break\n");
	    		refine_count_current=0;
	   	}




     
			//printf("--------------------- number of states =1/2  is  %d at iteration %d ---------------------\n",num_dc,count);


		// for(i=0;i<x_ref_new.size();i++)
		// 	{
		// 		printf("x_ref_new[%d]=",i);
		// 	for(j=0;j<nbitsx;j++)
		// 	{
		// 			printf("%d",x_ref_new[i][j]);
		// 	}
		// 	printf("\n");
		// 	}

		i=0;
     	refine_count_current=num_dc;

     	if(num_dc>refine_count)
     	{
     		//num_dc=refine_count;
     		refine_count_current=refine_count;
		 }

		//	printf("num_dc=%d\n",num_dc);
   		while((num_dc>1)  &&(i<x_ref_new.size())&& (num_cubes !=0))//
   		{
        	//printf("after expanding state %d  and now num_dc=%d\n",i,num_dc);
     		for (j=0;j<nbitsx;j++)
   			{
     			if(x_ref_new[i][j]==2)
     				{
			         for(k=0;k<nbitsx;k++)
			             x_ref[k]=x_ref_new[i][k];
			         x_ref_new.erase(x_ref_new.begin()+i);
			         x_ref_new.push_back(vector<int>());
			         x_ref_new.push_back(vector<int>());
			         i--;
			         break;
     				}
         
   			}
     		for(k=0;k<nbitsx;k++)
     		{
         		if(k==j)
         		{
           			x_ref_new[x_ref_new.size()-2].push_back(0);
           			x_ref_new[x_ref_new.size()-1].push_back(1);
         		}
         		else
         		{
            		x_ref_new[x_ref_new.size()-2].push_back(x_ref[k]);
            		x_ref_new[x_ref_new.size()-1].push_back(x_ref[k]);
         		}
     		}
    
			// for(k=0;k<x_ref_new.size();k++)
			// 	{
			// 		printf("x_ref_new[%d]=",k);
			// 	for(j=0;j<nbitsx;j++)
			// 	{
			// 			printf("%d",x_ref_new[k][j]);
			// 	}
			// 	printf("\n");
			// 	}
	     
	   		//  num_dc--;
     		i++;
    	}	

		for(i=0;i<x_ref_new.size();i++)
			for(j=0;j<nbitsx;j++)
				if(x_ref_new[i][j]==2)
					x_ref_new[i][j]=0;

		Cudd_RecursiveDeref(ddman,tmp2);

		///// debug code




		// printf("after expanding the states \n");
		// for(i=0;i<x_ref_new.size();i++)
		// 	{
		// 		printf("x_ref_new[%d]=",i);
		// 	for(j=0;j<nbitsx;j++)
		// 	{
		// 			printf("%d",x_ref_new[i][j]);
		// 	}
		// 	printf("\n");
		// 	}



		////// end debug ///////

		/// old code due to using -1 0 1 in Trans

		// change T(x,u,x') values from {1,0,-1} to {1,1/2,0}

		//tmp1=Cudd_addApply(ddman,Cudd_addMinus,trans_prev,cons);
		//Cudd_Ref(tmp1);

		//Cudd_RecursiveDeref(ddman,trans_prev);

		//trans_prev=Cudd_addApply(ddman,Cudd_addTimes,tmp1,cons2);
		//Cudd_Ref(trans_prev);

		//Cudd_RecursiveDeref(ddman,tmp1);

		//// end of old code

		//printf("just started refine T\n");



		//// debug ////////


		// if(count==1)
		// {
		// fp1=fopen("trans_prev1_test.dot", "w");
		// 			Cudd_DumpDot(ddman, numroots_design, &trans_prev,NULL,NULL,fp1); 
		// 			fclose(fp1);
		// Dddmp_cuddAddStore(ddman, NULL, trans_prev, NULL,
		// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
		// 					"trans_prev1_test.bdd", NULL);
		// }
		// if(count==2)
		// {
		// fp1=fopen("trans_prev2_test.dot", "w");
		// 			Cudd_DumpDot(ddman, numroots_design, &trans_prev,NULL,NULL,fp1); 
		// 			fclose(fp1);
		// Dddmp_cuddAddStore(ddman, NULL, trans_prev, NULL,
		// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
		// 					"trans_prev2_test.bdd", NULL);
		// }

		////// end debug ///////



		///////////////////  Refine T ///////////////////
		//printf("Refine started\n");

		for(i=0;i<refine_count_current;i++)
		{

			//printf("Refine started, this is state number %d from %d states\n",i+1,refine_count_current);

			if(i>0)
			trans_prev=trans;
			break_refine=0;
			for (j=0; j<totbits-nbitsloop; j++)
				u_check[j]=x_ref_new[i][j];



			tmp1 = Cudd_CubeArrayToBdd(ddman, u_check);
			Cudd_Ref(tmp1);

			tmp2 =Cudd_BddToAdd(ddman,tmp1);
			Cudd_Ref(tmp2);

			Cudd_RecursiveDeref(ddman,tmp1);

			tmp1=Cudd_addApply(ddman,Cudd_addTimes,trans_prev,tmp2);
			Cudd_Ref(tmp1);

			Cudd_RecursiveDeref(ddman,tmp2);


			gen=Cudd_FirstCube(ddman,tmp1,&u_ref,&value);
			if(Cudd_IsGenEmpty(gen))
		    {
		    	break_refine=1;
		    }
			else
				while(value != 0.5 && break_refine !=1)
				{
					if(Cudd_IsGenEmpty(gen))
			    	{
			    		break_refine=1;
			    	}
			    	else
						suc=Cudd_NextCube(gen,&u_ref,&value);

				}
			Cudd_RecursiveDeref(ddman,tmp1);

			// for(k=0;k<totbits;k++)
			// 	printf("u_ref[%d]=%d \n",k,u_ref[k]);

			// if(count>0)
			// {	printf("break refine = %d\n",break_refine);
			// for(k=totbits-nbitsloop;k<totbits-nbitsx;k++)
			// 	printf("picked u to refine is u[%d]=%d \n",k,u_ref[k]);
			// }
			//call function to compute post(x_ref,u_ref)
			// note u_ref might have dont cares   -----fixed 
			for (k=totbits-nbitsloop; k<totbits-nbitsx; k++)
				if(u_ref[k] == 2)
					u_ref[k]=0;
			j=0;
			for(k=0;k<m;k++)
			{
				u_oob[k]=min[k+n];
				for(jj=0;jj<nbits[k+n];jj++)
				u_oob[k]+=u_ref[j+nbitsx+jj]<<(nbits[k+n]-jj-1);
				j+=nbits[k+n];
				u_oob[k]=u_oob[k]*mu;
			}

			for(k=0;k<m;k++)
			{
				//	printf("u_oob[%d]=%f \n",k,u_oob[k]);
				if(u_oob[k]>max[k]*mu)
				{
					
					//	printf("u picked outside the range of u\n");
					break_refine=1;
				}				
			}

			///// commented for now, instead of breaking don't call computereach **************
			//if(break_refine==1)
			//	break;



			// nullifying all transtions from this state for now.
			tmp1 = Cudd_CubeArrayToBdd(ddman, u_ref);
			Cudd_Ref(tmp1);

			// if(count>8)
			// {	printf("break refine = %d\n",break_refine);
			// for(k=0;k<totbits;k++)
			// 	printf("picked u to nullify is u[%d]=%d \n",k,u_ref[k]);
			// }

			tmp2 =Cudd_BddToAdd(ddman,tmp1);
			Cudd_Ref(tmp2);

			Cudd_RecursiveDeref(ddman,tmp1);

			tmp1=Cudd_addApply(ddman,Cudd_addTimes,tmp2,cons2m);
			Cudd_Ref(tmp1);
			Cudd_RecursiveDeref(ddman,tmp2);

			if(break_refine==0)
			{
				tmp2=Cudd_addApply(ddman,Cudd_addPlus,trans_prev,tmp1);
				Cudd_Ref(tmp2);
			}
			else
			{
				tmp2=Cudd_addApply(ddman,Cudd_addMinus,trans_prev,tmp1);
				Cudd_Ref(tmp2);
			}

			Cudd_RecursiveDeref(ddman,tmp1);
			//Cudd_RecursiveDeref(ddman,trans_prev);


			//printf("calling computereach for in iteration number %d \n",count);
			for(k=0;k<totbits-nbitsloop;k++)
			x_ref[k]=x_ref_new[i][k];

			if(break_refine==0)
				tmp1=ComputeReach(ddman,x_ref,u_ref);
			else
			{
				tmp1=Cudd_ReadZero(ddman);
				Cudd_Ref(tmp1);
			}
			//Cudd_Ref(tmp1);   already ref inside the function....


			trans=Cudd_addApply(ddman,Cudd_addPlus,tmp2,tmp1);
			Cudd_Ref(trans);

			Cudd_RecursiveDeref(ddman,tmp1);
			Cudd_RecursiveDeref(ddman,tmp2);



		}
		if(refine_count_current==0)
			trans=trans_prev;

		//printf("Refine done\n");
		///////////////////  End of Refine T ///////////////////	




		//// debugging code

		// if(count==1)
		// {
		// fp1=fopen("t1_test.dot", "w");
		// 			Cudd_DumpDot(ddman, numroots_design, &trans,NULL,NULL,fp1); 
		// 			fclose(fp1);
		// Dddmp_cuddAddStore(ddman, NULL, trans, NULL,
		// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
		// 					"t1_test.bdd", NULL);
		// }
		// if(count==2)
		// {
		// fp1=fopen("t2_test.dot", "w");
		// 			Cudd_DumpDot(ddman, numroots_design, &trans,NULL,NULL,fp1); 
		// 			fclose(fp1);
		// Dddmp_cuddAddStore(ddman, NULL, trans, NULL,
		// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
		// 					"t2_test.bdd", NULL);
		// }


		/////end of debugging code


		// change T(x,u,x') values from  {1,1/2,0}  to {1,0,-1} 

		/// old code due to using -1 0 1

		//tmp1=Cudd_addApply(ddman,Cudd_addDivide,trans,cons2);
		//Cudd_Ref(tmp1);

		//Cudd_RecursiveDeref(ddman,trans);

		//trans=Cudd_addApply(ddman,Cudd_addPlus,tmp1,cons);
		//Cudd_Ref(trans);

		//Cudd_RecursiveDeref(ddman,tmp1);


		/// end of old code

		//////// test iteration

		//	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~We are in iteration %d and comparison leads to  %d %d~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",count,Cudd_addLeq(ddman,v,v_prev) ,Cudd_addLeq(ddman,v_prev,v));
		/////debug code ///////////
		if(count>0)
		{	




			///////// debug code ////////////


			// //printf("just saved v\n");

			// printf("$$$$$$ printing trans inside trans_prev  at iteration %d $$$$ \n",count);
			// gen11=Cudd_FirstCube(ddman,trans_prev,&testu,&value);
			// if(!Cudd_IsGenEmpty(gen11))
			//     {
			//     	if(value>0.2)
			//     	{
			// 	for(jj=0;jj<totbits;jj++)
			//     printf("testu[%d]=%d\n",jj,testu[jj]);
			// printf("value=%f\n",value);
			// }
			// }


			// for(cubes=0;cubes<30;cubes++)
			// if(!(Cudd_IsGenEmpty(gen11)))
			// {
			// suc11=Cudd_NextCube(gen11,&testu,&value);
			// if(value>0.2){
			// 	for(jj=0;jj<totbits;jj++)
			//     printf("testu[%d]=%d\n",jj,testu[jj]);
			// printf("value=%f\n",value);
			// }
			// }

			// 	////// end of debugging code ///////////



			// tmp1=Cudd_addApply(ddman,Cudd_addMinus,trans,cons);
			// Cudd_Ref(tmp1);
			// printf("$$$$$$ printing trans inside trans  at iteration %d $$$$ \n",count);
			// gen11=Cudd_FirstCube(ddman,tmp1,&testu,&value);
			// if(!Cudd_IsGenEmpty(gen11))
			//     {
			//     	if(value>0.5)
			//     	{
			// 	for(jj=0;jj<totbits;jj++)
			//     printf("testu[%d]=%d\n",jj,testu[jj]);
			// printf("value=%f\n",value);
			// }
			// }


			// for(cubes=0;cubes<30;cubes++)
			// if(!(Cudd_IsGenEmpty(gen11)))
			// {
			// suc11=Cudd_NextCube(gen11,&testu,&value);
			// if(value>0.5){


			// 	for(jj=0;jj<totbits;jj++)
			//     printf("testu[%d]=%d\n",jj,testu[jj]);
			// printf("value=%f\n",value);


			// // for(k=0;k<n;k++){
			// // 	test_x[k]=min[k];
			// // 	for(jj=0;jj<nbits[k];jj++)
			// // 		test_x[k]+=testu[k*nbits[0]+jj]<<jj;
			// // 	test_x[k]=test_x[k]*eta;
			// // 	printf("x[%d]=%f\n",k,test_x[k]);
			// // 	}
			// // 	for(k=0;k<m;k++){
			// // 	u_oob[k]=min[k+n];
			// // 	for(jj=0;jj<nbits[k+n];jj++)
			// // 		u_oob[k]+=testu[nbitsx+jj]<<jj;
			// // 	u_oob[k]=u_oob[k]*mu;
			// // 	printf("u[%d]=%f\n",k,u_oob[k]);
			// // }
			// // printf("value=%f\n",value);
			// }
			// }
			// Cudd_RecursiveDeref(ddman,tmp1);



				////// end of debugging code ///////////



		}





		////// end of debug code //////////





		if(Cudd_addLeq(ddman,v,v_prev) && Cudd_addLeq(ddman,v_prev,v))// && Cudd_addLeq(ddman,trans,trans_prev) && Cudd_addLeq(ddman,trans_prev,trans))
		{	
			fixed=1;
			printf("number of iterations = %d so far *****\n",count);
			//printf("number of iterations = %d\n",count);
		}		
	}while(fixed != 1);
	printf("number of iterations = %d and fixed = %d\n",count,fixed);

	// 	fp1=fopen("v_test.dot", "w");
	// 			Cudd_DumpDot(ddman, numroots_design, &v,NULL,NULL,fp1); 
	// 			fclose(fp1);
	// Dddmp_cuddAddStore(ddman, NULL, v, NULL,
	// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
	// 					"v_test.bdd", NULL);

	////////// debug code



	// printf("####################################### printing states inside V  finally !!!!!!!!!!! \n");



	// gen11=Cudd_FirstCube(ddman,v,&testu,&value);
	// if(!Cudd_IsGenEmpty(gen11))
	//     {
	//     	if(value>0.5)
	//     	{

	// 	for(k=0;k<n;k++){
	// 	test_x[k]=min[k];
	// 	for(jj=0;jj<nbits[k];jj++)
	// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
	// 	test_x[k]=test_x[k]*eta;
	// 	printf("x[%d]=%f\n",k,test_x[k]);
	// 	}
	    
	// printf("value=%f\n",value);
	// }
	// for(jj=0;jj<nbitsx;jj++)
	// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
	// printf("value=%f\n",value);
	// }


	// for(cubes=0;cubes<200;cubes++)
	// if(!(Cudd_IsGenEmpty(gen11)))
	// {
	// suc11=Cudd_NextCube(gen11,&testu,&value);
	// if(value>0.5){
	// for(k=0;k<n;k++){
	// 	test_x[k]=min[k];
	// 	for(jj=0;jj<nbits[k];jj++)
	// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
	// 	test_x[k]=test_x[k]*eta;
	// 	printf("x[%d]=%f\n",k,test_x[k]);
	// 	}
	// printf("value=%f\n",value);
	// }
	// for(jj=0;jj<nbitsx;jj++)
	// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
	// printf("value=%f\n",value);
	// }



	////// end of debugging code ///////////



	//printf("just finished safety1\n");

	delete u_oob;
	delete permutation;
	delete existential;
	delete existential2;
	delete u_check;

	delete test_x;

	Cudd_RecursiveDeref(ddman,cons);
	Cudd_RecursiveDeref(ddman,cons2);
	Cudd_RecursiveDeref(ddman,cons2m);
	Cudd_RecursiveDeref(ddman,v_prev);
	Cudd_RecursiveDeref(ddman,trans_prev);
	printf("just finished safety1\n");

}

void System::SynSafety2(DdManager* ddman,DdNode* targetset){






}

void System::SynReach1(DdManager* ddman,DdNode* targetset)
{

	//// vars for debugging

	int cubes;
	int* testu;
	DdGen *gen11;
	int suc11,cubes11;
	double* test_x;
	/////  end of vars for debugging	

	int fixed=0;
	int break_refine=0;
	int refine_count_current;
	int numroots_design=1; 
	FILE *fp1,*fp2;
	double* u_oob;
	int *permutation, *existential, *existential2, *u_check, *boundu;
	int i,j,suc,k,jj;
	long count;
	DdNode *tmp1, *tmp2, *tmp3, *tmp4, *tmp5, *tmp6, *cons, *cons2, *cons2m, *trans_prev, *v_prev, *tcube_u, *tcube_xp;
	DdGen *gen;
	CUDD_VALUE_TYPE value;
	int *x_ref, *u_ref;
	unsigned short num_dc;
	trans=Cudd_ReadOne(ddman);
	Cudd_Ref(trans);

	u_oob=new double[m];

	test_x=new double[n];



	///////// debug code ////////////


	//printf("just saved v\n");

	// printf("$$$$$$ printing states inside targetset  at the beginning  $$$$ \n");
	// gen11=Cudd_FirstCube(ddman,targetset,&testu,&value);
	// if(!Cudd_IsGenEmpty(gen11))
	//     {
	//     	if(value>0.5)
	//     	{
	// for(k=0;k<n;k++){
	// 	test_x[k]=min[k];
	// 	for(jj=0;jj<nbits[k];jj++)
	// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
	// 	test_x[k]=test_x[k]*eta;
	// 	printf("x[%d]=%f\n",k,test_x[k]);
	// 	}
	// 	if(count==1)
	// 		for(jj=0;jj<nbitsx;jj++)
	// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
	// printf("value=%f\n",value);
	// }
	// }


	// for(cubes=0;cubes<30;cubes++)
	// if(!(Cudd_IsGenEmpty(gen11)))
	// {
	// suc11=Cudd_NextCube(gen11,&testu,&value);
	// if(value>0.5){
	// for(k=0;k<n;k++){
	// 	test_x[k]=min[k];
	// 	for(jj=0;jj<nbits[k];jj++)
	// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
	// 	test_x[k]=test_x[k]*eta;
	// 	printf("x[%d]=%f\n",k,test_x[k]);
	// 	}
	// 		if(count==1)
	// 			for(jj=0;jj<nbitsx;jj++)
	// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
	// printf("value=%f\n",value);
	// }
	// }

	////// end of debugging code ///////////


	////////////////////////////////  preprocessing the target set from {0,1} to {-1,0,1}  ///////////////




	tmp1=Cudd_BddToAdd(ddman,targetset);
	Cudd_Ref(tmp1);

	v = tmp1;

	trans_prev=Cudd_ReadOne(ddman);
	Cudd_Ref(trans_prev);

	v_prev=Cudd_ReadZero(ddman);
	Cudd_Ref(v_prev);

	////////////////////////////////////////////////////////////////////////////////////////
	u_check= new int[totbits];

	for (i=totbits-nbitsloop; i<totbits; i++)
		u_check[i] = 2;	

	//Build for permutation of x,u,x' to x',u,x
	permutation =  new int[totbits];
	
	// nbitsloop=bits(x)+bits(u)

	//initialize u
	for (i=totbits-nbitsloop; i<nbitsloop; i++)
		permutation[i] = i;

	//reorder initial x (u remains the same)
	for (i=0; i<totbits-nbitsloop; i++)
		permutation[i] = nbitsloop+i;

	//reorder final x
	for (i=nbitsloop; i<totbits; i++)
		permutation[i] = i-nbitsloop;	

	////////////////////////////////////////////////////////////////////////////////////////

	// Build cube for existential check of u and x'
	existential = new int[totbits];

	//initialize everything to 2
	for (i=0; i<totbits-nbitsloop; i++)
		existential[i] = 2;

	//final x and u set to 1 -- rest is ignored
	for (i=totbits-nbitsloop; i<totbits; i++)
		existential[i] = 1;	

	// another existential for x' only
	existential2 = new int[totbits];

	//initialize everything to 2 
	for (i=0; i<totbits-nbitsx; i++)
		existential2[i] = 2;

	//final x to 1 -- rest is ignored
	for (i=totbits-nbitsx; i<totbits; i++)
		existential2[i] = 1;

	tmp1 = Cudd_CubeArrayToBdd(ddman, existential);
	Cudd_Ref(tmp1);

	tcube_u =Cudd_BddToAdd(ddman,tmp1);
	Cudd_Ref(tcube_u);

	Cudd_RecursiveDeref(ddman,tmp1);

	tmp1 = Cudd_CubeArrayToBdd(ddman, existential2);
	Cudd_Ref(tmp1);

	tcube_xp =Cudd_BddToAdd(ddman,tmp1);
	Cudd_Ref(tcube_xp);

	Cudd_RecursiveDeref(ddman,tmp1);

	////////////////////////////////////////////////////////////////////////////////////////


	/////////  Cons is used to get the negation of T(x,u,x') //////////////////////////

	cons=Cudd_addConst(ddman,1);
	Cudd_Ref(cons);

	cons2=Cudd_addConst(ddman,0.5);
	Cudd_Ref(cons2);

	cons2m=Cudd_addConst(ddman,-0.5);
	Cudd_Ref(cons2m);


	tmp1=Cudd_addApply(ddman,Cudd_addTimes,trans,cons2);
	Cudd_Ref(tmp1);

	Cudd_RecursiveDeref(ddman,trans);

	trans=tmp1;


 	suc=nbitsx;

 	boundu = new int[totbits];
 	for(i=0;i<m;i++)
 	{

 		for(j=num[n+i]+1;j<(1<<nbits[n+i]);j++)
 		{

 			for (k=0; k<totbits; k++)
 				boundu[k]=2;

			// construct the out of bound u
			for(k=0;k<nbits[n+i];k++)
				boundu[suc+k]=(j>>(nbits[n+i]-k-1))& 1;

			// printf("for input %d boundu[%d]= ",i,j);
			// for(k=0;k<totbits;k++)
			// 	printf("%d",boundu[k]);
			// printf("\n");


			// convert to bdd using cubearraytobdd function
			tmp1 = Cudd_CubeArrayToBdd(ddman, boundu);
			Cudd_Ref(tmp1);
			// invert then convert to ADD	
			tmp2=Cudd_BddToAdd(ddman,tmp1);
			Cudd_Ref(tmp2);

			Cudd_RecursiveDeref(ddman,tmp1);

			tmp1=Cudd_addApply(ddman,Cudd_addTimes,cons2,tmp2);
			Cudd_Ref(tmp1);

			Cudd_RecursiveDeref(ddman,tmp2);

			// multiply with trans to get rid of this input

			tmp2=Cudd_addApply(ddman,Cudd_addPlus,tmp1,trans);
			Cudd_Ref(tmp2);

			Cudd_RecursiveDeref(ddman,tmp1);
			Cudd_RecursiveDeref(ddman,trans);

			trans=tmp2;

 		}
 		suc+=nbits[n+i];	
 	}

	delete boundu;

	count=0;






	do
	{
		//printf("Update started\n");
		count++;
		refine_count_current=refine_count;

	

		Cudd_RecursiveDeref(ddman, trans_prev);
		trans_prev=trans;

		Cudd_RecursiveDeref(ddman, v_prev);

		tmp1=Cudd_addBddInterval(ddman,v, 0.6,1);
		Cudd_Ref(tmp1);

		Cudd_RecursiveDeref(ddman,v);



		v_prev=Cudd_BddToAdd(ddman,tmp1);
		Cudd_Ref(v_prev);

		Cudd_RecursiveDeref(ddman,tmp1);

		//printf("start the loop\n");

		////////// debug code



		// printf("----------------------------------------- printing states inside v_prev at the begining of iteration %d !!!!!!!!!!! \n",count);
		// gen11=Cudd_FirstCube(ddman,v_prev,&testu,&value);
		// if(!Cudd_IsGenEmpty(gen11))
		//     {
		// for(jj=0;jj<totbits;jj++)
		//     printf("testu[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }


		// for(cubes=0;cubes<6;cubes++)
		// if(!(Cudd_IsGenEmpty(gen11)))
		// {
		// suc11=Cudd_NextCube(gen11,&testu,&value);
		// for(jj=0;jj<totbits;jj++)
		//     printf("testu[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		//}

		////// end of debugging code ///////////


		///////////////////  Update V ///////////////////


		//printf("just started update V\n");
		tmp1=Cudd_addPermute(ddman, v_prev, permutation);
		Cudd_Ref(tmp1);

		tmp2=Cudd_addApply(ddman,Cudd_addMinimum,trans_prev,tmp1);
		Cudd_Ref(tmp2);

		Cudd_RecursiveDeref(ddman,tmp1);


		tmp3=Cudd_addApply(ddman,Cudd_addMinus,cons,trans_prev);
		Cudd_Ref(tmp3);


		tmp4=Cudd_addApply(ddman,Cudd_addMaximum,tmp2,tmp3);
		Cudd_Ref(tmp4);

		Cudd_RecursiveDeref(ddman,tmp2);
		Cudd_RecursiveDeref(ddman,tmp3);

		tmp5=Cudd_addUnivAbstract(ddman,tmp4,tcube_xp);
		Cudd_Ref(tmp5);

		Cudd_RecursiveDeref(ddman,tmp4);

		tmp6=Cudd_addExistAbstract(ddman,tmp5,tcube_u);
		Cudd_Ref(tmp6);

		Cudd_RecursiveDeref(ddman,tmp5);

		v=Cudd_addApply(ddman,Cudd_addMaximum,v_prev,tmp6);
		Cudd_Ref(v);

		Cudd_RecursiveDeref(ddman,tmp6);

		//printf("Update done\n");


		/////////////////  End of Update V ///////////////////



		///////// debug code ////////////



		// printf("$$$$$$ printing states inside V_prev  at iteration %d $$$$ \n",count);
		// gen11=Cudd_FirstCube(ddman,v_prev,&testu,&value);
		// if(!Cudd_IsGenEmpty(gen11))
		//     {
		//     	if(value>0.2)
		//     	{

		// 	for(k=0;k<n;k++){
		// 	test_x[k]=min[k];
		// 	for(jj=0;jj<nbits[k];jj++)
		// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
		// 	test_x[k]=test_x[k]*eta;
		// 	printf("x[%d]=%f\n",k,test_x[k]);
		// 	}
		    
		// printf("value=%f\n",value);
		// }
		// for(jj=0;jj<nbitsx;jj++)
		// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }


		// for(cubes=0;cubes<30;cubes++)
		// if(!(Cudd_IsGenEmpty(gen11)))
		// {
		// suc11=Cudd_NextCube(gen11,&testu,&value);
		// if(value>0.2){
		// for(k=0;k<n;k++){
		// 	test_x[k]=min[k];
		// 	for(jj=0;jj<nbits[k];jj++)
		// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
		// 	test_x[k]=test_x[k]*eta;
		// 	printf("x[%d]=%f\n",k,test_x[k]);
		// 	}
		// printf("value=%f\n",value);
		// }
		// for(jj=0;jj<nbitsx;jj++)
		// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }









		// printf("$$$$$$ printing states inside V  at iteration %d $$$$ \n",count);
		// gen11=Cudd_FirstCube(ddman,v,&testu,&value);
		// if(!Cudd_IsGenEmpty(gen11))
		//     {
		//     	if(value>0.2)
		//     	{

		// 	for(k=0;k<n;k++){
		// 	test_x[k]=min[k];
		// 	for(jj=0;jj<nbits[k];jj++)
		// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
		// 	test_x[k]=test_x[k]*eta;
		// 	printf("x[%d]=%f\n",k,test_x[k]);
		// 	}
		    
		// printf("value=%f\n",value);
		// }
		// for(jj=0;jj<nbitsx;jj++)
		// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }


		// for(cubes=0;cubes<30;cubes++)
		// if(!(Cudd_IsGenEmpty(gen11)))
		// {
		// suc11=Cudd_NextCube(gen11,&testu,&value);
		// if(value>0.2){
		// for(k=0;k<n;k++){
		// 	test_x[k]=min[k];
		// 	for(jj=0;jj<nbits[k];jj++)
		// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
		// 	test_x[k]=test_x[k]*eta;
		// 	printf("x[%d]=%f\n",k,test_x[k]);
		// 	}
		// printf("value=%f\n",value);
		// }
		// for(jj=0;jj<nbitsx;jj++)
		// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }




			////// end of debugging code ///////////




		if(Cudd_addLeq(ddman,v,v_prev) && Cudd_addLeq(ddman,v_prev,v))
		{
			printf("now v and v_prev are the same\n");
			fixed=1;
			break;
		}


		// printf("just finished update v in iteration number %d\n",count);


		// 	///// debugging code to be removed /////
		// if(count==1)
		// {
		// fp1=fopen("v1_test.dot", "w");
		// 			Cudd_DumpDot(ddman, numroots_design, &v,NULL,NULL,fp1); 
		// 			fclose(fp1);
		// Dddmp_cuddAddStore(ddman, NULL, v, NULL,
		// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
		// 					"v1_test.bdd", NULL);
		// }
		// if(count==2)
		// {
		// fp1=fopen("v2_test.dot", "w");
		// 			Cudd_DumpDot(ddman, numroots_design, &v,NULL,NULL,fp1); 
		// 			fclose(fp1);
		// Dddmp_cuddAddStore(ddman, NULL, v, NULL,
		// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
		// 					"v2_test.bdd", NULL);
		// }

		// printf("####################################### printing states inside V  before changing from -1 0 1 for iteration number %d \n",count);
		// gen11=Cudd_FirstCube(ddman,v,&testu,&value);
		// if(!Cudd_IsGenEmpty(gen11))
		//     {
		// for(jj=0;jj<totbits;jj++)
		//     printf("testu[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }


		// for(cubes=0;cubes<6;cubes++)
		// if(!(Cudd_IsGenEmpty(gen11)))
		// {
		// suc11=Cudd_NextCube(gen11,&testu,&value);
		// for(jj=0;jj<totbits;jj++)
		//     printf("testu[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }

			////// end of debugging code ///////////

		/// old code due to using -1 0 1

		tmp1=Cudd_addApply(ddman,Cudd_addMinus,v,cons);
		Cudd_Ref(tmp1);

		tmp2=Cudd_addApply(ddman,Cudd_addTimes,tmp1,cons2);
		Cudd_Ref(tmp2);

		Cudd_RecursiveDeref(ddman,tmp1);

		// end of old code

		//printf("about to get the set of states need to be discovered\n");


		////////// debug code



		// printf("####  printing states 1 1/2 0 inside V  at iteration %d ######## \n",count);
		// gen11=Cudd_FirstCube(ddman,v,&testu,&value);
		// if(!Cudd_IsGenEmpty(gen11))
		//     {
		//     	if(value>0.2)
		//     	{

		// 	for(k=0;k<n;k++){
		// 	test_x[k]=min[k];
		// 	for(jj=0;jj<nbits[k];jj++)
		// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
		// 	test_x[k]=test_x[k]*eta;
		// 	printf("x[%d]=%f\n",k,test_x[k]);
		// 	}
		    
		// printf("value=%f\n",value);
		// }
		// for(jj=0;jj<nbitsx;jj++)
		// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }


		// for(cubes=0;cubes<30;cubes++)
		// if(!(Cudd_IsGenEmpty(gen11)))
		// {
		// suc11=Cudd_NextCube(gen11,&testu,&value);
		// if(value>0.2){
		// for(k=0;k<n;k++){
		// 	test_x[k]=min[k];
		// 	for(jj=0;jj<nbits[k];jj++)
		// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
		// 	test_x[k]=test_x[k]*eta;
		// 	printf("x[%d]=%f\n",k,test_x[k]);
		// 	}
		// printf("value=%f\n",value);
		// }
		// for(jj=0;jj<nbitsx;jj++)
		// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
		// printf("value=%f\n",value);
		// }





		////// end of debugging code ///////////

		int num_cubes=0;
		vector<vector<int> > x_ref_new;
		int tmp_num_dc;
		int test_endof_cubes;
		num_dc=0;

		gen=Cudd_FirstCube(ddman,v,&x_ref,&value);
		if(Cudd_IsGenEmpty(gen))
		{
			//	printf("breaking because no more states to verify is empty 1st break\n");
			refine_count_current=0;
		}
		while(1)    //value != 0.5
		{
			//	printf("######getting states that needs to be explored, so far we have %d states after discovering %d cubes####\n",num_dc,num_cubes);
			if(num_dc>=refine_count)
			{
			//printf("breaking because we got the states needed to be discovered\n");
			break;
			}
			
			if(value != 0.5)
			{
				if(Cudd_IsGenEmpty(gen))
				{
					//printf("breaking because no more states to verify is empty 1st break\n");
		    		refine_count_current=0;
					break;
				}
				suc=Cudd_NextCube(gen,&x_ref,&value);
				//printf("we have the next cube now!!!!\n");
			}
			if(value==0.5)
			{
				test_endof_cubes=1;
				for(i=0;i<nbitsx;i++)
					if(x_ref[i]!=2)
						test_endof_cubes=0;
				if(test_endof_cubes==0)
				{
					// 	printf("currentcube");
					// 		for(i=0;i<nbitsx;i++)
					// 			{
					// printf("%d",x_ref[i]);//x_ref_new[i][j];
					// 			}
					// 		printf("\n");
					tmp_num_dc=1;
		       		x_ref_new.push_back(vector<int>());  
		        	for (i=0;i<nbitsx;i++)
					{	
		     			x_ref_new[num_cubes].push_back(x_ref[i]);
		     			if(x_ref[i]==2)
		         			tmp_num_dc=tmp_num_dc<<1;
		   			}
		   			num_dc+= tmp_num_dc;
		   			num_cubes++;
		   		}
		   	}
		  value=0;

		}
       	
       
   
		if(num_cubes==0)
	    {
	   			// 	printf("breaking because no more states to verify is empty 2nd break\n");
	    		refine_count_current=0;
	   	}




     
			//printf("--------------------- number of states =1/2  is  %d at iteration %d ---------------------\n",num_dc,count);


		// for(i=0;i<x_ref_new.size();i++)
		// 	{
		// 		printf("x_ref_new[%d]=",i);
		// 	for(j=0;j<nbitsx;j++)
		// 	{
		// 			printf("%d",x_ref_new[i][j]);
		// 	}
		// 	printf("\n");
		// 	}

		i=0;
     	refine_count_current=num_dc;

     	if(num_dc>refine_count)
     	{
     		//num_dc=refine_count;
     		refine_count_current=refine_count;
		 }

		//	printf("num_dc=%d\n",num_dc);
   		while((num_dc>1)  &&(i<x_ref_new.size())&& (num_cubes !=0))//
   		{
        	//printf("after expanding state %d  and now num_dc=%d\n",i,num_dc);
     		for (j=0;j<nbitsx;j++)
   			{
     			if(x_ref_new[i][j]==2)
     				{
			         for(k=0;k<nbitsx;k++)
			             x_ref[k]=x_ref_new[i][k];
			         x_ref_new.erase(x_ref_new.begin()+i);
			         x_ref_new.push_back(vector<int>());
			         x_ref_new.push_back(vector<int>());
			         i--;
			         break;
     				}
         
   			}
     		for(k=0;k<nbitsx;k++)
     		{
         		if(k==j)
         		{
           			x_ref_new[x_ref_new.size()-2].push_back(0);
           			x_ref_new[x_ref_new.size()-1].push_back(1);
         		}
         		else
         		{
            		x_ref_new[x_ref_new.size()-2].push_back(x_ref[k]);
            		x_ref_new[x_ref_new.size()-1].push_back(x_ref[k]);
         		}
     		}
    
			// for(k=0;k<x_ref_new.size();k++)
			// 	{
			// 		printf("x_ref_new[%d]=",k);
			// 	for(j=0;j<nbitsx;j++)
			// 	{
			// 			printf("%d",x_ref_new[k][j]);
			// 	}
			// 	printf("\n");
			// 	}
	     
	   		//  num_dc--;
     		i++;
    	}	

		for(i=0;i<x_ref_new.size();i++)
			for(j=0;j<nbitsx;j++)
				if(x_ref_new[i][j]==2)
					x_ref_new[i][j]=0;

		Cudd_RecursiveDeref(ddman,tmp2);

		///// debug code




		// printf("after expanding the states \n");
		// for(i=0;i<x_ref_new.size();i++)
		// 	{
		// 		printf("x_ref_new[%d]=",i);
		// 	for(j=0;j<nbitsx;j++)
		// 	{
		// 			printf("%d",x_ref_new[i][j]);
		// 	}
		// 	printf("\n");
		// 	}



		////// end debug ///////

		/// old code due to using -1 0 1 in Trans

		// change T(x,u,x') values from {1,0,-1} to {1,1/2,0}

		//tmp1=Cudd_addApply(ddman,Cudd_addMinus,trans_prev,cons);
		//Cudd_Ref(tmp1);

		//Cudd_RecursiveDeref(ddman,trans_prev);

		//trans_prev=Cudd_addApply(ddman,Cudd_addTimes,tmp1,cons2);
		//Cudd_Ref(trans_prev);

		//Cudd_RecursiveDeref(ddman,tmp1);

		//// end of old code

		//printf("just started refine T\n");



		//// debug ////////


		// if(count==1)
		// {
		// fp1=fopen("trans_prev1_test.dot", "w");
		// 			Cudd_DumpDot(ddman, numroots_design, &trans_prev,NULL,NULL,fp1); 
		// 			fclose(fp1);
		// Dddmp_cuddAddStore(ddman, NULL, trans_prev, NULL,
		// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
		// 					"trans_prev1_test.bdd", NULL);
		// }
		// if(count==2)
		// {
		// fp1=fopen("trans_prev2_test.dot", "w");
		// 			Cudd_DumpDot(ddman, numroots_design, &trans_prev,NULL,NULL,fp1); 
		// 			fclose(fp1);
		// Dddmp_cuddAddStore(ddman, NULL, trans_prev, NULL,
		// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
		// 					"trans_prev2_test.bdd", NULL);
		// }

		////// end debug ///////



		///////////////////  Refine T ///////////////////
		//printf("Refine started\n");

		for(i=0;i<refine_count_current;i++)
		{

			//printf("Refine started, this is state number %d from %d states\n",i+1,refine_count_current);

			if(i>0)
			trans_prev=trans;
			break_refine=0;
			for (j=0; j<totbits-nbitsloop; j++)
				u_check[j]=x_ref_new[i][j];



			tmp1 = Cudd_CubeArrayToBdd(ddman, u_check);
			Cudd_Ref(tmp1);

			tmp2 =Cudd_BddToAdd(ddman,tmp1);
			Cudd_Ref(tmp2);

			Cudd_RecursiveDeref(ddman,tmp1);

			tmp1=Cudd_addApply(ddman,Cudd_addTimes,trans_prev,tmp2);
			Cudd_Ref(tmp1);

			Cudd_RecursiveDeref(ddman,tmp2);


			gen=Cudd_FirstCube(ddman,tmp1,&u_ref,&value);
			if(Cudd_IsGenEmpty(gen))
		    {
		    	break_refine=1;
		    }
			else
				while(value != 0.5 && break_refine !=1)
				{
					if(Cudd_IsGenEmpty(gen))
			    	{
			    		break_refine=1;
			    	}
			    	else
						suc=Cudd_NextCube(gen,&u_ref,&value);

				}
			Cudd_RecursiveDeref(ddman,tmp1);

			// for(k=0;k<totbits;k++)
			// 	printf("u_ref[%d]=%d \n",k,u_ref[k]);

			// if(count>0)
			// {	printf("break refine = %d\n",break_refine);
			// for(k=totbits-nbitsloop;k<totbits-nbitsx;k++)
			// 	printf("picked u to refine is u[%d]=%d \n",k,u_ref[k]);
			// }
			//call function to compute post(x_ref,u_ref)
			// note u_ref might have dont cares   -----fixed 
			for (k=totbits-nbitsloop; k<totbits-nbitsx; k++)
				if(u_ref[k] == 2)
					u_ref[k]=0;
			j=0;
			for(k=0;k<m;k++)
			{
				u_oob[k]=min[k+n];
				for(jj=0;jj<nbits[k+n];jj++)
				u_oob[k]+=u_ref[j+nbitsx+jj]<<(nbits[k+n]-jj-1);
				j+=nbits[k+n];
				u_oob[k]=u_oob[k]*mu;
			}

			for(k=0;k<m;k++)
			{
				//	printf("u_oob[%d]=%f \n",k,u_oob[k]);
				if(u_oob[k]>max[k]*mu)
				{
					
					//	printf("u picked outside the range of u\n");
					break_refine=1;
				}				
			}

			///// commented for now, instead of breaking don't call computereach **************
			//if(break_refine==1)
			//	break;



			// nullifying all transtions from this state for now.
			tmp1 = Cudd_CubeArrayToBdd(ddman, u_ref);
			Cudd_Ref(tmp1);

			// if(count>8)
			// {	printf("break refine = %d\n",break_refine);
			// for(k=0;k<totbits;k++)
			// 	printf("picked u to nullify is u[%d]=%d \n",k,u_ref[k]);
			// }

			tmp2 =Cudd_BddToAdd(ddman,tmp1);
			Cudd_Ref(tmp2);

			Cudd_RecursiveDeref(ddman,tmp1);

			tmp1=Cudd_addApply(ddman,Cudd_addTimes,tmp2,cons2m);
			Cudd_Ref(tmp1);
			Cudd_RecursiveDeref(ddman,tmp2);

			if(break_refine==0)
			{
				tmp2=Cudd_addApply(ddman,Cudd_addPlus,trans_prev,tmp1);
				Cudd_Ref(tmp2);
			}
			else
			{
				tmp2=Cudd_addApply(ddman,Cudd_addMinus,trans_prev,tmp1);
				Cudd_Ref(tmp2);
			}

			Cudd_RecursiveDeref(ddman,tmp1);
			//Cudd_RecursiveDeref(ddman,trans_prev);


			//printf("calling computereach for in iteration number %d \n",count);
			for(k=0;k<totbits-nbitsloop;k++)
			x_ref[k]=x_ref_new[i][k];

			if(break_refine==0)
				tmp1=ComputeReach(ddman,x_ref,u_ref);
			else
			{
				tmp1=Cudd_ReadZero(ddman);
				Cudd_Ref(tmp1);
			}
			//Cudd_Ref(tmp1);   already ref inside the function....


			trans=Cudd_addApply(ddman,Cudd_addPlus,tmp2,tmp1);
			Cudd_Ref(trans);

			Cudd_RecursiveDeref(ddman,tmp1);
			Cudd_RecursiveDeref(ddman,tmp2);



		}
		if(refine_count_current==0)
			trans=trans_prev;

		//printf("Refine done\n");
		///////////////////  End of Refine T ///////////////////	




		//// debugging code

		// if(count==1)
		// {
		// fp1=fopen("t1_test.dot", "w");
		// 			Cudd_DumpDot(ddman, numroots_design, &trans,NULL,NULL,fp1); 
		// 			fclose(fp1);
		// Dddmp_cuddAddStore(ddman, NULL, trans, NULL,
		// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
		// 					"t1_test.bdd", NULL);
		// }
		// if(count==2)
		// {
		// fp1=fopen("t2_test.dot", "w");
		// 			Cudd_DumpDot(ddman, numroots_design, &trans,NULL,NULL,fp1); 
		// 			fclose(fp1);
		// Dddmp_cuddAddStore(ddman, NULL, trans, NULL,
		// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
		// 					"t2_test.bdd", NULL);
		// }


		/////end of debugging code


		// change T(x,u,x') values from  {1,1/2,0}  to {1,0,-1} 

		/// old code due to using -1 0 1

		//tmp1=Cudd_addApply(ddman,Cudd_addDivide,trans,cons2);
		//Cudd_Ref(tmp1);

		//Cudd_RecursiveDeref(ddman,trans);

		//trans=Cudd_addApply(ddman,Cudd_addPlus,tmp1,cons);
		//Cudd_Ref(trans);

		//Cudd_RecursiveDeref(ddman,tmp1);


		/// end of old code

		//////// test iteration

		//	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~We are in iteration %d and comparison leads to  %d %d~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",count,Cudd_addLeq(ddman,v,v_prev) ,Cudd_addLeq(ddman,v_prev,v));
		/////debug code ///////////
		if(count>0)
		{	




			///////// debug code ////////////


			// //printf("just saved v\n");

			// printf("$$$$$$ printing trans inside trans_prev  at iteration %d $$$$ \n",count);
			// gen11=Cudd_FirstCube(ddman,trans_prev,&testu,&value);
			// if(!Cudd_IsGenEmpty(gen11))
			//     {
			//     	if(value>0.2)
			//     	{
			// 	for(jj=0;jj<totbits;jj++)
			//     printf("testu[%d]=%d\n",jj,testu[jj]);
			// printf("value=%f\n",value);
			// }
			// }


			// for(cubes=0;cubes<30;cubes++)
			// if(!(Cudd_IsGenEmpty(gen11)))
			// {
			// suc11=Cudd_NextCube(gen11,&testu,&value);
			// if(value>0.2){
			// 	for(jj=0;jj<totbits;jj++)
			//     printf("testu[%d]=%d\n",jj,testu[jj]);
			// printf("value=%f\n",value);
			// }
			// }

			// 	////// end of debugging code ///////////



			// tmp1=Cudd_addApply(ddman,Cudd_addMinus,trans,cons);
			// Cudd_Ref(tmp1);
			// printf("$$$$$$ printing trans inside trans  at iteration %d $$$$ \n",count);
			// gen11=Cudd_FirstCube(ddman,tmp1,&testu,&value);
			// if(!Cudd_IsGenEmpty(gen11))
			//     {
			//     	if(value>0.5)
			//     	{
			// 	for(jj=0;jj<totbits;jj++)
			//     printf("testu[%d]=%d\n",jj,testu[jj]);
			// printf("value=%f\n",value);
			// }
			// }


			// for(cubes=0;cubes<30;cubes++)
			// if(!(Cudd_IsGenEmpty(gen11)))
			// {
			// suc11=Cudd_NextCube(gen11,&testu,&value);
			// if(value>0.5){


			// 	for(jj=0;jj<totbits;jj++)
			//     printf("testu[%d]=%d\n",jj,testu[jj]);
			// printf("value=%f\n",value);


			// // for(k=0;k<n;k++){
			// // 	test_x[k]=min[k];
			// // 	for(jj=0;jj<nbits[k];jj++)
			// // 		test_x[k]+=testu[k*nbits[0]+jj]<<jj;
			// // 	test_x[k]=test_x[k]*eta;
			// // 	printf("x[%d]=%f\n",k,test_x[k]);
			// // 	}
			// // 	for(k=0;k<m;k++){
			// // 	u_oob[k]=min[k+n];
			// // 	for(jj=0;jj<nbits[k+n];jj++)
			// // 		u_oob[k]+=testu[nbitsx+jj]<<jj;
			// // 	u_oob[k]=u_oob[k]*mu;
			// // 	printf("u[%d]=%f\n",k,u_oob[k]);
			// // }
			// // printf("value=%f\n",value);
			// }
			// }
			// Cudd_RecursiveDeref(ddman,tmp1);



				////// end of debugging code ///////////



		}





		////// end of debug code //////////





		if(Cudd_addLeq(ddman,v,v_prev) && Cudd_addLeq(ddman,v_prev,v))// && Cudd_addLeq(ddman,trans,trans_prev) && Cudd_addLeq(ddman,trans_prev,trans))
		{	
			fixed=1;
			printf("number of iterations = %d so far *****\n",count);
			//printf("number of iterations = %d\n",count);
		}		
	}while(fixed != 1);
	printf("number of iterations = %d and fixed = %d\n",count,fixed);

	// 	fp1=fopen("v_test.dot", "w");
	// 			Cudd_DumpDot(ddman, numroots_design, &v,NULL,NULL,fp1); 
	// 			fclose(fp1);
	// Dddmp_cuddAddStore(ddman, NULL, v, NULL,
	// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
	// 					"v_test.bdd", NULL);

	////////// debug code



	// printf("####################################### printing states inside V  finally !!!!!!!!!!! \n");



	// gen11=Cudd_FirstCube(ddman,v,&testu,&value);
	// if(!Cudd_IsGenEmpty(gen11))
	//     {
	//     	if(value>0.5)
	//     	{

	// 	for(k=0;k<n;k++){
	// 	test_x[k]=min[k];
	// 	for(jj=0;jj<nbits[k];jj++)
	// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
	// 	test_x[k]=test_x[k]*eta;
	// 	printf("x[%d]=%f\n",k,test_x[k]);
	// 	}
	    
	// printf("value=%f\n",value);
	// }
	// for(jj=0;jj<nbitsx;jj++)
	// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
	// printf("value=%f\n",value);
	// }


	// for(cubes=0;cubes<200;cubes++)
	// if(!(Cudd_IsGenEmpty(gen11)))
	// {
	// suc11=Cudd_NextCube(gen11,&testu,&value);
	// if(value>0.5){
	// for(k=0;k<n;k++){
	// 	test_x[k]=min[k];
	// 	for(jj=0;jj<nbits[k];jj++)
	// 		test_x[k]+=testu[k*nbits[0]+jj]<<(nbits[k]-jj-1);
	// 	test_x[k]=test_x[k]*eta;
	// 	printf("x[%d]=%f\n",k,test_x[k]);
	// 	}
	// printf("value=%f\n",value);
	// }
	// for(jj=0;jj<nbitsx;jj++)
	// 		printf("x_ref[%d]=%d\n",jj,testu[jj]);
	// printf("value=%f\n",value);
	// }



	////// end of debugging code ///////////



	//printf("just finished safety1\n");

	delete u_oob;
	delete permutation;
	delete existential;
	delete existential2;
	delete u_check;

	delete test_x;

	Cudd_RecursiveDeref(ddman,cons);
	Cudd_RecursiveDeref(ddman,cons2);
	Cudd_RecursiveDeref(ddman,cons2m);
	Cudd_RecursiveDeref(ddman,v_prev);
	Cudd_RecursiveDeref(ddman,trans_prev);
	printf("just finished safety1\n");


}



int System::SaveController(DdManager* ddman, string filename) {

DdNode *tmp1,*tmp2;
FILE *fp1;
int numroots_design=1; 
tmp1=Cudd_addBddInterval(ddman,v, 1,1);
Cudd_Ref(tmp1);

fp1=fopen("cont_dom.dot", "w");
			Cudd_DumpDot(ddman, numroots_design, &tmp1,NULL,NULL,fp1); 
			fclose(fp1);
Dddmp_cuddAddStore(ddman, NULL, tmp1, NULL,
					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
					"cont_dom.bdd", NULL);

printf("saved the domain\n");
Cudd_RecursiveDeref(ddman,tmp1);



// tmp2=Cudd_addBddInterval(ddman,trans, 1,1);
// Cudd_Ref(tmp2);

// fp1=fopen("cont.dot", "w");
// 			Cudd_DumpDot(ddman, numroots_design, &tmp2,NULL,NULL,fp1); 
// 			fclose(fp1);
// Dddmp_cuddAddStore(ddman, NULL, tmp2, NULL,
// 					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
// 					"cont.bdd", NULL);

// printf("saved the controller\n");
// Cudd_RecursiveDeref(ddman,tmp2);


printf("saving complete now deleting the trans and v ADDs\n");

Cudd_RecursiveDeref(ddman,trans);
	Cudd_RecursiveDeref(ddman,v);

return 0;
}

int System::SaveParams(string filename){

int i;
ofstream Output(filename + ".txt");
	if (!Output.is_open()) {
		cerr << "Error opening file:"<<filename + ".txt\n";
		return 0;
	}
Output<<n<<"\t"<<m<<"\n";
Output<<tau<<"\t"<<eta<<"\t"<<mu<<"\t"<<epsilon<<"\n";
for(i=0;i<n+m;i++)
Output<<min[i]<<"\t";
Output<<"\n";
for(i=0;i<n+m;i++)
Output<<max[i]<<"\t";
Output<<"\n";

for(i=0;i<n+m;i++)
Output<<num[i]<<"\t";
Output<<"\n";

for(i=0;i<2*n+m;i++)
Output<<nume[i]<<"\t";
Output<<"\n";

for(i=0;i<n;i++)
Output<<nxbits[i]<<"\t";
Output<<"\n";
for(i=0;i<m;i++)
Output<<nubits[i]<<"\t";
Output<<"\n";

Output<<nbitsx<<"\n";
Output<<nbitsloop<<"\n";

for(i=0;i<2*n+m;i++)
Output<<nbits[i]<<"\t";
Output<<"\n";

Output<<totbits<<"\n";

Output<<type<<"\n";

Output.close();

//printf("saving controller done successfuly!! \n");
return 1;


}

DdNode* System::ComputeReach(DdManager* ddman,int* x_ref,int* u_ref){


DdNode *added_trans, *tmp1, *tmp2;

added_trans=Cudd_ReadZero(ddman);
	Cudd_Ref(added_trans);


int i,j,k,ntrans;
int jj,kk;
double temp;
// to deal with the mod function problem
double prec=0.001;

double etah=0.5*eta;
int* xend_temp=new int[nbitsx];

int* current_trans=new int[totbits];

double* x=new double[n];

double* u=new double[m];
// double check this part, and make sure it is not the reversed order of bits
jj=0;
for(i=0;i<n;i++){
	x[i]=min[i];
	for(j=0;j<nbits[i];j++)
		x[i]+=x_ref[jj+j]<<(nbits[i]-j-1);
	jj+=nbits[i];
	x[i]=x[i]*eta;
}
jj=nbitsx;
for(i=0;i<m;i++){
	u[i]=min[i+n];
	for(j=0;j<nbits[i+n];j++)
		u[i]+=u_ref[jj+j]<<(nbits[i+n]-j-1);
	jj+=nbits[i+n];
	u[i]=u[i]*mu;
}

jj=0;



// for(i=0;i<n;i++)
// 	printf("x[%d]=%f\n",i,x[i]);
// for(i=0;i<m;i++)
// 	printf("u[%d]=%f\n",i,u[i]);

double* xend=new double[n];

double* max_val=new double[n];

double* min_val=new double[n];

int* max_label=new int[n];

int* min_label=new int[n];

int* min_label_temp=new int[n];

double* eta_new=new double[n];

///////////// simulate the point /////////


// double integrator



// xend[0]=x[0]+x[1]*tau+u[0]*tau*tau;
// xend[1]=x[1]+u[0]*tau;

// eta_new[0]=0.5*eta*(1+tau);
// eta_new[1]=0.5*eta;

// for(i=0;i<n;i++){
// 	max_val[i]=xend[i]+eta_new[i];
// 	min_val[i]=xend[i]-eta_new[i];
// }



///// end of double integrator


// unicycle



double wt=u[1]*tau;
double wth=0.5*wt;
double kvw;
if(u[1] !=0)
	kvw=2*(u[0]/u[1])*sin(wth);
else
	kvw=u[0]*tau;

double x1plus1=x[0]+kvw*cos(x[2]+wth+etah);
double x1plus2=x[0]+kvw*cos(x[2]+wth-etah);


double x2plus1=x[1]+kvw*sin(x[2]+wth+etah);
double x2plus2=x[1]+kvw*sin(x[2]+wth-etah);

double x1plus3,x2plus3;
double pi=3.14159265;

if(abs(x[2]+wth)<etah)
	x1plus3=x[0]+kvw;
else if ((abs(x[2]+wth-pi)<etah) ||(abs(x[2]+wth+pi)<etah))
	x1plus3=x[0]-kvw;
else
	x1plus3=x1plus2;


if(abs(x[2]+wth-pi/2)<etah)
	x2plus3=x[1]+kvw;
else if ((abs(x[2]+wth+pi/2)<etah))
	x2plus3=x[1]-kvw;
else
	x2plus3=x2plus2;

eta_new[0]=etah;
eta_new[1]=etah;
eta_new[2]=etah;

double temp_min,temp_max;

temp_min=x1plus1;
temp_max=x1plus1;

if(x1plus2>temp_max)
temp_max=x1plus2;
if(x1plus2<temp_min)
temp_min=x1plus2;

if(x1plus3>temp_max)
temp_max=x1plus3;
if(x1plus3<temp_min)
temp_min=x1plus3;


max_val[0]=temp_max+etah;
min_val[0]=temp_min-etah;

temp_min=x2plus1;
temp_max=x2plus1;

if(x2plus2>temp_max)
temp_max=x2plus2;
if(x2plus2<temp_min)
temp_min=x2plus2;

if(x2plus3>temp_max)
temp_max=x2plus3;
if(x2plus3<temp_min)
temp_min=x2plus3;


max_val[1]=temp_max+etah;
min_val[1]=temp_min-etah;

max_val[2]=x[2]+wt+etah;
min_val[2]=x[2]+wt-etah;


/// end of unicycle

//for(i=0;i<n;i++)
//	printf("xend[%d]=%f\n",i,xend[i]);
/////////////// end of simulate //////////



// for(i=0;i<n;i++)
// 	printf("min[%d]=%f\n",i,min_val[i]);

// for(i=0;i<n;i++)
// 	printf("max[%d]=%f\n",i,max_val[i]);


for(i=0;i<n;i++){

	temp=fmod(max_val[i],eta);
	//printf("temp=%f when i=%d and max[%d]=%f\n",temp,i,i,max_val[i]);
	if((temp > (etah+prec)) || (temp < (etah-prec)))
		max_label[i]=round(max_val[i]/eta)-min[i];
	else
		max_label[i]=floor(max_val[i]/eta)-min[i];
	
	temp=fmod(min_val[i],eta);
	if(temp<0)
		temp=temp*-1;
	//printf("temp=%f when i=%d and min[%d]=%f\n",temp,i,i,min_val[i]);
	if((temp > (etah+prec)) || (temp < (etah-prec)))
		min_label[i]=round(min_val[i]/eta)-min[i];
	else
		min_label[i]=ceil(min_val[i]/eta)-min[i];
}

// for(i=0;i<n;i++)
// 	printf("min_label[%d]=%d\n",i,min_label[i]);

// for(i=0;i<n;i++)
// 	printf("max_label[%d]=%d\n",i,max_label[i]);


// if out of domain for now set it to the boundary
for(i=0;i<n;i++)
	if(min_label[i]<0)
		min_label[i]=0;

for(i=0;i<n;i++)
	if(max_label[i]>num[i])
		max_label[i]=num[i];


for(i=0;i<nbitsx;i++)
current_trans[i]=x_ref[i];
for(i=nbitsx;i<totbits-nbitsx;i++)
current_trans[i]=u_ref[i];


for(i=0;i<n;i++)
	min_label_temp[i]=min_label[i];


// for(i=0;i<n;i++)
// 	printf("min_label[%d]=%d \n",i,min_label[i]);

// for(i=0;i<n;i++)
// 	printf("max_label[%d]=%d \n",i,max_label[i]);


k=1;
ntrans=0;

//printf("start adding labels\n");

while(k<n && min_label_temp[n-1]<=max_label[n-1]){


k=0;
kk=0;
for(i=0;i<n;i++)
{
	for(j=0;j<nbits[i];j++)
	{
		jj=(min_label_temp[i]>>(nbits[i]-j-1))& 1;
		xend_temp[kk+j]=jj;
		//printf("------------ xend_temp[%d]=%d, jj =%d for shift %d-------------\n",j,xend_temp[j],jj,j);
	}
	kk+=nbits[i];
}
	// for(i=nbitsloop;i<totbits;i++)
	// printf("xend_temp[%d]=%d \n",i-nbitsloop,xend_temp[i-nbitsloop]);

for(i=nbitsloop;i<totbits;i++)
current_trans[i]=xend_temp[i-nbitsloop];

//for(i=0;i<n;i++)
//	printf("min_label[%d]=%d \n",i,min_label[i]);

//for(i=0;i<n;i++)
//	printf("max_label[%d]=%d \n",i,max_label[i]);
//printf("k=%d n=%d min_label_temp[n]=%d  max_label[n]=%d \n",k,n,min_label_temp[n-1],max_label[n-1]);


//for(i=0;i<n;i++)
//	printf("min_label_temp[%d]=%d \n",i,min_label_temp[i]);
// for(i=0;i<totbits;i++)
// 	printf("current_trans[%d]=%d \n",i,current_trans[i]);
//printf("creating a cube for the transition\n");
	tmp1 = Cudd_CubeArrayToBdd(ddman,current_trans);
	Cudd_Ref(tmp1);

	tmp2 =Cudd_BddToAdd(ddman,tmp1);
	Cudd_Ref(tmp2);

	Cudd_RecursiveDeref(ddman,tmp1);
//printf("oring the two bdds\n");
	tmp1=Cudd_addApply(ddman,Cudd_addPlus,added_trans,tmp2);
	Cudd_Ref(tmp1);

	Cudd_RecursiveDeref(ddman,tmp2);
	Cudd_RecursiveDeref(ddman,added_trans);
//printf("adding a label to the bdd\n");
	added_trans=tmp1;

	while(k<(n-1) && min_label_temp[k]>=max_label[k]){
		min_label_temp[k]=min_label[k];
		k=k+1;
	}
	min_label_temp[k]+=1;
	ntrans++;
	//for(i=0;i<n;i++)
	//printf("min_label_temp[%d]=%d \n",i,min_label_temp[i]);

}
//printf("number of transitions added %d\n",ntrans);
if(ntrans==0)
{

// 	printf("number of transitions added %d\n",ntrans);
// printf("k=%d n=%d min_label_temp[n]=%d  max_label[n]=%d \n",k,n,min_label_temp[n-1],max_label[n-1]);
// for(i=0;i<n;i++)
// 	printf("x[%d]=%f\n",i,x[i]);
// for(i=0;i<n;i++){
// 	printf("x_ref[%d]=",i);
// 	for(j=0;j<nbits[i];j++)
// 		printf("%d",x_ref[i*nbits[0]+j]);
// 	printf("\n");
// }

// for(i=0;i<m;i++)
// 	printf("u[%d]=%f\n",i,u[i]);

// for(i=0;i<n;i++)
// 	printf("min_label[%d]=%d \n",i,min_label[i]);

// for(i=0;i<n;i++)
// 	printf("max_label[%d]=%d \n",i,max_label[i]);

// for(i=0;i<n;i++)
// 	printf("min[%d]=%f\n",i,min_val[i]);

// for(i=0;i<n;i++)
// 	printf("max[%d]=%f\n",i,max_val[i]);



}

// printf("number of transitions added %d\n",ntrans);
// printf("k=%d n=%d min_label_temp[n]=%d  max_label[n]=%d \n",k,n,min_label_temp[n-1],max_label[n-1]);
// for(i=0;i<n;i++)
// 	printf("min_global[%d]=%d\n",i,min[i]);
// 	for(i=0;i<n;i++)
// 	printf("max_global[%d]=%d\n",i,max[i]);
// for(i=0;i<n;i++)
// 	printf("x[%d]=%f\n",i,x[i]);
// for(i=0;i<n;i++){
// 	printf("x_ref[%d]=",i);
// 	for(j=0;j<nbits[i];j++)
// 		printf("%d",x_ref[i*nbits[0]+j]);
// 	printf("\n");
// 	}
// 	k=nbitsx;
// for(i=0;i<m;i++){
// 	printf("u_ref[%d]=",i);

// 	for(j=0;j<nbits[i+n];j++)
// 		printf("%d",u_ref[k+j]);
// 		// for(j=0;j<totbits;j++)
// 		// printf("%d",u_ref[j]);
// 	printf("\n");
// 	k+=nbits[i+n];
// 	}



// for(i=0;i<m;i++)
// 	printf("u[%d]=%f\n",i,u[i]);

// for(i=0;i<n;i++)
// 	printf("min_label[%d]=%d \n",i,min_label[i]);

// for(i=0;i<n;i++)
// 	printf("max_label[%d]=%d \n",i,max_label[i]);

// for(i=0;i<n;i++)
// 	printf("min[%d]=%f\n",i,min_val[i]);

// for(i=0;i<n;i++)
// 	printf("max[%d]=%f\n",i,max_val[i]);

// for(i=0;i<n;i++)
// 	printf("mod function max[%d]=%f\n",i,fmod(max_val[i],eta));

// for(i=0;i<n;i++)
// 	printf("mod function min[%d]=%f\n",i,fmod(min_val[i],eta));

// printf("number of transitions added %d\n",ntrans);


delete xend;
delete max_val;
delete min_val;
delete max_label;
delete min_label;
delete eta_new;
delete x;
delete u;
return added_trans;

}

int System::GetTotbits(){

	return totbits;
}

int System::CmpControllers(DdManager* ddman,DdNode* cont_old)
{

	int isequal=0;
DdNode *tmp1;
tmp1=Cudd_addBddInterval(ddman,v, 1,1);
Cudd_Ref(tmp1);


if(Cudd_bddLeq(ddman,tmp1,cont_old) && Cudd_bddLeq(ddman,cont_old,tmp1))// && Cudd_addLeq(ddman,trans,trans_prev) && Cudd_addLeq(ddman,trans_prev,trans))
	{	
		isequal=1;
		printf("********equal controllers*****\n");
	}
else	
	printf("******** controllers are not the same !!!!!!!!!!!*****\n");
	Cudd_RecursiveDeref(ddman,tmp1);
	return isequal;
}

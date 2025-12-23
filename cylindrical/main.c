#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include <time.h>

void checkEnergyCons(Domain *D,int iteration,double *sumB,double *sumBz,double *sumBr,double *sumBp);

int main(int argc, char *argv[])
{
  	int i,j,k,n,s,iteration=0,boost,filterStep,labSaveStep;
  	int rnk,suddenDump=OFF,shiftIteration;
  	double factor,time_spent,t,x;
  	clock_t begin,end;
  	struct tm *t_now;
  	time_t timer; 	//measure time
  	char name[100];
  	FILE *out;
  	Domain D;  
  	LaserList *L;
  	PlasmaLens *PL;
  	LoadList *LL;
  	External Ext;
  	int myrank, nTasks;
  	MPI_Status status; 

  	MPI_Init(&argc,&argv);
  	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


  	begin=clock();

  	if(argc < 2) 
  	{  
  	  	printf("mpirun -np N jopic [inputFile] [dumpNum]\n"); 
  	  	exit(0); 
  	}

  	timer=time(NULL);
  	t_now=localtime(&timer);
  	if(myrank==0) {  
  	  	sprintf(name,"report");
  	  	out = fopen(name,"a");
  	  	fprintf(out,"simulation start.\n");
  	  	fprintf(out,"%d-%d-%d %d:%d:%d\n",t_now->tm_year+1900,t_now->tm_mon+1,t_now->tm_mday,t_now->tm_hour,t_now->tm_min,t_now->tm_sec);
  	} else ;

  	//parameter setting
  	parameterSetting(&D,&Ext,argv[1]);
  	if(argc >= 3) { 
  	  	D.dumpStep = atoi(argv[2]); 
  	  	if(D.dumpStart==D.dumpStep) D.dumpStart+=1; else;
  	} else;
	
	if(myrank==0)	    saveDensityProfile(&D); else ;
  	MPI_Barrier(MPI_COMM_WORLD);

  	//create mesh
  	boundary(&D,&Ext);
  	MPI_Barrier(MPI_COMM_WORLD);

  	//load plasma or load dump file
  	if(argc >= 3)  {   
  	  	iteration=D.dumpStep;
  	  	restoreDump(D,iteration);
  	  	t=D.dt*iteration;
  	  	sprintf(name,"dumpField%d.h5",iteration);
  	  	if(myrank==0) restoreIntMeta(name,"/minXDomain",&(D.minXDomain),1);
  	  	else ;
  	  	MPI_Bcast(&(D.minXDomain),1,MPI_INT,0,MPI_COMM_WORLD);
  	  	MPI_Barrier(MPI_COMM_WORLD);
  	  	D.maxXDomain=D.minXDomain+D.nx;
  	  	D.minXSub+=D.minXDomain;
  	}  else   {
  	  	LL=D.loadList;
  	  	s=0;
  	  	while(LL->next)      {
  	  	  	loadPlasma(&D,LL,s,iteration,D.istart,D.iend,D.jstart,D.jend);
  	  	  	LL=LL->next;
  	  	  	s++;
  	  	}
  	  	t=0;

  	}
  	MPI_Barrier(MPI_COMM_WORLD);

  	//pair charge
  	if(iteration==0) {
  	  	LL=D.loadList; s=0;
  	  	while(LL->next)      {
  	  	  	if(LL->pair==ON) solveCharge(&D,LL,D.RhoPairR,D.RhoPairI,D.istart,D.iend,D.jstart,D.jend,s,-1.0);  else ;
  	  	  	LL=LL->next; s++;
  	  	}
  	  	D.shareF[0]=D.RhoPairR;
  	  	D.shareF[1]=D.RhoPairI;
   	MPI_TransferJNew_Xplus(&D,2,D.nySub+5,3);
   	MPI_TransferJNew_Xminus(&D,2,D.nySub+5,3);
	 	if(D.Period==ON) { MPI_TransferJNew_Period_X(&D,2,D.nySub+5,3); } else ;   
  	  	D.RhoPairR=D.shareF[0];
  	  	D.RhoPairI=D.shareF[1];

  	  	D.shareF[0]=D.RhoPairR;
  	  	D.shareF[1]=D.RhoPairI;
   	MPI_TransferFNew_Xplus(&D,2,D.nySub+5,3);
   	MPI_TransferFNew_Xminus(&D,2,D.nySub+5,3);
	 	if(D.Period==ON) { MPI_TransferFNew_Period_X(&D,2,D.nySub+5,3); } else ;   
  	  	D.RhoPairR=D.shareF[0];
  	  	D.RhoPairI=D.shareF[1];
  	  	// if(D.L>1)  {
  	  	//   MPI_TransferDen_Xplus(&D,D.RhoPairR,D.RhoPairI,D.nySub+5,3);
  	  	//   MPI_TransferDen_Xminus(&D,D.RhoPairR,D.RhoPairI,D.nySub+5,3);
  	  	//   if(D.Period==ON) MPI_TransferDen_Period_X(&D,D.RhoNoPairR,D.RhoNoPairI,D.nySub+5,3); 
  	  	//   // MPI_Transfer2F_Xplus(&D,D.RhoPairR,D.RhoPairI,D.nySub+5,3);
  	  	//   // MPI_Transfer2F_Xminus(&D,D.RhoPairR,D.RhoPairI,D.nySub+5,3);
  	  	// }  else ;
  	} else ;

  	//shot laser
  	L=D.laserList;
  	while(L->next)  {
  	  if(L->loadMethod==Shot) {
  	    shotLaser(&D,L);
  	  } else ;
  	  L=L->next;
  	}

  	double sumB,sumBz,sumBx,sumBy;
	sumB=sumBz=sumBx=sumBy=0.0;

  	//rooping time 
  	while(iteration<=D.maxStep)
  	{

  	  	LL=D.loadList;
  	  	s=0;
  	  	while(LL->next)      {
  	  	  	loadBeam(&D,LL,s,iteration);
  	  	  	LL=LL->next;
  	  	  	s++;
  	  	}

  	  	//save dump File
  	  	if(D.dumpSave==ON && iteration>=D.dumpStart && iteration%D.dumpSaveStep==0) {
  	  	  	saveDump(D,iteration); 
  	  	  	end=clock();
  	  	  	time_spent=(end-begin)/CLOCKS_PER_SEC/60.0;
  	  	  	if(myrank==0) {
  	  	  	  	fprintf(out,"Time duration at %ddump:%.4gmin\n",iteration,time_spent);
  	  	  	  	printf("Time duration at %ddump:%.4gmin\n",iteration,time_spent);
  	  	  	} else ;
  	  	} else	;


  	  	//save File      
  	  	if(iteration%D.saveStep==0 && iteration>=D.saveStart) {
  	  	  	saveFile(D,iteration);
  	  	  	end=clock();
  	  	  	time_spent=(end-begin)/CLOCKS_PER_SEC/60.0;
  	  	  	if(myrank==0) {
  	  	  	  fprintf(out,"Time duration at %dth saveFile:%.4gmin\n",iteration,time_spent);
  	  	  	  printf("Time duration at %dth saveFile:%.4gmin\n",iteration,time_spent);
  	  	  	} else ;
		} else	;

  	  	//save center field      
  	  	if(iteration%D.centerStep==0) saveCenterField(&D,iteration); else	;
  	  	MPI_Barrier(MPI_COMM_WORLD);
  	  	//if(myrank==0) printf("saveCenterField at iteration=%d\n",iteration);


  	  	// redistributing particles lala
  	  	//if(D.redist==ON && iteration>0) particle_redist(&D,iteration,&Ext); else;

  	  	fieldSolve2(D,t,iteration);
  	  	//if(myrank==0) printf("fieldSolve2,iteration=%d\n",iteration);
		
  	  	solveF(D,iteration);
  	  	//if(myrank==0) printf("solveF,iteration=%d\n",iteration);


  	  	interpolation(&D,&Ext,iteration);
  	  	//if(myrank==0) printf("interpolation,iteration=%d\n",iteration);

  	  	// Plasma lens
		PL=D.lensList;
		while(PL->next)      {
			plasmaLens(&D,PL,iteration);
			PL=PL->next;
		}
  	  	//if(myrank==0) printf("plasmaLens,iteration=%d\n",iteration);

  	  	if(D.fieldIonization==ON) fieldIonization(&D,iteration); else;
  	  	//if(myrank==0) printf("fieldIonization,iteration=%d\n",iteration);

  	  	particlePush(&D,iteration);
  	  	//if(myrank==0) printf("particle Push,iteration=%d\n",iteration);


  	  	if(D.consCheck==ON) checkEnergyCons(&D,iteration,&sumB,&sumBz,&sumBx,&sumBy); else ;

  	  	updateCurrent(D,iteration);
  	  	//if(myrank==0) printf("current,iteration=%d\n",iteration);
  	  	//calConservation(D,iteration);

  	  	// Moving domain calculation
  	  	x=D.movingV*iteration+D.shiftStart;
  	  	if(D.moving==ON && x>D.maxXDomain-1)    
  	  	{
  	  	  	movingDomain(&D,iteration);
  	  	   //if(myrank==0) printf("movingDomain,iteration=%d\n",iteration);

  	  	  	if(myrank==D.L-1) {
  	  	  	  	LL=D.loadList; s=0;
  	  	  	  	while(LL->next)      {
  	  	  	  	  	loadPlasma(&D,LL,s,iteration,D.iend-1,D.iend,D.jstart,D.jend);
  	  	  	  	  	LL=LL->next; s++;
  	  	  	  	}
  	  	  	} else ; 
  	  	   //if(myrank==0) printf("loadPlamsa,iteration=%d\n",iteration);
         
  	  	  	rearrangeParticles(&D);
  	  	   //if(myrank==0) printf("rearrngeParticle,iteration=%d\n",iteration);
  	  	  	particleShareX(D);
  	  	   //if(myrank==0) printf("particleShareX,iteration=%d\n",iteration);
  	  	  	// if(D.Period==ON)   {
						//   MPI_TransferP_Period_X(&D,D.istart,D.iend,D.jstart-1,D.jend+1);
					  // } else ;

  	  	  	if(myrank==D.L-1) {
  	  	  	  	LL=D.loadList; s=0;
  	  	  	  	while(LL->next)      {
  	  	  	  	  	if(LL->pair==ON) 	solveCharge(&D,LL,D.RhoPairR,D.RhoPairI,D.iend-1,D.iend,D.jstart,D.jend,s,-1.0);	else ;
  	  	  	  	  	LL=LL->next;  s++;
  	  	  	  	}
  	  	      //if(myrank==0) printf("solveChare for RhoPair,iteration=%d\n",iteration);
  	  	  	} else ;
  	  	  	removeEdge(&D);
  	  	   //if(myrank==0) printf("removeEdge,iteration=%d\n",iteration);
  	  	} 

  	  	// Non-Moving domain calculation
  	  	else       
  	  	{
  	  	  	rearrangeParticles(&D);        
  	  	  	//printf("rearrange ,iteration=%d\n",iteration);

  	  	  	particleShareX(D);
  	  	  	if(D.Period==ON) {
  	  	  	  	MPI_TransferP_Period_X(&D,D.istart,D.iend,D.jstart-1,D.jend+1); 
  	  	  	  	rearrangeParticles_Period_Y(&D);
  	  	  	} else ;

  	  	  	removeEdge(&D);
  	  	  	//printf("removeEdge ,iteration=%d\n",iteration);

  	  	}  


  	  	fieldSolve1(D,t,iteration);
  	  	//printf("fieldSolve1,iteration=%d\n",iteration);

  	  	//time update
  	  	if(iteration%10==0 && myrank==0) printf("iteration = %d\n",iteration); else ;

  	  	iteration+=1;
  	  	t=D.dt*iteration;  

  	}     //end of time roop                  
  	// if(D.tracking==ON)  saveTracking(&D);

  	end=clock();
  	time_spent=(end-begin)/CLOCKS_PER_SEC;

  	//make 'report' file
  	if(myrank==0) {
  	  fprintf(out,"nx=%d, ",D.nx);
  	  fprintf(out,"ny=%d, ",D.ny);
  	  fprintf(out,"cores=%d, \n",nTasks);
  	  fprintf(out,"nSpecies=%d\n",D.nSpecies);
  	  fprintf(out,"running time=%.4gm\n",time_spent/60.0);
  	  fprintf(out,"\n");
  	  fclose(out);
  	} else	;

  	cleanMemory(&D);
	
  	MPI_Finalize();

  	return 0;
}

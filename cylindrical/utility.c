#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "mesh.h"
#include <math.h>
#include <mpi.h>


void checkEnergyCons(Domain *D,int iteration,double *sumB,double *sumBz,double *sumBx,double *sumBy)
{
	int i,j,s,m,ii,jj,nSpecies,numMode,istart,iend,jstart,jend;
	double tmp,z,x,y,r,R,invR,index,factor,weight,alpha,gamma,recv;
	double sumP,sumPz,sumPx,sumPy;
	double sumF,sumFz,sumFx,sumFy;
	double sinTh,cosTh;
	double EzR[D->numMode],ErR[D->numMode],EpR[D->numMode],BzR[D->numMode],BrR[D->numMode],BpR[D->numMode];
	double EzI[D->numMode],ErI[D->numMode],EpI[D->numMode],BzI[D->numMode],BrI[D->numMode],BpI[D->numMode];
	ptclList *p;
	LoadList *LL;
	FILE *out;

	istart=D->istart; iend=D->iend;
	jstart=D->jstart; jend=D->jend;
	nSpecies = D->nSpecies;
	numMode = D->numMode;

	int myrank, nTasks;

	MPI_Status status;

	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	sumP=sumPz=sumPx=sumPy=0.0;
	LL=D->loadList;
	s=0;
	while(LL->next)  {
      factor=M_PI*D->dr*D->dr*D->dz*D->lambda*D->lambda*D->lambda*LL->density;		

		for(i=istart; i<iend; i++)
			for(j=jstart+1; j<jend; j++)
			{
				p=D->particle[i][j].head[s]->pt;
				while(p) {
					gamma=sqrt(1.0+p->pz*p->pz+p->px*p->px+p->py*p->py);
					weight=p->weight*factor;

					sumP+=weight*(gamma-1)*LL->mass;
					sumPz+=weight*p->pz*LL->mass;
					sumPx+=weight*LL->mass*p->px;
					sumPy+=weight*LL->mass*p->py;

					p=p->next;
				}
			}

      LL=LL->next;
		s++;
	}
	
	sumF=sumFz=sumFx=sumFy=0.0;
 	factor=M_PI*D->dr*D->dr*D->dz*D->omega*D->omega*eMass*eps0/eCharge/eCharge*D->lambda*D->lambda*D->lambda;		
	if(D->fieldType==Split) {
		for(i=istart; i<iend; i++)
			for(j=jstart+1; j<jend; j++)
			{
				r=j-jstart;
				R=r+0.5;
				m=0;	
				ErR[m]=(D->PrR[m][i][j]+D->PlR[m][i][j])*0.5;
				EpR[m]=(D->SlR[m][i][j]+D->SrR[m][i][j])*0.5;
				BrR[m]=(D->SlR[m][i][j]-D->SrR[m][i][j])*0.5;
				BpR[m]=(D->PrR[m][i][j]-D->PlR[m][i][j])*0.5;
//				EzR[m]=D->EzR[m][i][j]*(r+0.75)/(2*r+1.0)+D->EzR[m][i][j+1]*(r+0.25)/(2*r+1.0);
//				BzR[m]=D->BzR[m][i][j]*(r+0.75)/(2*r+1.0)+D->BzR[m][i][j+1]*(r+0.25)/(2*r+1.0);
				EzR[m]=0.25*(D->EzR[m][i][j]+D->EzR[m][i][j+1]+D->EzR[m][i+1][j]+D->EzR[m][i+1][j+1]);
				BzR[m]=0.25*(D->BzR[m][i][j]+D->BzR[m][i][j+1]+D->BzR[m][i+1][j]+D->BzR[m][i+1][j+1]);
				for(m=1; m<D->numMode; m++) {
					ErR[m]=(D->PrR[m][i][j]+D->PlR[m][i][j])*0.5;
					EpR[m]=(D->SlR[m][i][j]+D->SrR[m][i][j])*0.5;
					BrR[m]=(D->SlR[m][i][j]-D->SrR[m][i][j])*0.5;
					BpR[m]=(D->PrR[m][i][j]-D->PlR[m][i][j])*0.5;
					EzR[m]=0.25*(D->EzR[m][i][j]+D->EzR[m][i][j+1]+D->EzR[m][i+1][j]+D->EzR[m][i+1][j+1]);
					BzR[m]=0.25*(D->BzR[m][i][j]+D->BzR[m][i][j+1]+D->BzR[m][i+1][j]+D->BzR[m][i+1][j+1]);

					ErI[m]=(D->PrI[m][i][j]+D->PlI[m][i][j])*0.5;
					EpI[m]=(D->SlI[m][i][j]+D->SrI[m][i][j])*0.5;
					BrI[m]=(D->SlI[m][i][j]-D->SrI[m][i][j])*0.5;
					BpI[m]=(D->PrI[m][i][j]-D->PlI[m][i][j])*0.5;
					EzI[m]=0.25*(D->EzI[m][i][j]+D->EzI[m][i][j+1]+D->EzI[m][i+1][j]+D->EzI[m][i+1][j+1]);
					BzI[m]=0.25*(D->BzI[m][i][j]+D->BzI[m][i][j+1]+D->BzI[m][i+1][j]+D->BzI[m][i+1][j+1]);
				}
				m=0;
				sumF+=2*(ErR[m]*ErR[m]+EpR[m]*EpR[m]+EzR[m]*EzR[m]+BrR[m]*BrR[m]+BpR[m]*BpR[m]+BzR[m]*BzR[m])*R*factor*0.5;
				sumFz+=2*(ErR[m]*BpR[m]-EpR[m]*BrR[m])*R*factor;
				for(m=1; m<D->numMode; m++) {
					sumF+=(ErR[m]*ErR[m]+EpR[m]*EpR[m]+EzR[m]*EzR[m]+BrR[m]*BrR[m]+BpR[m]*BpR[m]+BzR[m]*BzR[m])*R*factor*0.5;
					sumF+=(ErI[m]*ErI[m]+EpI[m]*EpI[m]+EzI[m]*EzI[m]+BrI[m]*BrI[m]+BpI[m]*BpI[m]+BzI[m]*BzI[m])*R*factor*0.5;
					sumFz+=(ErR[m]*BpR[m]-EpR[m]*BrR[m])*R*factor;
					sumFz+=(ErI[m]*BpI[m]-EpI[m]*BrI[m])*R*factor;
				}
				sumFx+=(ErI[2]*BzI[1]+ErI[1]*BzI[2]+(2*ErR[0]+ErR[2])*BzR[1]+(2*BzR[0]+BzR[2])*ErR[1])*R*factor;
				sumFx-=(EpR[1]*BzI[2]+EpI[2]*BzR[1]+(2*EpR[0]-EpR[2])*BzI[1]+(2*BzR[0]-BzR[2])*EpI[1])*R*factor;
				sumFx-=(EzI[2]*BrI[1]+EzI[1]*BrI[2]+(2*EzR[0]+EzR[2])*BrR[1]+(2*BrR[0]+BrR[2])*EzR[1])*R*factor;
				sumFx+=(EzR[1]*BpI[2]+EzI[2]*BpR[1]+(2*EzR[0]-EzR[2])*BpI[1]+(2*BpR[0]-BpR[2])*EzI[1])*R*factor;

				sumFy+=(BrI[2]*EzI[1]+BrI[1]*EzI[2]+(2*BrR[0]+BrR[2])*EzR[1]+(2*EzR[0]+EzR[2])*BrR[1])*R*factor;
				sumFy-=(BpR[1]*EzI[2]+BpI[2]*EzR[1]+(2*BpR[0]-BpR[2])*EzI[1]+(2*EzR[0]-EzR[2])*BpI[1])*R*factor;
				sumFy-=(BzI[2]*ErI[1]+BzI[1]*ErI[2]+(2*BzR[0]+BzR[2])*ErR[1]+(2*ErR[0]+ErR[2])*BzR[1])*R*factor;
				sumFy+=(BzR[1]*EpI[2]+BzI[2]*EpR[1]+(2*BzR[0]-BzR[2])*EpI[1]+(2*EpR[0]-EpR[2])*BzI[1])*R*factor;
			}

		if(myrank==0) {
			i=istart;
			for(j=jstart+1; j<jend; j++)
			{
				r=j-jstart;
				R=r+0.5;
				m=0;	
				ErR[m]=(0+D->PlR[m][i][j])*0.5;
				EpR[m]=(D->SlR[m][i][j]+0)*0.5;
				BrR[m]=(D->SlR[m][i][j]+0)*0.5;
				BpR[m]=(0-D->PlR[m][i][j])*0.5;
				EzR[m]=0.25*(D->EzR[m][i][j]+D->EzR[m][i][j+1]+D->EzR[m][i+1][j]+D->EzR[m][i+1][j+1]);
				BzR[m]=0.25*(D->BzR[m][i][j]+D->BzR[m][i][j+1]+D->BzR[m][i+1][j]+D->BzR[m][i+1][j+1]);
				for(m=1; m<D->numMode; m++) {
					ErR[m]=(0+D->PlR[m][i][j])*0.5;
					EpR[m]=(D->SlR[m][i][j]+0)*0.5;
					BrR[m]=(D->SlR[m][i][j]-0)*0.5;
					BpR[m]=(0-D->PlR[m][i][j])*0.5;
					EzR[m]=0.25*(D->EzR[m][i][j]+D->EzR[m][i][j+1]+D->EzR[m][i+1][j]+D->EzR[m][i+1][j+1]);
					BzR[m]=0.25*(D->BzR[m][i][j]+D->BzR[m][i][j+1]+D->BzR[m][i+1][j]+D->BzR[m][i+1][j+1]);

					ErI[m]=(0+D->PlI[m][i][j])*0.5;
					EpI[m]=(D->SlI[m][i][j]+0)*0.5;
					BrI[m]=(D->SlI[m][i][j]-0)*0.5;
					BpI[m]=(0-D->PlI[m][i][j])*0.5;
					EzI[m]=0.25*(D->EzI[m][i][j]+D->EzI[m][i][j+1]+D->EzI[m][i+1][j]+D->EzI[m][i+1][j+1]);
					BzI[m]=0.25*(D->BzI[m][i][j]+D->BzI[m][i][j+1]+D->BzI[m][i+1][j]+D->BzI[m][i+1][j+1]);
				}
				m=0;
				*sumB+=2*(ErR[m]*ErR[m]+EpR[m]*EpR[m]+EzR[m]*EzR[m]+BrR[m]*BrR[m]+BpR[m]*BpR[m]+BzR[m]*BzR[m])*R*factor*0.5;
				*sumBz+=2*(ErR[m]*BpR[m]-EpR[m]*BrR[m])*R*factor;
				for(m=1; m<D->numMode; m++) {
					*sumB+=(ErR[m]*ErR[m]+EpR[m]*EpR[m]+EzR[m]*EzR[m]+BrR[m]*BrR[m]+BpR[m]*BpR[m]+BzR[m]*BzR[m])*R*factor*0.5;
					*sumB+=(ErI[m]*ErI[m]+EpI[m]*EpI[m]+EzI[m]*EzI[m]+BrI[m]*BrI[m]+BpI[m]*BpI[m]+BzI[m]*BzI[m])*R*factor*0.5;
					*sumBz+=(ErR[m]*BpR[m]-EpR[m]*BrR[m])*R*factor;
					*sumBz+=(ErI[m]*BpI[m]-EpI[m]*BrI[m])*R*factor;
				}
				*sumBx+=(ErI[2]*BzI[1]+ErI[1]*BzI[2]+(2*ErR[0]+ErR[2])*BzR[1]+(2*BzR[0]+BzR[2])*ErR[1])*R*factor;
				*sumBx-=(EpR[1]*BzI[2]+EpI[2]*BzR[1]+(2*EpR[0]-EpR[2])*BzI[1]+(2*BzR[0]-BzR[2])*EpI[1])*R*factor;
				*sumBx-=(EzI[2]*BrI[1]+EzI[1]*BrI[2]+(2*EzR[0]+EzR[2])*BrR[1]+(2*BrR[0]+BrR[2])*EzR[1])*R*factor;
				*sumBx+=(EzR[1]*BpI[2]+EzI[2]*BpR[1]+(2*EzR[0]-EzR[2])*BpI[1]+(2*BpR[0]-BpR[2])*EzI[1])*R*factor;

				*sumBy+=(BrI[2]*EzI[1]+BrI[1]*EzI[2]+(2*BrR[0]+BrR[2])*EzR[1]+(2*EzR[0]+EzR[2])*BrR[1])*R*factor;
				*sumBy-=(BpR[1]*EzI[2]+BpI[2]*EzR[1]+(2*BpR[0]-BpR[2])*EzI[1]+(2*EzR[0]-EzR[2])*BpI[1])*R*factor;
				*sumBy-=(BzI[2]*ErI[1]+BzI[1]*ErI[2]+(2*BzR[0]+BzR[2])*ErR[1]+(2*ErR[0]+ErR[2])*BzR[1])*R*factor;
				*sumBy+=(BzR[1]*EpI[2]+BzI[2]*EpR[1]+(2*BzR[0]-BzR[2])*EpI[1]+(2*EpR[0]-EpR[2])*BzI[1])*R*factor;
			}
		}

	}	
	
	else if(D->fieldType==Yee) {
		for(i=istart; i<iend; i++)
			for(j=jstart+1; j<jend; j++)
			{
				r=j-jstart;
				m=0;	
				EzR[m]=(D->EzR[m][i][j]+D->EzR[m][i][j+1])*0.5;
				ErR[m]=(D->ErR[m][i][j]+D->ErR[m][i+1][j])*0.5;				
				EpR[m]=(D->EpR[m][i][j]+D->EpR[m][i+1][j]+D->EpR[m][i][j+1]+D->EpR[m][i+1][j+1])*0.25;				
				BzR[m]=(D->BzR[m][i][j]+D->BzR[m][i+1][j])*0.5;
				BrR[m]=(D->BrR[m][i][j]+D->BrR[m][i][j+1])*0.5;				
				BpR[m]=D->BpR[m][i][j];
				for(m=1; m<D->numMode; m++) {
					EzR[m]=(D->EzR[m][i][j]+D->EzR[m][i][j+1])*0.5;
					ErR[m]=(D->ErR[m][i][j]+D->ErR[m][i+1][j])*0.5;				
					EpR[m]=(D->EpR[m][i][j]+D->EpR[m][i+1][j]+D->EpR[m][i][j+1]+D->EpR[m][i+1][j+1])*0.25;				
					BzR[m]=(D->BzR[m][i][j]+D->BzR[m][i+1][j])*0.5;
					BrR[m]=(D->BrR[m][i][j]+D->BrR[m][i][j+1])*0.5;				
					BpR[m]=D->BpR[m][i][j];
					EzI[m]=(D->EzI[m][i][j]+D->EzI[m][i][j+1])*0.5;
					ErI[m]=(D->ErI[m][i][j]+D->ErI[m][i+1][j])*0.5;				
					EpI[m]=(D->EpI[m][i][j]+D->EpI[m][i+1][j]+D->EpI[m][i][j+1]+D->EpI[m][i+1][j+1])*0.25;				
					BzI[m]=(D->BzI[m][i][j]+D->BzI[m][i+1][j])*0.5;
					BrI[m]=(D->BrI[m][i][j]+D->BrI[m][i][j+1])*0.5;				
					BpI[m]=D->BpI[m][i][j];
				}
				m=0;
				sumF+=2*(ErR[m]*ErR[m]+EpR[m]*EpR[m]+EzR[m]*EzR[m]+BrR[m]*BrR[m]+BpR[m]*BpR[m]+BzR[m]*BzR[m])*R*factor*0.5;
				sumFz+=2*(ErR[m]*BpR[m]-EpR[m]*BrR[m])*R*factor;
				for(m=1; m<D->numMode; m++) {
					sumF+=(ErR[m]*ErR[m]+EpR[m]*EpR[m]+EzR[m]*EzR[m]+BrR[m]*BrR[m]+BpR[m]*BpR[m]+BzR[m]*BzR[m])*R*factor*0.5;
					sumF+=(ErI[m]*ErI[m]+EpI[m]*EpI[m]+EzI[m]*EzI[m]+BrI[m]*BrI[m]+BpI[m]*BpI[m]+BzI[m]*BzI[m])*R*factor*0.5;
					sumFz+=(ErR[m]*BpR[m]-EpR[m]*BrR[m])*R*factor;
					sumFz+=(ErI[m]*BpI[m]-EpI[m]*BrI[m])*R*factor;
				}
				sumFx+=(ErI[2]*BzI[1]+ErI[1]*BzI[2]+(2*ErR[0]+ErR[2])*BzR[1]+(2*BzR[0]+BzR[2])*ErR[1])*R*factor;
				sumFx-=(EpR[1]*BzI[2]+EpI[2]*BzR[1]+(2*EpR[0]-EpR[2])*BzI[1]+(2*BzR[0]-BzR[2])*EpI[1])*R*factor;
				sumFx-=(EzI[2]*BrI[1]+EzI[1]*BrI[2]+(2*EzR[0]+EzR[2])*BrR[1]+(2*BrR[0]+BrR[2])*EzR[1])*R*factor;
				sumFx+=(EzR[1]*BpI[2]+EzI[2]*BpR[1]+(2*EzR[0]-EzR[2])*BpI[1]+(2*BpR[0]-BpR[2])*EzI[1])*R*factor;

				sumFy+=(BrI[2]*EzI[1]+BrI[1]*EzI[2]+(2*BrR[0]+BrR[2])*EzR[1]+(2*EzR[0]+EzR[2])*BrR[1])*R*factor;
				sumFy-=(BpR[1]*EzI[2]+BpI[2]*EzR[1]+(2*BpR[0]-BpR[2])*EzI[1]+(2*EzR[0]-EzR[2])*BpI[1])*R*factor;
				sumFy-=(BzI[2]*ErI[1]+BzI[1]*ErI[2]+(2*BzR[0]+BzR[2])*ErR[1]+(2*ErR[0]+ErR[2])*BzR[1])*R*factor;
				sumFy+=(BzR[1]*EpI[2]+BzI[2]*EpR[1]+(2*BzR[0]-BzR[2])*EpI[1]+(2*EpR[0]-EpR[2])*BzI[1])*R*factor;

			}
	}


	if(myrank!=0) {
		for(i=1; i<nTasks; i++) {
			if(myrank==i)	MPI_Send(&sumF,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
			else ;
		}
	} else {
		for(i=1; i<nTasks; i++) {
			MPI_Recv(&recv,1,MPI_DOUBLE,i,i, MPI_COMM_WORLD,&status);
			sumF+=recv;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0) {
		for(i=1; i<nTasks; i++) {
			if(myrank==i)	MPI_Send(&sumFx,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
			else ;
		}
	} else {
		for(i=1; i<nTasks; i++) {
			MPI_Recv(&recv,1,MPI_DOUBLE,i,i, MPI_COMM_WORLD,&status);
			sumFx+=recv;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0) {
		for(i=1; i<nTasks; i++) {
			if(myrank==i)	MPI_Send(&sumFy,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
			else ;
		}
	} else {
		for(i=1; i<nTasks; i++) {
			MPI_Recv(&recv,1,MPI_DOUBLE,i,i, MPI_COMM_WORLD,&status);
			sumFy+=recv;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0) {
		for(i=1; i<nTasks; i++) {
			if(myrank==i)	MPI_Send(&sumFz,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
			else ;
		}
	} else {
		for(i=1; i<nTasks; i++) {
			MPI_Recv(&recv,1,MPI_DOUBLE,i,i, MPI_COMM_WORLD,&status);
			sumFz+=recv;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if(myrank!=0) {
		for(i=1; i<nTasks; i++) {
			if(myrank==i)	MPI_Send(&sumP,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
			else ;
		}
	} else {
		for(i=1; i<nTasks; i++) {
			MPI_Recv(&recv,1,MPI_DOUBLE,i,i, MPI_COMM_WORLD,&status);
			sumP+=recv;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0) {
		for(i=1; i<nTasks; i++) {
			if(myrank==i)	MPI_Send(&sumPx,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
			else ;
		}
	} else {
		for(i=1; i<nTasks; i++) {
			MPI_Recv(&recv,1,MPI_DOUBLE,i,i, MPI_COMM_WORLD,&status);
			sumPx+=recv;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0) {
		for(i=1; i<nTasks; i++) {
			if(myrank==i)	MPI_Send(&sumPy,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
			else ;
		}
	} else {
		for(i=1; i<nTasks; i++) {
			MPI_Recv(&recv,1,MPI_DOUBLE,i,i, MPI_COMM_WORLD,&status);
			sumPy+=recv;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0) {
		for(i=1; i<nTasks; i++) {
			if(myrank==i)	MPI_Send(&sumPz,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
			else ;
		}
	} else {
		for(i=1; i<nTasks; i++) {
			MPI_Recv(&recv,1,MPI_DOUBLE,i,i, MPI_COMM_WORLD,&status);
			sumPz+=recv;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);



	double unit,unitP;

	unit=eMass*velocityC*velocityC;
	unitP=eMass*velocityC;
	if(myrank==0) {
		out=fopen("energyCons","a");
		fprintf(out,"%d %g %g %g %g %g %g %g %g %g %g %g %g\n",iteration
			,sumP*unit		,sumF*unit		,*sumB*unit
			,sumPx*unitP	,sumFx*unitP	,*sumBx*unitP
			,sumPy*unitP	,sumFy*unitP	,*sumBy*unitP
			,sumPz*unitP	,sumFz*unitP	,*sumBz*unitP);
		fclose(out);
	} else ;

}




void saveDensityProfile(Domain *D)
{
   int i,j,NUM=1000,numX;
   double MIN=1e7,MAX=-1e7,minX,maxX,x,dx,posX,posY,neX,neFuncX;
   double realDz,realDr;
   double *density,*dataX;
   FILE *out;
   LoadList *LL;

   minX=MIN; maxX=MAX;
   realDz=D->dz*D->lambda;
   realDr=D->dr*D->lambda;

   LL=D->loadList;
   while(LL->next)      {
      switch (LL->type) {
         case Polygon :
            if(LL->xnodes>0) {
               for(i=0; i<LL->xnodes; i++) {
                  x=LL->xpoint[i];
                     if(x<minX)	minX=x;	else;
                     if(x>maxX)	maxX=x;	else;
               }
            }	else ;
         break;
      }
      LL=LL->next;
   }

   if(maxX==MAX || minX==MIN)	{
      maxX=0.0;
      minX=0.0;
   }	else ;

   numX=NUM;
   dx=(maxX-minX)/(1.0*numX);
   dataX=(double *)malloc((numX+1)*sizeof(double ));
   density=(double *)malloc((numX+1)*sizeof(double ));
   for(i=0; i<=numX; i++)	{
      dataX[i]=minX+i*dx;
      density[i]=0.0;
   }

   posY=0.0;
   for(j=0; j<=numX; j++) {
      posX=dataX[j];

      LL=D->loadList;
      while(LL->next)	{
         switch (LL->type)	{
         case Polygon :
            if(LL->xnodes>0) {			
               for(i=0; i<LL->xnodes-1; i++) {
                  if(posX>=LL->xpoint[i] && posX<LL->xpoint[i+1]) {
                     neX=((LL->xn[i+1]-LL->xn[i])/(LL->xpoint[i+1]-LL->xpoint[i])*(posX-LL->xpoint[i])+LL->xn[i]);
                     neFuncX=evaluate_rpn(LL->rpn_x[i],LL->rpn_size_x[i],posX*realDz,posY*realDr);
                     density[j]+=neX*neFuncX*LL->density;
			            i=LL->xnodes;
                  }	else ;
               }
            }	else ;
            break;
         }
         LL=LL->next;
      }
   }

   out=fopen("densityProfile","w");
   for(i=0; i<=numX; i++) 
      fprintf(out,"%g %g\n",dataX[i]*D->lambda*D->dz,density[i]);
   fclose(out);
   printf("\n\ndensityProfile is made.\n\n");

   free(dataX);
   free(density);
}


void RungeKutta(double *W,double *prob,int iter,double dt,int start,int end,int flag)
{
  int n,Z;
  double ddt,N,w,oldN,oldW,oldWN,k1,k2,k3,k4;

  ddt=dt/((double)iter);
  //initialization of prob[Z]
  prob[start]=1.0; for(Z=start+1; Z<end; Z++) prob[Z]=0.0;

  n=0;
  while(n<iter)
  {
    //at start
    N=prob[start], w=W[start];
    k1=-w*N;
    k2=-w*(N+0.5*ddt*k1);
    k3=-w*(N+0.5*ddt*k2);
    k4=-w*(N+ddt*k3);   
    prob[start]=N+ddt/6.0*(k1+2.0*k2+2.0*k3+k4); 
    //rest of ..
    Z=start+1;
    while(Z<end) {
      oldWN=N*w; N=prob[Z], w=W[Z];
      k1=oldWN-w*N;
      k2=oldWN-w*(N+0.5*ddt*k1);
      k3=oldWN-w*(N+0.5*ddt*k2);
      k4=oldWN-w*(N+ddt*k3);
      prob[Z]=N+ddt/6.0*(k1+2.0*k2+2.0*k3+k4); 
      Z++;
    }
//    //last
//    k1=w*N;
//    k2=w*(N+0.5*ddt*k1);
//    k3=w*(N+0.5*ddt*k2);
//    k4=w*(N+ddt*k3);   
//    prob[end]=N+ddt/6.0*(k1+2.0*k2+2.0*k3+k4); 
    n++;
  }

  Z=start+1;
  while(Z<end) {
   prob[Z]+=prob[Z-1];
   Z++;
  }
//  if(flag==1) {
//    for(Z=start; Z<=end; Z++) 
//    printf("prob[%d]=%g, ",Z,prob[s][Z]);
//    printf("\n");
//  }
}

double randomValue(double beta) 
{ 
   double r; 
   int intRand, randRange=10000000, rangeDev; 
 
   rangeDev=(int)(randRange*(1.0-beta)); 
   intRand = rand() % (randRange-rangeDev); 
   r = ((double)intRand)/randRange+(1.0-beta); 
 
   return r; 
} 


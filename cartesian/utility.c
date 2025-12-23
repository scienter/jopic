#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>




void checkEnergyCons(Domain *D,int iteration)
{
	int i,j,k,s,nSpecies,istart,iend,jstart,jend,kstart,kend;
	double tmp,z,x,y,index,factor,weight,alpha,gamma,recv;
	double sumP,sumPz,sumPx,sumPy;
	double sumF,sumFz,sumFx,sumFy;
	double Ex,Ey,Ez,Bx,By,Bz;
	ptclList *p;
	LoadList *LL;
	FILE *out;

	istart=D->istart; iend=D->iend;
	jstart=D->jstart; jend=D->jend;
	kstart=D->kstart; kend=D->kend;
	nSpecies = D->nSpecies;

	int myrank, nTasks;

	MPI_Status status;

	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	sumP=sumPz=sumPx=sumPy=0.0;
	LL=D->loadList;
	s=0;
	while(LL->next)  {
      factor=D->dx*D->dy*D->dz*D->lambda*D->lambda*D->lambda*LL->density;		

		for(i=istart; i<iend; i++)
			for(j=jstart; j<jend; j++)
				for(k=kstart; k<kend; k++)
				{
					p=D->particle[i][j][k].head[s]->pt;
					while(p) {
						gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
						weight=p->weight*factor;

						sumP+=weight*(gamma-1)*LL->mass;
						sumPx+=weight*LL->mass*p->p1;
						sumPy+=weight*LL->mass*p->p2;
						sumPz+=weight*LL->mass*p->p3;

						p=p->next;
					}
				}

      LL=LL->next;
		s++;
	}
	
	sumF=0.0;	
	factor=D->dx*D->dy*D->dz*D->omega*D->omega*eMass*eps0/eCharge/eCharge*D->lambda*D->lambda*D->lambda;		
	if(D->fieldType==Split) {
		for(i=istart; i<iend; i++)
			for(j=jstart; j<jend; j++)
				for(k=kstart; k<kend; k++)
				{
					Ey=0.5*(D->Pr[i][j][k]+D->Pl[i][j][k]);
					Ez=0.5*(D->Sl[i][j][k]+D->Sr[i][j][k]);
					Ex=D->Ex[i][j][k];
					Bz=0.5*(D->Pr[i][j][k]-D->Pl[i][j][k]);
					By=0.5*(D->Sl[i][j][k]-D->Sr[i][j][k]);
					Bx=D->Bx[i][j][k];
					sumF+=(Ey*Ey+Ez*Ez+Ex*Ex+By*By+Bz*Bz+Bx*Bx)*factor;
				}
	}	
	else if(D->fieldType==Yee) {
		for(i=istart; i<iend; i++)
			for(j=jstart; j<jend; j++)
				for(k=kstart; k<kend; k++)
				{
					Ey=D->Ey[i][j][k];
					Ez=D->Ez[i][j][k];
					Ex=D->Ex[i][j][k];
					Bz=D->Bz[i][j][k];
					By=D->By[i][j][k];
					Bx=D->Bx[i][j][k];
					sumF+=(Ey*Ey+Ez*Ez+Ex*Ex+By*By+Bz*Bz+Bx*Bx)*factor;
				}
	}


	if(nTasks>1) {
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
	} else ;

	double unit,unitP;

	unit=eMass*velocityC*velocityC;
	unitP=eMass*velocityC;

   if(myrank==0) {
      out=fopen("energyCons","a");
      fprintf(out,"%d %g %g\n",iteration
         ,sumP*unit     ,sumF*unit);
		fclose(out);
	} else ;	

//	D->consVal[iteration][0]=sumP*unit;
//	D->consVal[iteration][1]=sumF;

}









void saveDensityProfile(Domain *D)
{
	int i,j,NUM=100,numX;
	double MIN=1e7,MAX=-1e7,minX,maxX,x,dx,posX,ne,neFunc,realDx;
	double *density,*dataX;
	LoadList *LL;
   FILE *out;

	minX=MIN; maxX=MAX;
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
   for(i=0; i<=numX; i++)  {
      dataX[i]=minX+i*dx;
      density[i]=0.0;
   }
   realDx=D->dx*D->lambda;

   for(j=0; j<=numX; j++) {
      posX=dataX[j];

      LL=D->loadList;
      while(LL->next)   {
         switch (LL->type) {
         case Polygon :
            if(LL->xnodes>0) {
               for(i=0; i<LL->xnodes-1; i++) {
                  if(posX>=LL->xpoint[i] && posX<LL->xpoint[i+1]) {
                  	ne=((LL->xn[i+1]-LL->xn[i])/(LL->xpoint[i+1]-LL->xpoint[i])*(posX-LL->xpoint[i])+LL->xn[i]);
                     neFunc = evaluate_rpn(LL->rpn_x[i], LL->rpn_size_x[i], posX*realDx);

                     density[j]+=ne*neFunc*LL->density;
                     i=LL->xnodes;
                  }  else ;
               }
            }  else ;
            break;
         }
         LL=LL->next;
      }
   }

   out=fopen("densityProfile","w");
   for(i=0; i<=numX; i++)
      fprintf(out,"%g %g\n",dataX[i]*realDx,density[i]);
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


#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <complex.h>
#include "mpi.h"

void loadLaser_Gaussian_2D_Yee_Pukhov(Domain *D,LaserList *L,double t);
void loadLaser_Gaussian_3D_Yee_Pukhov(Domain *D,LaserList *L,double t);
void loadLaser_Gaussian_2D_Split(Domain *D,LaserList *L,double t);
void loadLaser_Gaussian_3D_Split(Domain *D,LaserList *L,double t);
void loadLaser_SSTF_2D_Yee_Pukhov(Domain *D,LaserList *L,double t);
void loadLaser_SSTF_3D_Yee_Pukhov(Domain *D,LaserList *L,double t);
void loadLaser_SSTF_2D_Split(Domain *D,LaserList *L,double t);
void loadLaser_SSTF_3D_Split(Domain *D,LaserList *L,double t);
void loadLaser1D_Split(Domain *D,LaserList *L,double t);
void loadLaser1D_Yee_Pukhov(Domain *D,LaserList *L,double t);


void loadLaser(Domain *D,LaserList *L,double t)
{

  double laserDur;

  if(L->mode==Gaussian)
	  laserDur=(L->rU*2+L->rD*2+L->flat+L->retard)*D->divisionLambda;
   else if(L->mode==SSTF)
	  laserDur=(16.0/L->sigOmegaSSTF*velocityC/D->lambda+L->retard)*D->divisionLambda;
	else laserDur=0.0;


  //if(D->boostOn==OFF && t<D->nx*D->dt)
  if(D->boostOn==OFF && t<laserDur*D->dt)
  {
    switch((D->fieldType-1)*3+D->dimension)  {
    case (Split-1)*3+1 :
      loadLaser1D_Split(D,L,t);
      break;
    case (Split-1)*3+2 :
      if(L->mode==Gaussian) loadLaser_Gaussian_2D_Split(D,L,t);
      else if(L->mode==SSTF) loadLaser_SSTF_2D_Split(D,L,t);
      else ;
      break;
    case (Split-1)*3+3 :
      if(L->mode==Gaussian) loadLaser_Gaussian_3D_Split(D,L,t);
      else if(L->mode==SSTF) loadLaser_SSTF_3D_Split(D,L,t);
      else;
      break;
    case (Yee-1)*3+1 :
    case (Pukhov-1)*3+1 :
      loadLaser1D_Yee_Pukhov(D,L,t);
      break;
    case (Yee-1)*3+2 :
    case (Pukhov-1)*3+2 :
      if(L->mode==Gaussian) loadLaser_Gaussian_2D_Yee_Pukhov(D,L,t);
      else if(L->mode==SSTF) loadLaser_SSTF_2D_Yee_Pukhov(D,L,t);
      else ;
      break;
    case (Yee-1)*3+3 :
    case (Pukhov-1)*3+3 :
      if(L->mode==Gaussian) loadLaser_Gaussian_3D_Yee_Pukhov(D,L,t);
      else if(L->mode==SSTF) loadLaser_SSTF_3D_Yee_Pukhov(D,L,t);
      else ;
      break;
    default :
      printf("In loadLaser, what is field_type? and what is dimension?\n");
    }
  } else ;
}

void loadLaser1D_Yee_Pukhov(Domain *D,LaserList *L,double t)
{
   double rU,rD,longitudinal,t0,flat,omega,amp,dt,delT;
   int istart,iend,positionX,j,k,laserOK=0;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;
   iend=D->iend;
   dt=D->dt/D->dtRatio;  
 
   rU=L->rU*D->divisionLambda*dt;
   rD=L->rD*D->divisionLambda*dt;
   flat=L->flat*D->divisionLambda*L->lambda/D->lambda*dt;

   t0=2.0*rU;
   omega=2.0*M_PI*L->omega/D->omega*(1.0+L->gdd/t0*(t-t0));

   if(t<2.0*rU)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2.0*rU && t<2.0*rU+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2.0*rU+flat && t<2.0*rU+flat+2.0*rD) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2.0*rU+2.0*rD+flat) 
      longitudinal=0.0;

   positionX=L->loadPointX;
   if(positionX>=D->minXSub && positionX<D->maxXSub)
     laserOK=ON;
   else ;
   positionX=L->loadPointX+istart-D->minXSub;

   if(laserOK==1)   {
     switch (L->add)  {
     case OFF :
       j=k=0;
       if(L->polarity==2)     {
         amp=longitudinal*sin(omega*t);
         D->Ey[positionX][j][k]=amp;            
         D->Ey[positionX][j][k]=amp;           
       } else if(L->polarity==3)   {
         amp=longitudinal*sin(omega*t);
         D->Ez[positionX][j][k]=amp;            
         D->Ez[positionX][j][k]=amp;           
       } else ;
       break;

     case ON : //adding
       j=k=0;
       if(L->polarity==2)     {
         amp=longitudinal*sin(omega*t);
         D->Ey[positionX][j][k]+=amp;            
         D->Ey[positionX][j][k]+=amp;           
       } else if(L->polarity==3)   {
         amp=longitudinal*sin(omega*t);
         D->Ez[positionX][j][k]+=amp;            
         D->Ez[positionX][j][k]+=amp;           
       } else ;
       break;
     }
   }  else ;   //End of field is OK

/*
   delT=dt*0.5;
   if(t<2.0*rU)
      longitudinal=L->amplitude*exp(-(t-t0+delT)*(t-t0+delT)/rU/rU);
   else if(t>=2.0*rU && t<2.0*rU+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2.0*rU+flat && t<2.0*rU+flat+2.0*rD) 
      longitudinal=L->amplitude*exp(-(t-t0+delT)*(t-t0+delT)/rD/rD);
   else if(t>=2.0*rU+2.0*rD+flat) 
      longitudinal=0.0;

   positionX=L->loadPointX;
   if(positionX>=D->minXSub && positionX<D->maxXSub)
     laserOK=1;
   else ;
   positionX=L->loadPointX+istart-D->minXSub;

   if(laserOK==1)   {
     j=k=0;
     if(L->polarity==2)     {
       amp=longitudinal*sin(omega*(t+delT));
       D->Bz[positionX][j][k]=amp;            
       D->Bz[positionX][j][k]=amp;           
     } else if(L->polarity==3)   {
       amp=longitudinal*sin(omega*(t+delT));
       D->By[positionX][j][k]=amp;            
       D->By[positionX][j][k]=amp;           
     } else ;
   }  else ;   //End of field is OK
  */
}

void loadLaser1D_Split(Domain *D,LaserList *L,double t)
{
   double rU,rD,longitudinal,t0,flat,omega,amp,factR,factL;
   int istart,iend,positionX,j,k,laserOK=0;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;
   iend=D->iend;

   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;
   flat=L->flat*D->divisionLambda*D->dt*L->lambda/D->lambda;

   t0=2.0*rU;
   omega=2.0*M_PI*L->omega/D->omega;

   if(t<2.0*rU)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2.0*rU && t<2.0*rU+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2.0*rU+flat && t<2.0*rU+flat+2.0*rD) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2.0*rU+2.0*rD+flat) 
      longitudinal=0.0;

   positionX=L->loadPointX;
   if(positionX>=D->minXSub && positionX<D->maxXSub)
     laserOK=ON;
   else ;
   positionX=L->loadPointX+istart-D->minXSub;

   if(L->direction==RIGHT) { factR=1.0; factL=0.0; }
   else if(L->direction==LEFT) { factR=0.0; factL=1.0; }
   else if(L->direction==BOTH) { factR=1.0; factL=1.0; }
   else { factR=0.0; factL=0.0; }

   if(laserOK==ON)   {
     switch (L->add)  {
     case OFF :
       j=k=0;
       if(L->polarity==2)     {
         amp=longitudinal*sin(omega*t);
         D->Pr[positionX][j][k]=amp*factR;            
         D->Pl[positionX][j][k]=amp*factL;           
       } else if(L->polarity==3)   {
         amp=longitudinal*sin(omega*t);
         D->Sr[positionX][j][k]=amp*factR;            
         D->Sl[positionX][j][k]=amp*factL;           
       } else ;
       break;

     case ON : //adding
       j=k=0;
       if(L->polarity==2)     {
         amp=longitudinal*sin(omega*t);
         D->Pr[positionX][j][k]+=amp*factR;            
         D->Pl[positionX][j][k]+=amp*factL;           
       } else if(L->polarity==3)   {
         amp=longitudinal*sin(omega*t);
         D->Sr[positionX][j][k]+=amp*factR;            
         D->Sl[positionX][j][k]+=amp*factL;           
       } else ;
       break;
     }
   }  else ;   //End of field is OK
}


void loadLaser_Gaussian_2D_Yee_Pukhov(Domain *D,LaserList *L,double t)
{
   double rU,rD,longitudinal,t0,flat;
   double zR,w0,w,phi,omega,kx,pphi,amp,dt;
   double x,y,z,r2,w2;
   int istart,iend,jstart,jend,kstart,kend,minj,maxj;
   int positionX,rank,i,j,k,jC,kC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   dt=D->dt/D->dtRatio;

   rU=L->rU*D->divisionLambda*dt;
   rD=L->rD*D->divisionLambda*dt;
   flat=L->flat*D->divisionLambda*L->lambda/D->lambda*dt;

   jC=L->loadPointY;	//y position

   t0=2*rU+L->retard;
   zR=L->rayleighLength;
   w0=L->beamWaist;  
   if(L->direction==RIGHT)
     x=L->loadPointX*D->dx-L->focus;
   else if(L->direction==LEFT)
     x=L->focus-L->loadPointX*D->dx;
   else ;
   w=w0*sqrt(1.0+x*x/zR/zR);
   phi=atan(x/zR);
   omega=2*M_PI*L->omega/D->omega;
   kx=2*M_PI*D->lambda/L->lambda;

   if(t<t0)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=t0 && t<t0+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=t0+flat && t<t0+flat+2*rD) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=t0+2*rD+flat) 
      longitudinal=0.0;

//   positionX=L->loadPointX+istart-D->minXSub;
   positionX=L->loadPointX;
   if(positionX>=D->minXSub && positionX<D->maxXSub) {
     laserOK=1;
     positionX+=istart-D->minXSub;
   }   else ;

//   if(positionX>D->minXSub && positionX<=D->maxXSub &&
//      jC-D->minYSub>=jstart && jC-D->minYSub<jend &&
//      kC-D->minZSub>=kstart && kC-D->minZSub<kend)
//     laserOK=1;

   if(laserOK==1 && L->add==OFF)  {
     if(L->polarity==2) {
       k=0;
       w2=w*w;
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         r2=y*y;
         pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
         D->Ey[positionX][j][k]=amp;            
       }
     } else if(L->polarity==3) {
       k=0;
       w2=w*w;
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         r2=y*y;
         pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
         D->Ez[positionX][j][k]=amp;            
       }
     }  
   } else if(laserOK==1 && L->add==ON)  {
     if(L->polarity==2) {
       k=0;
       w2=w*w;
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         r2=y*y;
         pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
         D->Ey[positionX][j][k]+=amp;            
       }
     } else if(L->polarity==3) {
       k=0;
       w2=w*w;
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         r2=y*y;
         pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
         D->Ez[positionX][j][k]+=amp;            
       }
     }  
   }  else ;   //End of field is OK
}
void loadLaser_Gaussian_2D_Split(Domain *D,LaserList *L,double t)
{
   double rU,rD,longitudinal,t0,flat;
   double zR,w0,w,phi,omega,kx,pphi;
   double x,y,z,r2,w2,retard;
	double complex amp;
   int istart,iend,jstart,jend,kstart,kend,minj,maxj;
   int positionX,rank,j,k,jC,kC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;  k=0;

   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;
   flat=L->flat*D->divisionLambda*D->dt*L->lambda/D->lambda;
   retard=L->retard;

   jC=L->loadPointY;	//y position

   t0=2*rU+retard;
   zR=L->rayleighLength;
   w0=L->beamWaist;  
   if(L->direction==RIGHT)
     x=L->loadPointX*D->dx-L->focus;
   else if(L->direction==LEFT)
     x=L->focus-L->loadPointX*D->dx;
   else ;
   w=w0*sqrt(1.0+x*x/zR/zR);
   phi=atan(x/zR);
   omega=2*M_PI*L->omega/D->omega;
   kx=2*M_PI*D->lambda/L->lambda;

   if(t<2*rU+retard)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2*rU+retard && t<2*rU+flat+retard) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2*rU+flat+retard && t<2*rU+flat+2*rD+retard) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2*rU+2*rD+flat+retard) 
      longitudinal=0.0;

   positionX=L->loadPointX;
   if(positionX>=D->minXSub && positionX<D->maxXSub) {
     laserOK=ON;
     positionX+=istart-D->minXSub;
   }   else ;

   if(laserOK==ON)  {
     switch (L->add)  {
     case OFF :
       w2=w*w;
       if(L->polarity==2)       {
         for(j=0; j<jend+3; j++)         {
           y=(j-jstart+D->minYSub-jC)*D->dy;   r2=y*y;
           pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
           amp=sqrt(w0/w)*exp(-r2/w2)*cexp(I*pphi);
           D->Pr[positionX][j][k]=longitudinal*creal(amp);            
//           pphi=-x/zR*r2/w2+0.5*phi-kx*x-omega*t;
//           amp=longitudinal*sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
           D->Pl[positionX][j][k]=longitudinal*creal(amp);            
         }
       } 
       else if(L->polarity==3)       {
         for(j=jstart; j<jend; j++)         {
           y=(j-jstart+D->minYSub-jC)*D->dy;   r2=y*y;
           pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
           amp=sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
           D->Sr[positionX][j][k]=longitudinal*creal(amp);            
           D->Sl[positionX][j][k]=longitudinal*creal(amp);            
         }
       }
       break;

     case ON : //adding
       if(L->polarity==2)    {
         for(j=jstart; j<jend; j++)    {
           y=(j-jstart+D->minYSub-jC)*D->dy;    r2=y*y;
           pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
           amp=sqrt(w0/w)*exp(-r2/w2)*cexp(I*pphi);
           D->Pr[positionX][j][k]+=longitudinal*creal(amp);            
           D->Pl[positionX][j][k]+=longitudinal*creal(amp);            
         }
       }
       else if(L->polarity==3)       {
         for(j=jstart; j<jend; j++)    {
           y=(j-jstart+D->minYSub-jC)*D->dy;   r2=y*y;
           pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
           amp=sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
           D->Sr[positionX][j][k]+=longitudinal*creal(amp);            
           D->Sl[positionX][j][k]+=longitudinal*creal(amp);            
         }
       }
       break;
     default :
       printf("In loadLaser.c, what laser mode?\n");
       break;
     }		//End of switch
   }  else ;   	//End of field is OK
}


void loadLaser_SSTF_2D_Yee_Pukhov(Domain *D,LaserList *L,double t)
{
   double t0,omega,amp,dy,y,alpha,minY,f12,shiftT;
   double sigW,minT,sigY,f0,f,x,k0,a0,lambda0,phase;
   double complex a_k0,coef,arg1,compA,m,n,tmp,gaussTmp;
   int istart,iend,jstart,jend,kstart,kend,minj,maxj;
   int positionX,rank,i,j,k,jC,kC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;

   dy=D->dy;
   minY=D->minYDomain*dy*D->lambda;

   jC=L->loadPointY;	//y position

   sigW=L->sigOmegaSSTF;  
   sigY=L->sigYSSTF;  
   alpha=L->alphaYSSTF;
//   minT=-4.0/sigW*sqrt((1.0+alpha*alpha*sigW*sigW/sigY/sigY))-L->retard*D->lambda/velocityC;
//   minT=-8.0/sigW-L->retard*D->dt*D->divisionLambda*2.0*M_PI/D->omega;
   minT=-4.0/sigW*sqrt((1.0+alpha*alpha*sigW*sigW/sigY/sigY))-L->retard*D->dt*D->divisionLambda*2.0*M_PI/D->omega;

   f0=L->lensFocus;
   f=L->focus*D->lambda;
   x=f0-f;
   k0=2*M_PI/L->lambda;
   a0=L->amplitude;
   lambda0=D->lambda;
   omega=D->omega;

   f12=f0*f0*k0*k0*sigY*sigY*sigY*sigY/(4.0*f0*f0+k0*k0*sigY*sigY*sigY*sigY);

   // for normalization at the focus
   a_k0=f12/k0/k0/sigY/sigY-I*0.5/k0*(f0-f12/f0);
   gaussTmp=sigW*sigY*0.5*csqrt(M_PI/a_k0/(1-I*k0*sigY*sigY*0.5/f0));

   a_k0=f12/k0/k0/sigY/sigY-I*0.5/k0*(x-f12/f0);
   m=1.0+alpha*alpha*sigW*sigW*(x-f0)*(x-f0)/4.0/f0/f0/a_k0-I*k0*alpha*alpha*sigW*sigW*(x-f0)/2.0/f0/f0;
   n=k0*alpha*sigW/f0+I*alpha*sigW*(x-f0)/2.0/f0/a_k0;
   tmp=M_PI/m/a_k0/(1-I*k0*sigY*sigY*0.5/f0);
   tmp=csqrt(tmp)*sigY*sigW*0.5;
   coef=tmp/gaussTmp*a0;			//normalization factor at z=x
//	printf("tmp=%g,nncoef1=%g,coef=%g\n",cabs(tmp),cabs(coef1),coef);

   //time
   t0=minT+t*2.0*M_PI/D->omega;

//   positionX=L->loadPointX+istart-D->minXSub;
   positionX=L->loadPointX;
   if(positionX>=D->minXSub && positionX<D->maxXSub) {
     laserOK=1;
     positionX+=istart-D->minXSub;
   }   else ;

   if(laserOK==1 && L->add==OFF)  {
     if(L->polarity==2) {
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*dy*lambda0;

         arg1=-0.25*y*y/a_k0-0.25/m*(sigW*t0+n*y)*(sigW*t0+n*y)+I*omega*t0;
         compA=coef*cexp(arg1);
//         amp=cabs(compA);
//	 phase=carg(compA);
//         D->Ey[positionX][j][0]=amp*cos(phase); 
         D->Ey[positionX][j][0]=creal(compA); 
       }
     }
/*	  
	   else if(L->polarity==3) {
       k=0;
       w2=w*w;
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         r2=y*y;
         pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
         D->Ez[positionX][j][k]=amp;            
       }
     }  
*/
   } 
   else if(laserOK==1 && L->add==ON)  {
     if(L->polarity==2) {
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*dy*lambda0;

         arg1=-0.25*y*y/a_k0-0.25/m*(sigW*t0+n*y)*(sigW*t0+n*y)+I*omega*t0;
         compA=coef*cexp(arg1);
         amp=cabs(compA);
	 phase=carg(compA);
//			phase=arg1;
         D->Ey[positionX][j][0]+=amp*cos(phase); 
       }
     }
   } 
   else ;   //End of field is OK
}

void loadLaser_SSTF_2D_Split(Domain *D,LaserList *L,double t)
{
   double t0,omega,amp,dy,y,alpha,minY,f12,shiftT;
   double sigW,minT,sigY,f0,f,x,k0,a0,lambda0,phase;
   double complex a_k0,coef,arg1,compA,m,n,tmp,gaussTmp;
   int istart,iend,jstart,jend,kstart,kend,minj,maxj;
   int positionX,rank,i,j,k,jC,kC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;

   dy=D->dy;
   minY=D->minYDomain*dy*D->lambda;

   jC=L->loadPointY;	//y position

   sigW=L->sigOmegaSSTF;  
   sigY=L->sigYSSTF;  
   alpha=L->alphaYSSTF;
//   minT=-4.0/sigW*sqrt((1.0+alpha*alpha*sigW*sigW/sigY/sigY))-L->retard*D->lambda/velocityC;
//   minT=-8.0/sigW-L->retard*D->dt*D->divisionLambda*2.0*M_PI/D->omega;
   minT=-4.0/sigW*sqrt((1.0+alpha*alpha*sigW*sigW/sigY/sigY))-L->retard*D->dt*D->divisionLambda*2.0*M_PI/D->omega;

   f0=L->lensFocus;
   f=L->focus*D->lambda;
   x=f0-f;
   k0=2*M_PI/L->lambda;
   a0=L->amplitude;
   lambda0=D->lambda;
   omega=D->omega;

   f12=f0*f0*k0*k0*sigY*sigY*sigY*sigY/(4.0*f0*f0+k0*k0*sigY*sigY*sigY*sigY);

   // for normalization at the focus
   a_k0=f12/k0/k0/sigY/sigY-I*0.5/k0*(f0-f12/f0);
   gaussTmp=sigW*sigY*0.5*csqrt(M_PI/a_k0/(1-I*k0*sigY*sigY*0.5/f0));

   a_k0=f12/k0/k0/sigY/sigY-I*0.5/k0*(x-f12/f0);
   m=1.0+alpha*alpha*sigW*sigW*(x-f0)*(x-f0)/4.0/f0/f0/a_k0-I*k0*alpha*alpha*sigW*sigW*(x-f0)/2.0/f0/f0;
   n=k0*alpha*sigW/f0+I*alpha*sigW*(x-f0)/2.0/f0/a_k0;
   tmp=M_PI/m/a_k0/(1-I*k0*sigY*sigY*0.5/f0);
   tmp=csqrt(tmp)*sigY*sigW*0.5;
   coef=tmp/gaussTmp*a0;			//normalization factor at z=x

   //time
   t0=minT+t*2.0*M_PI/D->omega;

   positionX=L->loadPointX;
   if(positionX>=D->minXSub && positionX<D->maxXSub) {
     laserOK=1;
     positionX+=istart-D->minXSub;
   }   else ;

   if(laserOK==1 && L->add==OFF)  {
     if(L->polarity==2) {
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*dy*lambda0;

         arg1=-0.25*y*y/a_k0-0.25/m*(sigW*t0+n*y)*(sigW*t0+n*y)+I*omega*t0;
         compA=coef*cexp(arg1);
//         amp=cabs(compA);
//	 phase=carg(compA);
//         D->Ey[positionX][j][0]=amp*cos(phase); 
         D->Pr[positionX][j][0]=creal(compA); 
         D->Pl[positionX][j][0]=creal(compA); 
       }
     }
   } 
   else if(laserOK==1 && L->add==ON)  {
     if(L->polarity==2) {
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*dy*lambda0;

         arg1=-0.25*y*y/a_k0-0.25/m*(sigW*t0+n*y)*(sigW*t0+n*y)+I*omega*t0;
         compA=coef*cexp(arg1);
         amp=cabs(compA);
	 phase=carg(compA);
         D->Pr[positionX][j][0]+=creal(compA); 
         D->Pl[positionX][j][0]+=creal(compA); 
       }
     }
   } 
   else ;   //End of field is OK
}

void loadLaser_Gaussian_3D_Yee_Pukhov(Domain *D,LaserList *L,double t)
{
   double rU,rD,longitudinal,t0,flat;
   double zR,w0,w,phi,omega,kx,pphi,amp,dt;
   double x,y,z,r2,w2,w2p,elliptic;
   int istart,iend,jstart,jend,kstart,kend,minj,maxj;
   int positionX,rank,i,j,k,jC,kC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   dt=D->dt/D->dtRatio;

   rU=L->rU*D->divisionLambda*dt;
   rD=L->rD*D->divisionLambda*dt;
   flat=L->flat*D->divisionLambda*L->lambda/D->lambda*dt;

   jC=L->loadPointY;	//y position
   kC=L->loadPointZ;	//z position

   t0=2*rU+L->retard;
   x=-L->focus;   
   zR=L->rayleighLength;
   w0=L->beamWaist;  
   elliptic=L->elliptic;  
   w=w0*sqrt(1.0+x*x/zR/zR);
   omega=2*M_PI*L->omega/D->omega;
   kx=2*M_PI*D->lambda/L->lambda;
   phi=atan(x/zR);

   if(t<2*rU)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2*rU && t<2*rU+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2*rU+flat && t<2*rU+flat+2*rD) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2*rU+2*rD+flat) 
      longitudinal=0.0;

   positionX=L->loadPointX;
   if(L->loadPointX>=D->minXSub && L->loadPointX<D->maxXSub) {
     laserOK=ON;
     positionX=L->loadPointX+istart;
   } else {
     laserOK=OFF;
   }

//   positionX=L->loadPointX+istart-D->minXSub;
//   if(positionX>=D->minXSub && positionX<D->maxXSub)
//     laserOK=1;
//   if(positionX>D->minXSub && positionX<=D->maxXSub &&
//      jC-D->minYSub>=jstart && jC-D->minYSub<jend &&
//      kC-D->minZSub>=kstart && kC-D->minZSub<kend)
//     laserOK=1;

   if(laserOK==ON && L->add==OFF)
   {
     if(L->polarity==2) {
       w2=w*w;
       w2p=w*w*elliptic*elliptic;
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         for(k=0; k<kend+3; k++) {
           z=(k-kstart+D->minZSub-kC)*D->dz;

           r2=y*y+z*z;
           pphi=x/zR*r2/w2-phi+kx*x-omega*t;
           amp=longitudinal*w0/w*exp(-y*y/w2)*exp(-z*z/w2p)*sin(pphi);
           D->Ey[positionX][j][k]=amp;            
         }
       }
     } else if(L->polarity==3) {
       w2=w*w;
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         for(k=0; k<kend+3; k++) {
           z=(k-kstart+D->minZSub-kC)*D->dz;

           r2=y*y+z*z;
           pphi=x/zR*r2/w2-phi+kx*x-omega*t;
           amp=longitudinal*w0/w*exp(-r2/w2)*sin(pphi);
           D->Ez[positionX][j][k]=amp;            
         }
       }
     }  
   } else if(laserOK==ON && L->add==ON)  {
     if(L->polarity==2) {
       w2=w*w;
       w2p=w*w*elliptic*elliptic;
/*
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         for(k=0; k<kend+3; k++) {
           z=(k-kstart+D->minZSub-kC)*D->dz;

           r2=y*y+z*z;
           if(r2==0) {
             pphi=x/zR*r2/w2-phi+kx*x-omega*t;
             amp=longitudinal*w0/w*exp(-r2/w2)*sin(pphi);
             D->Ey[positionX][j][k]+=amp;
           }
           else {
             cosP2=y*y/r2;
             sinP2=z*z/r2;

             zR=sqrt(zRY*zRY*cosP2+zRZ*zRZ*sinP2);
             phi=atan(x/zR);
             w=sqrt(wY*wY*cosP2+wZ*wZ*sinP2);
             w0=sqrt(w0Y*w0Y*cosP2+w0Z*w0Z*sinP2);
             w2=w*w;
             pphi=x/zR*r2/w2-phi+kx*x-omega*t;
             amp=longitudinal*w0/w*exp(-r2/w2)*sin(pphi);
             D->Ey[positionX][j][k]+=amp;            
           }        
         }
       }
*/
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         for(k=0; k<kend+3; k++) {
           z=(k-kstart+D->minZSub-kC)*D->dz;

           r2=y*y+z*z;
           pphi=x/zR*r2/w2-phi+kx*x-omega*t;
           amp=longitudinal*w0/w*exp(-y*y/w2)*exp(-z*z/w2p)*sin(pphi);
           D->Ey[positionX][j][k]+=amp;            
         }
       }
     } else if(L->polarity==3) {
       w2=w*w;
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         for(k=0; k<kend+3; k++) {
           z=(k-kstart+D->minZSub-kC)*D->dz;

           r2=y*y+z*z;
           pphi=x/zR*r2/w2-phi+kx*x-omega*t;
           amp=longitudinal*w0/w*exp(-r2/w2)*sin(pphi);
           D->Ez[positionX][j][k]+=amp;            
         }
       }
     }  
   }     //End of field is OK
}

void loadLaser_Gaussian_3D_Split(Domain *D,LaserList *L,double t)
{
   double rU,rD,longitudinal,t0,flat;
   double zR,w0,w,phi,omega,kx,pphi,amp,Phi;
   double x,y,z,r,r2,w2;
	double ***Right,***Left;
   int istart,iend,jstart,jend,kstart,kend;
   int positionX,rank,j,k,jC,kC,laserOK=OFF;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;

   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;
   flat=L->flat*D->divisionLambda*D->dt*L->lambda/D->lambda;

   jC=L->loadPointY;	//y position
   kC=L->loadPointZ;	//z position

   t0=2*rU;
   zR=L->rayleighLength;
   w0=L->beamWaist;  
   x=-L->focus;   
   w=w0*sqrt(1.0+x*x/zR/zR);
   phi=atan(x/zR);
   omega=2*M_PI*L->omega/D->omega;
   kx=2*M_PI*D->lambda/L->lambda;

   if(t<2*rU)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2*rU && t<2*rU+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2*rU+flat && t<2*rU+flat+2*rD) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2*rU+2*rD+flat) 
      longitudinal=0.0;

//   positionX=L->loadPointX+istart-D->minXSub;
   positionX=L->loadPointX+istart;
   if(positionX>=D->minXSub+istart && positionX<D->maxXSub+istart) laserOK=ON; else ;
//	printf("myrank=%d, positionX=%d, minXSub=%d, maxXSub=%d, laserOK=%d\n",myrank,positionX,D->minXSub,D->maxXSub,laserOK);

   w2=w*w;
   if(L->polarity==2) { Right=D->Pr; Left=D->Pl; }
   else if(L->polarity==3) { Right=D->Sr; Left=D->Sl; }
	else ;

   if(laserOK==ON && L->add==OFF)
   {
     for(j=0; j<jend+3; j++)  {
       y=(j-jstart+D->minYSub-jC)*D->dy;
       for(k=0; k<kend+3; k++)  {
         z=(k-kstart+D->minZSub-kC)*D->dz;
         r2=y*y+z*z;
         pphi=x/zR*r2/w2-phi+kx*x-omega*t;
         amp=longitudinal*w0/w*exp(-r2/w2)*sin(pphi);
         Right[positionX][j][k]=amp;            
         Left[positionX][j][k]=amp;           
       } 
     }
	} else if(laserOK==ON && L->add==ON)
   {
     for(j=0; j<jend+3; j++)  {
       y=(j-jstart+D->minYSub-jC)*D->dy;
       for(k=0; k<kend+3; k++)  {
         z=(k-kstart+D->minZSub-kC)*D->dz;
         r2=y*y+z*z;
         pphi=x/zR*r2/w2-phi+kx*x-omega*t;
         amp=longitudinal*w0/w*exp(-r2/w2)*sin(pphi);
         Right[positionX][j][k]+=amp;            
         Left[positionX][j][k]+=amp;           
       } 
     }
   } else ;    //end of laserOK
}

void loadLaser_SSTF_3D_Yee_Pukhov(Domain *D,LaserList *L,double t)
{
   double t0,omega,dy,dz,y,z,alphaY,alphaZ,minY,minZ,f12Y,f12Z,tcoefY,tcoefZ,zR,w0y,w0z,wy,wz;
   double sigW,minT,minTY,minTZ,sigYSSTF,sigZSSTF,f0,f,x,k0,a0,lambda0,phase,phi,pphiY,pphiZ,gaussX;
   double complex a_k0Y,coefY,arg1Y,sstfCompAY,mY,nY,tmpY,gaussTmpY,ampY,sstfCompA;
   double complex a_k0Z,coefZ,arg1Z,sstfCompAZ,mZ,nZ,tmpZ,gaussTmpZ,ampZ,gaussAmp;
	double ***field;
   int istart,iend,jstart,jend,kstart,kend,minYSub,minZSub;
   int positionX,rank,i,j,k,jC,kC,laserOK=OFF;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;

   dy=D->dy;   dz=D->dz;
   minY=D->minYDomain*dy*D->lambda;
   minZ=D->minZDomain*dz*D->lambda;
	minYSub=D->minYSub;
	minZSub=D->minZSub;

   jC=L->loadPointY;	//y position
   kC=L->loadPointZ;	//y position

   sigW=L->sigOmegaSSTF;  
   sigYSSTF=L->sigYSSTF;  
   alphaY=L->alphaYSSTF;
   sigZSSTF=L->sigZSSTF;  
   alphaZ=L->alphaZSSTF;

   minTY=-4.0/sigW*sqrt(1.0+alphaY*alphaY*sigW*sigW/sigYSSTF/sigYSSTF)-L->retard*D->dt*D->divisionLambda*2.0*M_PI/D->omega;
   minTZ=-4.0/sigW*sqrt(1.0+alphaZ*alphaZ*sigW*sigW/sigZSSTF/sigZSSTF)-L->retard*D->dt*D->divisionLambda*2.0*M_PI/D->omega;
	if(minTY>=minTZ) minT=minTY;
	else minT=minTZ;
	if(alphaY==0.0 && alphaZ!=0.0) { tcoefY=1.0; tcoefZ=0.0; }
	else if(alphaY!=0.0 && alphaZ==0.0) { tcoefY=0.0; tcoefZ=1.0; }
	else ;

   f0=L->lensFocus;
   f=L->focus*D->lambda;
   x=f0-f;
   k0=2*M_PI/L->lambda;
   a0=L->amplitude;
   lambda0=D->lambda;
   omega=velocityC*k0;

   f12Y=f0*f0*k0*k0*sigYSSTF*sigYSSTF*sigYSSTF*sigYSSTF/(4.0*f0*f0+k0*k0*sigYSSTF*sigYSSTF*sigYSSTF*sigYSSTF);
   f12Z=f0*f0*k0*k0*sigZSSTF*sigZSSTF*sigZSSTF*sigZSSTF/(4.0*f0*f0+k0*k0*sigZSSTF*sigZSSTF*sigZSSTF*sigZSSTF);

   // for normalization at the focus
   a_k0Y=f12Y/k0/k0/sigYSSTF/sigYSSTF-I*0.5/k0*(f0-f12Y/f0);
   a_k0Z=f12Z/k0/k0/sigZSSTF/sigZSSTF-I*0.5/k0*(f0-f12Z/f0);
   gaussTmpY=sigW*sigYSSTF*0.5*csqrt(M_PI/a_k0Y/(1-I*k0*sigYSSTF*sigYSSTF*0.5/f0));
   gaussTmpZ=sigW*sigZSSTF*0.5*csqrt(M_PI/a_k0Z/(1-I*k0*sigZSSTF*sigZSSTF*0.5/f0));

   a_k0Y=f12Y/k0/k0/sigYSSTF/sigYSSTF-I*0.5/k0*(x-f12Y/f0);
   a_k0Z=f12Z/k0/k0/sigZSSTF/sigZSSTF-I*0.5/k0*(x-f12Z/f0);
   mY=1.0+alphaY*alphaY*sigW*sigW*(x-f0)*(x-f0)/4.0/f0/f0/a_k0Y-I*k0*alphaY*alphaY*sigW*sigW*(x-f0)/2.0/f0/f0;
   mZ=1.0+alphaZ*alphaZ*sigW*sigW*(x-f0)*(x-f0)/4.0/f0/f0/a_k0Z-I*k0*alphaZ*alphaZ*sigW*sigW*(x-f0)/2.0/f0/f0;
   nY=k0*alphaY*sigW/f0+I*alphaY*sigW*(x-f0)/2.0/f0/a_k0Y;
   nZ=k0*alphaZ*sigW/f0+I*alphaZ*sigW*(x-f0)/2.0/f0/a_k0Z;
   tmpY=M_PI/mY/a_k0Y/(1-I*k0*sigYSSTF*sigYSSTF*0.5/f0);
   tmpZ=M_PI/mZ/a_k0Z/(1-I*k0*sigZSSTF*sigZSSTF*0.5/f0);
   tmpY=csqrt(tmpY)*sigYSSTF*sigW*0.5;
   tmpZ=csqrt(tmpZ)*sigZSSTF*sigW*0.5;
   coefY=tmpY/gaussTmpY;			//normalization factor at z=x
   coefZ=tmpZ/gaussTmpZ;			//normalization factor at z=x

   // Gaussian property
   w0y=L->sigY_t;
   w0z=L->sigZ_t;
   zR=M_PI/L->lambda*(w0y*w0y*tcoefY+w0z*w0z*tcoefZ);
   gaussX=(L->focus-L->loadPointX*D->dx)*L->lambda;
   wy=w0y*sqrt(1.0+gaussX*gaussX/zR/zR);
   wz=w0z*sqrt(1.0+gaussX*gaussX/zR/zR);
   phi=atan(gaussX/zR);      
 

   //time
   t0=minT+t*2.0*M_PI/D->omega;

   positionX=L->loadPointX+istart;
   if(positionX>=D->minXSub+istart && positionX<D->maxXSub+istart) laserOK=ON; else ;

   if(L->polarity==2) { field=D->Ey; }
   else if(L->polarity==3) { field=D->Ez; }
	else ;
//printf("myrank=%d, positionX=%d, laserOK=%d\n",myrank,positionX,laserOK);

   if(laserOK==ON && L->add==OFF)  
	{
     for(j=0; j<jend+3; j++) {
       y=(j-jstart+minYSub-jC)*dy*lambda0;
       for(k=0; k<kend+3; k++) {
         z=(k-kstart+minZSub-kC)*dz*lambda0;
			arg1Y=-0.25*y*y/a_k0Y-0.25/mY*(sigW*t0+nY*y)*(sigW*t0+nY*y);
         arg1Z=-0.25*z*z/a_k0Z-0.25/mZ*(sigW*t0+nZ*z)*(sigW*t0+nZ*z);
         sstfCompAY=coefY*cexp(arg1Y+I*omega*t0);
         sstfCompAZ=coefZ*cexp(arg1Z+I*omega*t0);
			sstfCompA=sstfCompAY*(1.0-tcoefY)+sstfCompAZ*(1.0-tcoefZ);
//         compA=coefY*coefZ*sqrt(w0/w)*cexp(arg1Y+arg1Z+I*omega*t0);
         pphiY=x/zR*y*y/wy/wy-0.5*phi;
         pphiZ=x/zR*z*z/wz/wz-0.5*phi;
         ampY=sqrt(w0y/wy)*exp(-y*y/wy/wy)*cexp(I*pphiY);
         ampZ=sqrt(w0z/wz)*exp(-z*z/wz/wz)*cexp(I*pphiZ);
         gaussAmp=ampY*tcoefY+ampZ*tcoefZ;
         field[positionX][j][k]=a0*creal(sstfCompA*gaussAmp); 
       }
     }
   }
   else if(laserOK==ON && L->add==ON)  
	{
     for(j=0; j<jend+3; j++) {
       y=(j-jstart+minYSub-jC)*dy*lambda0;
       for(k=0; k<kend+3; k++) {
         z=(k-kstart+minZSub-kC)*dz*lambda0;
			arg1Y=-0.25*y*y/a_k0Y-0.25/mY*(sigW*t0+nY*y)*(sigW*t0+nY*y);
         arg1Z=-0.25*z*z/a_k0Z-0.25/mZ*(sigW*t0+nZ*z)*(sigW*t0+nZ*z);
         sstfCompAY=coefY*cexp(arg1Y+I*omega*t0);
         sstfCompAZ=coefZ*cexp(arg1Z+I*omega*t0);
			sstfCompA=sstfCompAY*(1.0-tcoefY)+sstfCompAZ*(1.0-tcoefZ);

         pphiY=x/zR*y*y/wy/wy-0.5*phi;
         pphiZ=x/zR*z*z/wz/wz-0.5*phi;
         ampY=sqrt(w0y/wy)*exp(-y*y/wy/wy)*cexp(I*pphiY);
         ampZ=sqrt(w0z/wz)*exp(-z*z/wz/wz)*cexp(I*pphiZ);
         gaussAmp=ampY*tcoefY+ampZ*tcoefZ;
         field[positionX][j][k]+=a0*creal(sstfCompA*gaussAmp); 
       }
     }
   }
   else ;   //End of field is OK
}

void loadLaser_SSTF_3D_Split(Domain *D,LaserList *L,double t)
{
   double t0,omega,dy,dz,y,z,alphaY,alphaZ,minY,minZ,f12Y,f12Z,tcoefY,tcoefZ,pphiY,pphiZ;
   double sigW,minT,minTY,minTZ,sigYSSTF,sigZSSTF,f0,f,x,k0,a0,lambda0,phase,w0y,w0z,wy,wz,zR,phi,gaussX;
   double complex a_k0Y,coefY,arg1Y,sstfCompAY,mY,nY,tmpY,gaussTmpY,ampY,sstfCompA;
   double complex a_k0Z,coefZ,arg1Z,sstfCompAZ,mZ,nZ,tmpZ,gaussTmpZ,ampZ,gaussAmp;
	double ***Right,***Left;
   int istart,iend,jstart,jend,kstart,kend,minYSub,minZSub;
   int positionX,rank,i,j,k,jC,kC,laserOK=OFF;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;

   dy=D->dy;   dz=D->dz;
   minY=D->minYDomain*dy*D->lambda;
   minZ=D->minZDomain*dz*D->lambda;
	minYSub=D->minYSub;
	minZSub=D->minZSub;

   jC=L->loadPointY;	//y position
   kC=L->loadPointZ;	//y position

   sigW=L->sigOmegaSSTF;  
   sigYSSTF=L->sigYSSTF;  
   alphaY=L->alphaYSSTF;
   sigZSSTF=L->sigZSSTF;  
   alphaZ=L->alphaZSSTF;

   minTY=-4.0/sigW*sqrt((1.0+alphaY*alphaY*sigW*sigW/sigYSSTF/sigYSSTF))-L->retard*D->dt*D->divisionLambda*2.0*M_PI/D->omega;
   minTZ=-4.0/sigW*sqrt((1.0+alphaZ*alphaZ*sigW*sigW/sigZSSTF/sigZSSTF))-L->retard*D->dt*D->divisionLambda*2.0*M_PI/D->omega;
	if(minTY>=minTZ) minT=minTY;
	else minT=minTZ;
	if(alphaY==0.0 && alphaZ!=0.0) { tcoefY=1.0; tcoefZ=0.0; }
	else if(alphaY!=0.0 && alphaZ==0.0) { tcoefY=0.0; tcoefZ=1.0; }
	else  ;

   f0=L->lensFocus;
   f=L->focus*D->lambda;
   x=f0-f;
   k0=2*M_PI/L->lambda;
   a0=L->amplitude;
   lambda0=D->lambda;
   omega=velocityC*k0;

   f12Y=f0*f0*k0*k0*sigYSSTF*sigYSSTF*sigYSSTF*sigYSSTF/(4.0*f0*f0+k0*k0*sigYSSTF*sigYSSTF*sigYSSTF*sigYSSTF);
   f12Z=f0*f0*k0*k0*sigZSSTF*sigZSSTF*sigZSSTF*sigZSSTF/(4.0*f0*f0+k0*k0*sigZSSTF*sigZSSTF*sigZSSTF*sigZSSTF);

   // for normalization at the focus
   a_k0Y=f12Y/k0/k0/sigYSSTF/sigYSSTF-I*0.5/k0*(f0-f12Y/f0);
   a_k0Z=f12Z/k0/k0/sigZSSTF/sigZSSTF-I*0.5/k0*(f0-f12Z/f0);
   gaussTmpY=sigW*sigYSSTF*0.5*csqrt(M_PI/a_k0Y/(1-I*k0*sigYSSTF*sigYSSTF*0.5/f0));
   gaussTmpZ=sigW*sigZSSTF*0.5*csqrt(M_PI/a_k0Z/(1-I*k0*sigZSSTF*sigZSSTF*0.5/f0));

   a_k0Y=f12Y/k0/k0/sigYSSTF/sigYSSTF-I*0.5/k0*(x-f12Y/f0);
   a_k0Z=f12Z/k0/k0/sigZSSTF/sigZSSTF-I*0.5/k0*(x-f12Z/f0);
   mY=1.0+alphaY*alphaY*sigW*sigW*(x-f0)*(x-f0)/4.0/f0/f0/a_k0Y-I*k0*alphaY*alphaY*sigW*sigW*(x-f0)/2.0/f0/f0;
   mZ=1.0+alphaZ*alphaZ*sigW*sigW*(x-f0)*(x-f0)/4.0/f0/f0/a_k0Z-I*k0*alphaZ*alphaZ*sigW*sigW*(x-f0)/2.0/f0/f0;
   nY=k0*alphaY*sigW/f0+I*alphaY*sigW*(x-f0)/2.0/f0/a_k0Y;
   nZ=k0*alphaZ*sigW/f0+I*alphaZ*sigW*(x-f0)/2.0/f0/a_k0Z;
   tmpY=M_PI/mY/a_k0Y/(1-I*k0*sigYSSTF*sigYSSTF*0.5/f0);
   tmpZ=M_PI/mZ/a_k0Z/(1-I*k0*sigZSSTF*sigZSSTF*0.5/f0);
   tmpY=csqrt(tmpY)*sigYSSTF*sigW*0.5;
   tmpZ=csqrt(tmpZ)*sigZSSTF*sigW*0.5;
   coefY=tmpY/gaussTmpY;			//normalization factor at z=x
   coefZ=tmpZ/gaussTmpZ;			//normalization factor at z=x

   // Gaussian property 
   w0y=L->sigY_t;
   w0z=L->sigZ_t;
   zR=M_PI/L->lambda*(w0y*w0y*tcoefY+w0z*w0z*tcoefZ);
   gaussX=-(L->loadPointX*D->dx-L->focus)*D->lambda;
   wy=w0y*sqrt(1.0+gaussX*gaussX/zR/zR);
   wz=w0z*sqrt(1.0+gaussX*gaussX/zR/zR);
   phi=atan(gaussX/zR);      
 

   //time
   t0=minT+t*2.0*M_PI/D->omega;

   positionX=L->loadPointX+istart;
   if(positionX>=D->minXSub+istart && positionX<D->maxXSub+istart) laserOK=ON; else ;

   if(L->polarity==2) { Right=D->Pr; Left=D->Pl; }
   else if(L->polarity==3) { Right=D->Sr; Left=D->Sl; }
	else ;
//printf("myrank=%d, positionX=%d, laserOK=%d\n",myrank,positionX,laserOK);

   if(laserOK==ON && L->add==OFF)  
	{
     for(j=0; j<jend+3; j++) {
       y=(j-jstart+minYSub-jC)*dy*lambda0;
       for(k=0; k<kend+3; k++) {
         z=(k-kstart+minZSub-kC)*dz*lambda0;
			arg1Y=-0.25*y*y/a_k0Y-0.25/mY*(sigW*t0+nY*y)*(sigW*t0+nY*y);
         arg1Z=-0.25*z*z/a_k0Z-0.25/mZ*(sigW*t0+nZ*z)*(sigW*t0+nZ*z);
         sstfCompAY=coefY*cexp(arg1Y+I*omega*t0);
         sstfCompAZ=coefZ*cexp(arg1Z+I*omega*t0);
			sstfCompA=sstfCompAY*(1.0-tcoefY)+sstfCompAZ*(1.0-tcoefZ);
//         compA=coefY*coefZ*sqrt(w0/w)*cexp(arg1Y+arg1Z+I*omega*t0);
         pphiY=gaussX/zR*y*y/wy/wy-0.5*phi;
         pphiZ=gaussX/zR*z*z/wz/wz-0.5*phi;
         ampY=sqrt(w0y/wy)*exp(-y*y/wy/wy)*cexp(I*pphiY);
         ampZ=sqrt(w0z/wz)*exp(-z*z/wz/wz)*cexp(I*pphiZ);
         gaussAmp=ampY*tcoefY+ampZ*tcoefZ;
         Right[positionX][j][k]=a0*creal(sstfCompA*gaussAmp); 
         Left[positionX][j][k]=a0*creal(sstfCompA*gaussAmp); 
       }
     }
   }
   else if(laserOK==ON && L->add==ON)  
	{
     for(j=0; j<jend+3; j++) {
       y=(j-jstart+minYSub-jC)*dy*lambda0;
       for(k=0; k<kend+3; k++) {
         z=(k-kstart+minZSub-kC)*dz*lambda0;
			arg1Y=-0.25*y*y/a_k0Y-0.25/mY*(sigW*t0+nY*y)*(sigW*t0+nY*y);
         arg1Z=-0.25*z*z/a_k0Z-0.25/mZ*(sigW*t0+nZ*z)*(sigW*t0+nZ*z);
         sstfCompAY=coefY*cexp(arg1Y+I*omega*t0);
         sstfCompAZ=coefZ*cexp(arg1Z+I*omega*t0);
			sstfCompA=sstfCompAY*(1.0-tcoefY)+sstfCompAZ*(1.0-tcoefZ);

         pphiY=gaussX/zR*y*y/wy/wy-0.5*phi;
         pphiZ=gaussX/zR*z*z/wz/wz-0.5*phi;
         ampY=sqrt(w0y/wy)*exp(-y*y/wy/wy)*cexp(I*pphiY);
         ampZ=sqrt(w0z/wz)*exp(-z*z/wz/wz)*cexp(I*pphiZ);
         gaussAmp=ampY*tcoefY+ampZ*tcoefZ;			
         Right[positionX][j][k]+=a0*creal(sstfCompA*gaussAmp); 
         Left[positionX][j][k]+=a0*creal(sstfCompA*gaussAmp); 
       }
     }
   }
   else ;   //End of field is OK
}

/*
void boostLoadLaser2D(Domain *D,LaserList *L)
{
   double x,x0,rU,rD,beta,gamma,labWaist,f,yy;
   double zR,w0,z,w,phi,omega,k,y,pphi,amp,unitCompen;
   int i,j,jC;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;

   jC=(int)(D->ny*0.5);	//center position
   gamma=D->gamma;
   beta=D->beta;

   labWaist=L->beamWaist*D->lambda;
   zR=L->rayleighLength/gamma;
//   zR=pi/(L->lambda/D->gamma/(1.0+D->beta))*labWaist*labWaist/gamma/D->lambda;
   w0=labWaist/D->lambda;  
   f=-L->focus/gamma;   
   omega=k=2*pi;
   x0=-2*rU;
   unitCompen=gamma*gamma*(1.0-beta*beta);

   if(D->fieldType==1)
   {
     FieldDSX **field;
     field=D->fieldDSX;

     if(L->polarity==2)
     {
       for(i=2; i<D->nxSub+2; i++)
         for(j=2; j<D->nySub+2; j++)
         {
           x=(i+D->minXSub-2)*D->dx;
           y=(j-2+D->minYSub-jC)*D->dy;
           z=(f+x+x0);          
           w=w0*sqrt(1.0+z*z/zR/zR);
           if(x>=2*x0)
           {
             phi=atan(z/zR);
             pphi=z/zR*y*y/w/w-0.5*phi+unitCompen*k*z;
             amp=L->amplitude*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi)*exp(-(x-x0)*(x-x0)/rU/rU);
//             amp=L->amplitude*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi);
             field[i][j].Pr=amp;            
//             field[i][j].Pl=amp;            
           }
           else      {
             field[i][j].Pr=0.0;            
//             field[i][j].Pl=0.0;            
           }
         }
     }
     else if(L->polarity==3)
     {
       for(i=1; i<=D->nxSub; i++)
         for(j=1; j<=D->nySub; j++)
         {
           x=(i+D->minXSub-1)*D->dx;
           z=f+x+x0;          
           if(x>=2*x0)
           {
             y=(j-1+D->minYSub-jC)*D->dy;
             w=w0*sqrt(1.0+z*z/zR/zR);
             phi=atan(z/zR);
             pphi=z/zR*y*y/w/w-0.5*phi+unitCompen*k*z;
             amp=L->amplitude*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi)*exp(-(x-x0)*(x-x0)/rU/rU);
             field[i][j].Sr=amp;            
             field[i][j].Sl=amp;            
           }
           else      {
             field[i][j].Sr=0.0;            
             field[i][j].Sl=0.0;            
           }
         }
     }
   }     //End of fieldType=1
}
*/

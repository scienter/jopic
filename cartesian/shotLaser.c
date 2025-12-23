#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <complex.h>
#include "mpi.h"

void shotLaser_Gaussian_1D(Domain *D,LaserList *L);
void shotLaser_Gaussian_2D(Domain *D,LaserList *L);
void shotLaser_SSTF_2D(Domain *D,LaserList *L);
void shotLaser_Gaussian_3D(Domain *D,LaserList *L);
void shotLaser_SSTF_3D(Domain *D,LaserList *L);


void shotLaser(Domain *D,LaserList *L)
{
    switch((D->fieldType-1)*3+D->dimension)  {
    case (Split-1)*3+1 :
      if(L->mode==Gaussian) shotLaser_Gaussian_1D(D,L);
      else ;
    case (Yee-1)*3+1 :
    case (Pukhov-1)*3+1 :
      break;
    case (Yee-1)*3+2 :
    case (Pukhov-1)*3+2 :
    case (Split-1)*3+2 :
      if(L->mode==Gaussian) shotLaser_Gaussian_2D(D,L);
      else if(L->mode==SSTF) shotLaser_SSTF_2D(D,L);
      else ;
      break;
    case (Yee-1)*3+3 :
    case (Pukhov-1)*3+3 :
    case (Split-1)*3+3 :
      if(L->mode==Gaussian) shotLaser_Gaussian_3D(D,L);
      else if(L->mode==SSTF) shotLaser_SSTF_3D(D,L);
      else ;
      break;
    default :
      printf("In shotLaser, what is field_type? and what is dimension?\n");
    }
}

void shotLaser_Gaussian_1D(Domain *D,LaserList *L)
{
   double rU,rD,t0,flat,x1,x2,factR,factL;
   double zR,w0,w,phi,omega,kx,pphi,amp,focus,a0;
   double x,y,z,r2,w2,retard,positionX,dx,dy,dt,dtOverdy;
   double ***field1,***field2;	
   int istart,iend,jstart,jend,kstart,kend,minj,maxj;
   int rank,i,j,k,jC,kC,minXSub,minYSub;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;  j=0;
   kstart=D->kstart;   kend=D->kend;  k=0;
   minXSub=D->minXSub; minYSub=D->minYSub;
   dx=D->dx;           dy=D->dy; dt=D->dt;

   rU=L->rU*L->lambda/D->lambda;  //D->lambda;
   rD=L->rD*L->lambda/D->lambda; //*D->lambda;
   flat=L->flat*L->lambda/D->lambda;
   retard=L->retard*L->lambda/D->lambda;

   a0=L->amplitude;
   focus=L->focus;
   zR=L->rayleighLength;	//normalized
   w0=L->beamWaist;  		//normalized
   positionX=L->loadPointX*D->dx;
   x1=positionX-flat*0.5-retard;
   x2=positionX+flat*0.5-retard;
   kx=2.0*M_PI/L->lambda*D->lambda;

   factR=1.0; factL=1.0;
   if(L->polarity==2) {
     switch(D->fieldType)  {
     case Split:
       field1=D->Pr; field2=D->Pl;
       if(L->direction==RIGHT) { factR=1.0; factL=0.0; }
       else if(L->direction==LEFT) { factR=0.0; factL=1.0; }
       else if(L->direction==BOTH) { factR=1.0; factL=1.0; }
       else { factR=0.0; factL=0.0; }       
       break;
     case Yee:
     case Pukhov:
       field1=D->Ey; field2=D->Bz;
       break;
     }
   } else if(L->polarity==3) {
     switch(D->fieldType)  {
     case Split:
       field1=D->Sr; field2=D->Sl;
       if(L->direction==RIGHT) { factR=1.0; factL=0.0; }
       else if(L->direction==LEFT) { factR=0.0; factL=1.0; }
       else if(L->direction==BOTH) { factR=1.0; factL=1.0; }
       else { factR=0.0; factL=0.0; }
       break;
     case Yee:
     case Pukhov:
       field1=D->Ez; field2=D->By;
       break;
     }
   } else ;

   switch (L->add)  {
   case OFF :
     for(i=1; i<iend+2; i++)  {
       x=(minXSub+i)*dx+retard;
       if(x<=x1)      amp=exp(-(x-positionX)*(x-positionX)/rU/rU)*a0*sin(kx*x);
       else if(x>x1 && x<=x2) amp=a0*sin(kx*x);
       else if(x>x2)  amp=exp(-(x-positionX)*(x-positionX)/rD/rD)*a0*sin(kx*x);
       field1[i][j][k]=amp*factR;
       field2[i][j][k]=amp*factL;
     }
     break;
   case ON :
     for(i=1; i<iend+3; i++)  {
       x=(minXSub+i)*dx+retard;
       if(x<=x1)      amp=exp(-(x-positionX)*(x-positionX)/rU/rU)*a0*sin(kx*x);
       else if(x>x1 && x<=x2) amp=a0*sin(kx*x);
       else if(x>x2)  amp=exp(-(x-positionX)*(x-positionX)/rD/rD)*a0*sin(kx*x);
       field1[i][j][k]+=amp*factR; 
       field2[i][j][k]+=amp*factL;
     }
     break; 
   default :
     printf("In shotLaser.c, what laser mode?\n");
     break;
   }		//End of switch
}


void shotLaser_Gaussian_2D(Domain *D,LaserList *L)
{
   double rU,rD,longitudinal,t0,flat,x1,x2;
   double zR,w0,w,phi,omega,kx,pphi,amp,focus,a0;
   double x,y,z,r2,w2,retard,positionX,dx,dy,dt,dtOverdy;
   double ***field1,***field2;	
   int istart,iend,jstart,jend,kstart,kend,minj,maxj;
   int rank,i,j,k,jC,kC,minXSub,minYSub;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;  k=0;
   minXSub=D->minXSub; minYSub=D->minYSub;
   dx=D->dx;           dy=D->dy; dt=D->dt;
   dtOverdy=D->dt/D->dy;

   rU=L->rU*L->lambda/D->lambda;  //D->lambda;
   rD=L->rD*L->lambda/D->lambda; //*D->lambda;
   flat=L->flat*L->lambda/D->lambda;
   retard=L->retard*L->lambda/D->lambda;

   jC=L->loadPointY;	//y position

   a0=L->amplitude;
   focus=L->focus;
   zR=L->rayleighLength;	//normalized
   w0=L->beamWaist;  		//normalized
   positionX=L->loadPointX*D->dx;
   x1=positionX-flat*0.5-retard;
   x2=positionX+flat*0.5-retard;
   kx=2.0*M_PI/L->lambda*D->lambda;


   if(L->polarity==2) {
     switch(D->fieldType)  {
     case Split:
       field1=D->Pr; field2=D->Pl;
       break;
     case Yee:
     case Pukhov:
       field1=D->Ey; field2=D->Bz;
       break;
	  }
	} else if(L->polarity==3) {
     switch(D->fieldType)  {
     case Split:
       field1=D->Sr; field2=D->Sl;
       break;
     case Yee:
     case Pukhov:
       field1=D->Ez; field2=D->By;
       break;
	  }
	} else ;

   switch (L->add)  {
   case OFF :
     for(i=1; i<iend+2; i++)  {
       x=(minXSub+i)*dx+retard;
       w=w0*sqrt(1.0+(x-focus)*(x-focus)/zR/zR);
       w2=w*w;
       if(x<=x1)      longitudinal=exp(-(x-positionX)*(x-positionX)/rU/rU)*a0;
       else if(x>x1 && x<=x2) longitudinal=a0;
       else if(x>x2)  longitudinal=exp(-(x-positionX)*(x-positionX)/rD/rD)*a0;
       phi=atan((x-focus)/zR);
       for(j=0; j<jend+3; j++)         {
         y=(j-jstart+minYSub-jC)*dy;
         r2=y*y;
         pphi=(x-focus)/zR*r2/w2-0.5*phi+kx*x;           
         amp=longitudinal*sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
         field1[i][j][k]=amp;            
         field2[i][j][k]=amp;            
       }
     }
     break;
   case ON :
     for(i=1; i<iend+3; i++)  {
       x=(minXSub+i)*dx+retard;
       w=w0*sqrt(1.0+(x-focus)*(x-focus)/zR/zR);
       w2=w*w;
       if(x<=x1)      longitudinal=exp(-(x-positionX)*(x-positionX)/rU/rU)*a0;
       else if(x>x1 && x<=x2) longitudinal=a0;
       else if(x>x2)  longitudinal=exp(-(x-positionX)*(x-positionX)/rD/rD)*a0;
       phi=atan((x-focus)/zR);
       for(j=0; j<jend+3; j++)         {
         y=(j-jstart+minYSub-jC)*dy;
         r2=y*y;
         pphi=(x-focus)/zR*r2/w2-0.5*phi+kx*x;           
         amp=longitudinal*sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
         field1[i][j][k]+=amp; 
         field2[i][j][k]+=amp;
       }
     }
     break; 
   default :
     printf("In shotLaser.c, what laser mode?\n");
     break;
   }		//End of switch
}

void shotLaser_SSTF_2D(Domain *D,LaserList *L)
{
   double dx,dy,x,y,alpha,minY,minXSub,minYSub,f12,positionX;
   double sigW,minT,sigY,f0,f,k0,a0,lambda0;
   double complex a_k0,amp,arg1,argY,compA,m,n,tmp,gaussTmp;
   double ***field1,***field2;	
   int istart,iend,jstart,jend,kstart,kend;
   int i,j,k,jC,kC;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;

   dx=D->dx;   dy=D->dy;
   minY=D->minYDomain*dy*D->lambda;
   minXSub=D->minXSub;   minYSub=D->minYSub;

   lambda0=D->lambda;
   jC=L->loadPointY;	//y position
   positionX=L->loadPointX*dx*lambda0;

   sigW=L->sigOmegaSSTF;  
   sigY=L->sigYSSTF;  
   alpha=L->alphaYSSTF;

   f0=L->lensFocus;
   f=L->focus*D->lambda;
   k0=2*M_PI/L->lambda;
   a0=L->amplitude;

   f12=f0*f0*k0*k0*sigY*sigY*sigY*sigY/(4.0*f0*f0+k0*k0*sigY*sigY*sigY*sigY);

   // for normalization at the focus
   a_k0=f12/k0/k0/sigY/sigY-I*0.5/k0*(f0-f12/f0);
   gaussTmp=sigW*sigY*0.5*csqrt(M_PI/a_k0/(1-I*k0*sigY*sigY*0.5/f0));


   if(L->polarity==2) {
     switch(D->fieldType)  {
     case Split:
       field1=D->Pr; field2=D->Pl;
       break;
     case Yee:
     case Pukhov:
       field1=D->Ey; field2=D->Bz;
       break;
	  }
	} else if(L->polarity==3) {
     switch(D->fieldType)  {
     case Split:
       field1=D->Sr; field2=D->Sl;
       break;
     case Yee:
     case Pukhov:
       field1=D->Ez; field2=D->By;
       break;
	  }
	} else ;

   switch (L->add)  {
   case OFF :
     for(i=0; i<iend+3; i++)  {
       x=(minXSub+i)*dx*lambda0+(f0-f);
       a_k0=f12/k0/k0/sigY/sigY-I*0.5/k0*(x-f12/f0);
       m=1.0+alpha*alpha*sigW*sigW*(x-f0)*(x-f0)/4.0/f0/f0/a_k0-I*k0*alpha*alpha*sigW*sigW*(x-f0)/2.0/f0/f0;
       n=k0*alpha*sigW/f0+I*alpha*sigW*(x-f0)/2.0/f0/a_k0;
       tmp=M_PI/m/a_k0/(1-I*k0*sigY*sigY*0.5/f0);
       tmp=csqrt(tmp)*sigY*sigW*0.5;
       amp=tmp/gaussTmp*a0;			//normalization factor at z=x
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+minYSub-jC)*dy*lambda0;
	      argY=sigW/velocityC*(x-positionX-(f0-f))+n*y;
         arg1=-0.25*y*y/a_k0-0.25/m*argY*argY-I*k0*(x-positionX-(f0-f));
         compA=amp*cexp(arg1);
         field1[i][j][0]=creal(compA); 
         field2[i][j][0]=creal(compA); 
		 }
     }
	  break;
   case ON :
     for(i=0; i<iend+3; i++)  {
       x=(minXSub+i)*dx*lambda0+(f0-f);
       a_k0=f12/k0/k0/sigY/sigY-I*0.5/k0*(x-f12/f0);
       m=1.0+alpha*alpha*sigW*sigW*(x-f0)*(x-f0)/4.0/f0/f0/a_k0-I*k0*alpha*alpha*sigW*sigW*(x-f0)/2.0/f0/f0;
       n=k0*alpha*sigW/f0+I*alpha*sigW*(x-f0)/2.0/f0/a_k0;
       tmp=M_PI/m/a_k0/(1-I*k0*sigY*sigY*0.5/f0);
       tmp=csqrt(tmp)*sigY*sigW*0.5;
       amp=tmp/gaussTmp*a0;			//normalization factor at z=x
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+minYSub-jC)*dy*lambda0;
	      argY=sigW/velocityC*(x-positionX-(f0-f))+n*y;
         arg1=-0.25*y*y/a_k0-0.25/m*argY*argY-I*k0*(x-positionX-(f0-f));
         compA=amp*cexp(arg1);
         field1[i][j][0]+=creal(compA); 
         field2[i][j][0]+=creal(compA); 
		 }
     }
	  break;
   } 
}


void shotLaser_Gaussian_3D(Domain *D,LaserList *L)
{
   double rU,rD,longitudinal,t0,flat,x1,x2,elliptic;
   double zR,w0,w,phi,omega,kx,pphi,amp,focus,a0;
   double x,y,z,r2,w2,w2p,retard,positionX,dx,dy,dz,dt;
   double ***field1,***field2;	
   int istart,iend,jstart,jend,kstart,kend,minj,maxj;
   int rank,i,j,k,jC,kC,minXSub,minYSub,minZSub;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   minXSub=D->minXSub; 
	minYSub=D->minYSub;
   minZSub=D->minZSub;
   dx=D->dx; dy=D->dy; dz=D->dz; dt=D->dt;

   rU=L->rU*L->lambda/D->lambda;  //D->lambda;
   rD=L->rD*L->lambda/D->lambda; //*D->lambda;
   flat=L->flat*L->lambda/D->lambda;
   retard=L->retard*L->lambda/D->lambda;

   jC=L->loadPointY;	//y position
   kC=L->loadPointZ;	//y position

   a0=L->amplitude;
   focus=L->focus;
   zR=L->rayleighLength;	//normalized
   w0=L->beamWaist;  		//normalized
	elliptic=L->elliptic;
   positionX=L->loadPointX*D->dx;
   x1=positionX-flat*0.5-retard;
   x2=positionX+flat*0.5-retard;
   kx=2.0*M_PI/L->lambda*D->lambda;

   if(L->polarity==2) {
     switch(D->fieldType)  {
     case Split:
       field1=D->Pr; field2=D->Pl;
       break;
     case Yee:
     case Pukhov:
       field1=D->Ey; field2=D->Bz;
       break;
	  }
	} else if(L->polarity==3) {
     switch(D->fieldType)  {
     case Split:
       field1=D->Sr; field2=D->Sl;
       break;
     case Yee:
     case Pukhov:
       field1=D->Ez; field2=D->By;
       break;
	  }
	} else ;

   switch (L->add)  {
   case OFF :
     for(i=0; i<iend+3; i++)  {
       x=(minXSub+i)*dx+retard;
       w=w0*sqrt(1.0+(x-focus)*(x-focus)/zR/zR);
       w2=w*w;
       w2p=w*w*elliptic*elliptic;
       if(x<=x1)      longitudinal=exp(-(x-positionX)*(x-positionX)/rU/rU)*a0;
       else if(x>x1 && x<=x2) longitudinal=a0;
       else if(x>x2)  longitudinal=exp(-(x-positionX)*(x-positionX)/rD/rD)*a0;
       phi=atan((x-focus)/zR);
       for(j=0; j<jend+3; j++)         {
         y=(j-jstart+minYSub-jC)*dy;
         for(k=0; k<kend+3; k++)         {
           z=(k-kstart+minZSub-kC)*dz;
           r2=y*y+z*z;
           pphi=(x-focus)/zR*r2/w2-0.5*phi+kx*x;           
           amp=longitudinal*w0/w*exp(-y*y/w2)*exp(-z*z/w2p)*sin(pphi);
           field1[i][j][k]=amp;            
           field2[i][j][k]=amp;           
         }
		 }
     }
     break;
   case ON :
     for(i=0; i<iend+3; i++)  {
       x=(minXSub+i)*dx+retard;
       w=w0*sqrt(1.0+(x-focus)*(x-focus)/zR/zR);
       w2=w*w;
       w2p=w*w*elliptic*elliptic;
       if(x<=x1)      longitudinal=exp(-(x-positionX)*(x-positionX)/rU/rU)*a0;
       else if(x>x1 && x<=x2) longitudinal=a0;
       else if(x>x2)  longitudinal=exp(-(x-positionX)*(x-positionX)/rD/rD)*a0;
       phi=atan((x-focus)/zR);
       for(j=0; j<jend+3; j++)         {
         y=(j-jstart+minYSub-jC)*dy;
         for(k=0; k<kend+3; k++)         {
           z=(k-kstart+minZSub-kC)*dz;
           r2=y*y+z*z;
           pphi=(x-focus)/zR*r2/w2-0.5*phi+kx*x;           
           amp=longitudinal*w0/w*exp(-y*y/w2)*exp(-z*z/w2p)*sin(pphi);
           field1[i][j][k]+=amp;            
           field2[i][j][k]+=amp;           
         }
		 }
     }
     break;
   default :
     printf("In shotLaser.c, what laser mode?\n");
     break;
   }		//End of switch
}

void shotLaser_SSTF_3D(Domain *D,LaserList *L)
{
   double dx,dy,dz,x,y,z,positionX;
   double sigW,minT,sigY,f0,f,k0,a0,lambda0;
	double alphaY,minY,sigYSSTF,f12Y,tcoefY,minYSub,minXSub;
	double alphaZ,minZ,sigZSSTF,f12Z,tcoefZ,minZSub;
   double complex coefY,a_k0Y,arg1Y,argY,mY,nY,tmpY,gaussTmpY,amp,compA;
   double complex coefZ,a_k0Z,arg1Z,argZ,mZ,nZ,tmpZ,gaussTmpZ;
   double ***field1,***field2;	
   int istart,iend,jstart,jend,kstart,kend;
   int i,j,k,jC,kC;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;

   dx=D->dx;   dy=D->dy;  dz=D->dz;
   minXSub=D->minXSub;   minYSub=D->minYSub; minZSub=D->minZSub;
   minY=D->minYDomain*dy*D->lambda;
   minZ=D->minZDomain*dz*D->lambda;

   jC=L->loadPointY;	//y position
   kC=L->loadPointZ;	//z position
   positionX=L->loadPointX*dx*D->lambda;

   sigW=L->sigOmegaSSTF;
	sigYSSTF=L->sigYSSTF;
	alphaY=L->alphaYSSTF;
	sigZSSTF=L->sigZSSTF;
	alphaZ=L->alphaZSSTF;

   if(L->polarity==2) {
     switch(D->fieldType)  {
     case Split:
       field1=D->Pr; field2=D->Pl;
       break;
     case Yee:
     case Pukhov:
       field1=D->Ey; field2=D->Bz;
       break;
	  }
	} else if(L->polarity==3) {
     switch(D->fieldType)  {
     case Split:
       field1=D->Sr; field2=D->Sl;
       break;
     case Yee:
     case Pukhov:
       field1=D->Ez; field2=D->By;
       break;
	  }
	} else ;


   f0=L->lensFocus;
   f=L->focus*D->lambda;
   k0=2*M_PI/L->lambda;
   a0=L->amplitude;
   lambda0=D->lambda;

   f12Y=f0*f0*k0*k0*sigYSSTF*sigYSSTF*sigYSSTF*sigYSSTF/(4.0*f0*f0+k0*k0*sigYSSTF*sigYSSTF*sigYSSTF*sigYSSTF);
   f12Z=f0*f0*k0*k0*sigZSSTF*sigZSSTF*sigZSSTF*sigZSSTF/(4.0*f0*f0+k0*k0*sigZSSTF*sigZSSTF*sigZSSTF*sigZSSTF);

   // for normalization at the focus
   a_k0Y=f12Y/k0/k0/sigYSSTF/sigYSSTF-I*0.5/k0*(f0-f12Y/f0);
   a_k0Z=f12Z/k0/k0/sigZSSTF/sigZSSTF-I*0.5/k0*(f0-f12Z/f0);
   gaussTmpY=sigW*sigYSSTF*0.5*csqrt(M_PI/a_k0Y/(1-I*k0*sigYSSTF*sigYSSTF*0.5/f0));
   gaussTmpZ=sigW*sigZSSTF*0.5*csqrt(M_PI/a_k0Z/(1-I*k0*sigZSSTF*sigZSSTF*0.5/f0));


   switch (L->add)  {
   case OFF :
     for(i=0; i<iend+3; i++)  {
       x=(minXSub+i)*dx*lambda0+(f0-f);
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
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+minYSub-jC)*dy*lambda0;
	      argY=tcoefY*sigW/velocityC*(-(x-positionX)+(f0-f))+nY*y;
         arg1Y=-0.25*y*y/a_k0Y-0.25/mY*argY*argY;
         for(k=0; k<kend+3; k++) {
           z=(k-kstart+minZSub-kC)*dz*lambda0;
	        argZ=tcoefZ*sigW/velocityC*(-(x-positionX)+(f0-f))+nZ*z;
           arg1Z=-0.25*z*z/a_k0Z-0.25/mZ*argZ*argZ;
           compA=a0*coefY*coefZ*cexp(arg1Y+arg1Z-I*k0*((x-positionX)-(f0-f)));
           field1[i][j][k]=creal(compA); 
           field2[i][j][k]=creal(compA); 
			}
		 }
     }
	  break;
   case ON :
     for(i=0; i<iend+3; i++)  {
       x=(minXSub+i)*dx*lambda0+(f0-f);
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
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+minYSub-jC)*dy*lambda0;
	      argY=tcoefY*sigW/velocityC*(-(x-positionX)+(f0-f))+nY*y;
         arg1Y=-0.25*y*y/a_k0Y-0.25/mY*argY*argY;
         for(k=0; k<kend+3; k++) {
           z=(k-kstart+minZSub-kC)*dz*lambda0;
	        argZ=tcoefZ*sigW/velocityC*(-(x-positionX)+(f0-f))+nZ*z;
           arg1Z=-0.25*z*z/a_k0Z-0.25/mZ*argZ*argZ;
           compA=a0*coefY*coefZ*cexp(arg1Y+arg1Z-I*k0*((x-positionX)-(f0-f)));
           field1[i][j][k]+=creal(compA); 
           field2[i][j][k]+=creal(compA); 
			}
		 }
     }
	  break;
   } 
}


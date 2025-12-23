#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_randist.h>
#include <fftw3.h>

double maxwellianVelocity(double temperature);
void random1D_sobol(double *x,gsl_qrng *q1);
void random2D_sobol(double *x,double *y,gsl_qrng *q2);
void random3D_sobol(double *x,double *y,double *z,gsl_qrng *q3);
void loadBeamPlasma1D(Domain *D,LoadList *LL,int s,int iteration);
void loadBeamPlasma2D(Domain *D,LoadList *LL,int s,int iteration);
void loadBeamPlasma3D(Domain *D,LoadList *LL,int s,int iteration);
double gaussian_dist(double sig);
void MPI_TransferBeamDen_Xminus(double **f1,int N,int istart,int iend);
void MPI_TransferBeamDen_Xplus(double **f1,int N,int istart,int iend);
void MPI_TransferBeamPhi_Xminus(double **f1,int N,int istart,int iend);
void MPI_TransferBeamPhi_Xplus(double **f1,int N,int istart,int iend);


void loadBeam(Domain *D,LoadList *LL,int s,int iteration)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  switch((LL->type-1)*3+D->dimension)  {
  //1D
  case ((Polygon-1)*3+1):
	  ;
    break;

  //2D
  case (Polygon-1)*3+2:
          ;
    break;

  case (Defined-1)*3+2:
    ;
    break;

  //3D
  case (Polygon-1)*3+3:
    ;
    break;
  case (Defined-1)*3+3:
    ;
    break;

  case (Beam-1)*3+1:
//    if(iteration==LL->loadingStep) {
      loadBeamPlasma1D(D,LL,s,iteration);
//    }    else ;
//    MPI_Barrier(MPI_COMM_WORLD);
    break;

  case (Beam-1)*3+2:
//    if(iteration==LL->loadingStep) {
      loadBeamPlasma2D(D,LL,s,iteration);
//    }    else ;
//    MPI_Barrier(MPI_COMM_WORLD);
    break;
  case (Beam-1)*3+3:
//    if(iteration==LL->loadingStep) {
      loadBeamPlasma3D(D,LL,s,iteration);
//    }    else ;
//    MPI_Barrier(MPI_COMM_WORLD);
    break;


  default:
    ;
  }
}


void loadBeamPlasma1D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,jC,kC;
   int minYSub,maxYSub,totalCnt;
   double posX,weight,charge,weightCoef,ne,nenergy;
   double px,py,pz,positionX,y,yPrime,z,zPrime,sigYPrime,sigZPrime,gammaY,gammaZ;
   double dt,dx,dy,dz,emitY,emitZ,betaY,betaZ,sigY,sigZ;
   double gamma0,dGam,gamma,sigGam,delGam,lowGam,upGam;
   double distance,vx,delT,sqrt2,x,phase;
   int myrank,nTasks,rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   Particle ***particle;
   particle=D->particle;

   ptclList *New,*p;   

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   minYSub=D->minYSub;
   maxYSub=D->maxYSub;
   j=k=0;

   dx=D->dx; dy=D->dy; dz=D->dz; dt=D->dt;
   jC=0;
   kC=0;

   sqrt2=1.414214;
   gamma0=LL->energy/mc2+1.0;
   dGam=LL->spread*gamma0;
   emitY=LL->emitY/gamma0;  
   emitZ=LL->emitZ/gamma0;
   gammaY=(1+LL->alphaY*LL->alphaY)/LL->betaY;
   gammaZ=(1+LL->alphaZ*LL->alphaZ)/LL->betaZ;
   sigY=sqrt(emitY/gammaY)/D->lambda;
   sigZ=sqrt(emitZ/gammaZ)/D->lambda;
   sigYPrime=sqrt(emitY*gammaY);
   sigZPrime=sqrt(emitZ*gammaZ);

   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;

   charge=LL->charge;
   vx=sqrt(gamma0*gamma0-1.0)/gamma0;				//normalized

   srand(myrank+1);
   intNum=(int)LL->numberInCell;
   lowGam=1e10; upGam=-1e10;

   //position define      
   gsl_qrng *q2 = gsl_qrng_alloc (gsl_qrng_sobol,2);
   for(i=istart; i<iend; i++)
   {
     x=(i-istart+D->minXDomain)*D->dx*D->lambda;
     ne=0.0;
     if(LL->gaussMode==OFF) {
       for(l=0; l<LL->xnodes-1; l++) {
         posX=(double)(i+D->minXSub-istart);
         if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1])  {
           ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
         } else ;
       }
     } else if(LL->gaussMode==ON) {
       posX=LL->posX;
       phase=pow((x-posX)/LL->sigX,LL->gaussPower);
       ne=exp(-phase);
     } else ;

     //energy chirp
     nenergy=1.0+(x-LL->posX)*LL->eChirp/LL->energy;

     totalCnt=ne*intNum;
     weight=1.0/LL->numberInCell*weightCoef;

     cnt=0;
     while(cnt<totalCnt)  {      
       positionX=randomValue(1.0);
       y=gaussian_dist(sigY*sqrt2);
       yPrime=gaussian_dist(sigYPrime*sqrt2);
       z=gaussian_dist(sigZ*sqrt2);
       zPrime=gaussian_dist(sigZPrime*sqrt2);
       sigGam=dGam*0.5*1.201122;		// 1.201122=1/sqrt(ln(2))
       delGam=gaussian_dist(sigGam);
       gamma=gamma0*nenergy+delGam;

       px=sqrt((gamma*gamma-1.0)/(1.0+yPrime*yPrime+zPrime*zPrime));
       py=yPrime*px;
       pz=zPrime*px;
       
       New = (ptclList *)malloc(sizeof(ptclList)); 
       New->next = particle[i][j][k].head[s]->pt;
       particle[i][j][k].head[s]->pt = New;
 
       New->x = positionX; New->y = 0.0; New->z = 0.0;
       New->oldX=i+positionX-px/gamma*dt/dx;
       New->oldY=0.0; New->oldZ=0.0;
  
       New->E1=New->E2=New->E3=0.0;
       New->B1=New->B2=New->B3=0.0;
       New->p1=px;      New->p2=py;       New->p3=pz;
       New->p1Old1=0.0; New->p2Old1=0.0; New->p3Old1=0.0;
       New->weight=weight;
       New->charge=charge;
       D->index+=1;
       New->index=D->index;            
       New->core=myrank;            

       gamma=sqrt(1.0+px*px+py*py+pz*pz);
       if(gamma<=lowGam) lowGam=gamma; else ;
       if(gamma>upGam) upGam=gamma; else ;

       cnt++; 	   
     }		//end of while(cnt)
   }			//End of for(i,j) 
   gsl_qrng_free(q2);

   double tmpMin[nTasks],tmpMax[nTasks],intervalGam;

   if(fabs(lowGam)==1e10) lowGam=0.0; else ;
   if(fabs(upGam)==1e10) upGam=0.0; else ;
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Gather(&lowGam,1,MPI_DOUBLE,tmpMin,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   MPI_Bcast(tmpMin,nTasks,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Gather(&upGam,1,MPI_DOUBLE,tmpMax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(tmpMax,nTasks,MPI_DOUBLE,0,MPI_COMM_WORLD);

   lowGam=1e10;      upGam=-1e10;
   for(rank=0; rank<nTasks; rank++) {
     if(lowGam>=tmpMin[rank] && tmpMin[rank]>0) lowGam=tmpMin[rank]; else ;
     if(upGam<tmpMax[rank])   upGam=tmpMax[rank]; else ;
   }
   intervalGam=(upGam-lowGam)/(LL->numGam*1.0);

   if(myrank==0) {
     printf("minGamma=%g, maxGamma=%g\n",lowGam,upGam);
   } else ;


}


void loadBeamPlasma2D(Domain *D,LoadList *LL,int s,int iteration)
{
   int ii,jj,i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,jC,kC,n;
   int minXSub,minYSub,maxYSub,totalCnt,i1,j1,intXc,intYc,flag,nxSub,nySub,nx,ny,nz,N;
   double posX,posY,posZ,tmp,weight,charge,weightCoef,ne,nenergy,testY,shiftY,shiftZ;
   double px,py,pz,positionX,positionY,yPrime,zPrime,sigYPrime,sigZPrime;
   double dt,dx,dy,dz,emitY,emitZ,sigY,sigZ,gammaY,gammaZ,gamma0,dGam,gamma,sigGam,delGam;
   double distanceY,distanceZ,delTY,delTZ,wx[2],wy[2],density;
   double alpha,beta,aa,bb,x,y,z,xPrime,absX,phase,vx,integralG,gam,dX,dY,beta2,norm,signAlphaY,signAlphaZ;
   double coeff,WX[2],WY[2],xc,yc,xcc,ycc,vy,vz;
   double **den,**phi,**fftDen,**fftPhi,*shareXZ,*sendXZ,*shareXY,*sendXY;
   FILE *out;
   char name[100];

   MPI_Status status;
   int myrank,nTasks,rankX,rank;
   double recv;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   int minmax[nTasks+1],coreCntZ[nTasks],corePositionZ[nTasks],coreCntY[nTasks],corePositionY[nTasks];

   Particle ***particle;
   particle=D->particle;
   ptclList *New,*p;   
   fftw_plan fftp;

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   minXSub=D->minXSub;
   minYSub=D->minYSub; maxYSub=D->maxYSub;
   nx=D->nx; ny=D->ny; nz=D->nz; k=0;
   dx=D->dx; dy=D->dy; dz=D->dz; dt=D->dt;
   dX=D->dx*D->lambda; dY=D->dy*D->lambda;
   shiftY=LL->posY/D->lambda;
   shiftZ=LL->posZ/D->lambda;

   coeff=LL->density/LL->criticalDensity;

   gamma0=LL->energy/mc2;
   dGam=LL->spread*gamma0;
   emitY=LL->emitY/gamma0;  
   emitZ=LL->emitZ/gamma0;
   gammaY=(1+LL->alphaY*LL->alphaY)/LL->betaY;   
   gammaZ=(1+LL->alphaZ*LL->alphaZ)/LL->betaZ;   
   signAlphaY=LL->alphaY/fabs(LL->alphaY);
   signAlphaZ=LL->alphaZ/fabs(LL->alphaZ);
   
   sigY=sqrt(emitY/gammaY)/D->lambda;
   sigZ=sqrt(emitZ/gammaZ)/D->lambda;
   sigYPrime=sqrt(emitY*gammaY);
   sigZPrime=sqrt(emitZ*gammaZ);
   
   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;

   charge=LL->charge;
//   distanceY=sqrt((LL->betaY-1.0/gammaY)/gammaY)/D->lambda;   
//   distanceZ=sqrt((LL->betaZ-1.0/gammaZ)/gammaZ)/D->lambda;   
   distanceY=LL->alphaY/gammaY/D->lambda+LL->focalLy/D->lambda;   
   distanceZ=LL->alphaZ/gammaZ/D->lambda+LL->focalLz/D->lambda;   
   vx=sqrt(gamma0*gamma0-1.0)/gamma0;				//normalized
   if(vx==0.0) { delTY=delTZ=0.0; }
   else  {
     delTY=distanceY/vx;						//normalized
     delTZ=distanceZ/vx;						//normalized
   }
   density=LL->density;


   srand(myrank+1);
   intNum=(int)LL->numberInCell;

   //position define     
   double sum,coef,v[4],tmpMin[nTasks],tmpMax[nTasks],lowGam,upGam,intervalGam;
	double minGam,maxGam;
	const gsl_rng_type * T;

   gsl_rng *ran;
   gsl_rng_env_setup();
   T = gsl_rng_default;
   ran = gsl_rng_alloc(T);
   gsl_qrng *q1=gsl_qrng_alloc(gsl_qrng_reversehalton,4);
   lowGam=1e10; upGam=-1e10;
   for(i=istart; i<iend; i++)
   {
     x=(i-istart+minXSub+D->minXDomain)*dX;
     ne=0.0;
     if(LL->gaussMode==OFF) {
       for(l=0; l<LL->xnodes-1; l++) {
         posX=(double)(i+D->minXSub-istart+D->minXDomain);
         if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1]) {
           ne=(LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l];
         } else ;
       }
     } else if(LL->gaussMode==ON) {
       posX=LL->posX;
       phase=pow((x-posX)/LL->sigX,LL->gaussPower);
       if(fabs(x-posX)>4*LL->sigX)  ne=0.0;
       else                         ne=exp(-phase*0.5);
     } else ;

     //energy chirp
     nenergy=1.0+(x-LL->posX)*LL->eChirp/LL->energy;       
    
     totalCnt=(int)(ne*sigY/dy*intNum*sqrt(2.0*M_PI));
     weight=1.0/LL->numberInCell*weightCoef;

     cnt=0;
     while(cnt<totalCnt)  {      
       flag=0;
       while(flag==0) {
         gsl_qrng_get(q1,v);
         sum=0.0;
         for(ii=0; ii<4; ii++) sum+=v[ii];
         if(sum==0.0) flag=0; else flag=1 ;
       }

       coef=sqrt(-2.0*log(v[0]));
       y=coef*cos(2*M_PI*v[1]);
       yPrime=coef*sin(2*M_PI*v[1]);
       y*=sigY;
       yPrime*=sigYPrime;

       coef=sqrt(-2.0*log(v[2]));
       z=coef*cos(2*M_PI*v[3]);
       zPrime=coef*sin(2*M_PI*v[3]);
       z*=sigZ;
       zPrime*=sigZPrime;	   

       positionX=randomValue(1.0);
       sigGam=dGam*0.5*1.201122;		// 1.201122=1/sqrt(ln(2))
       delGam=gaussian_dist(sigGam);
       gamma=gamma0*nenergy+delGam;

       px=sqrt((gamma*gamma-1.0)/(1.0+yPrime*yPrime+zPrime*zPrime));
       py=yPrime*px;
       pz=zPrime*px;
       
       y-=delTY*py/gamma;
       z-=delTZ*pz/gamma;
       y+=shiftY;
       z+=shiftZ;

       testY=y/dy;
       if(testY>=minYSub && testY<maxYSub) {
         j=(int)(y/dy-minYSub+jstart);
         positionY=y/dy-minYSub-(int)(y/dy-minYSub);

         New = (ptclList *)malloc(sizeof(ptclList)); 
         New->next = particle[i][j][k].head[s]->pt;
         particle[i][j][k].head[s]->pt = New;
 
         New->x = positionX; New->y = positionY; New->z = 0.0;
         New->oldX=i+positionX-px/gamma*dt/dx;
         New->oldY=j+positionY-py/gamma*dt/dy;
         New->oldZ=0.0;
  
         New->E1=New->E2=New->E3=0.0;
         New->B1=New->B2=New->B3=0.0;
         New->p1=px;      	New->p2=py;       New->p3=pz;
         New->p1Old1=px; 	New->p2Old1=py; 	New->p3Old1=pz;
         New->weight=weight;
         New->charge=charge;
         D->index+=1;
         New->index=D->index;            
         New->core=myrank;            

         gamma=sqrt(1.0+px*px+py*py+pz*pz);
         if(gamma<=lowGam) lowGam=gamma; else ;
         if(gamma>upGam) upGam=gamma; else ;

         //current
//         xc=i+positionX+0.5*px/gamma*dt/dx;
//         yc=j+positionY+0.5*py/gamma*dt/dy;
//         intXc=(int)xc;   intYc=(int)yc;
//         xcc=xc-intXc;    ycc=yc-intYc;
//         xc=i+positionX+0.5; yc=j+positionY+0.5;
//         intXc=(int)xc;      intYc=(int)yc;
//         xcc=xc-intXc;       ycc=yc-intYc;
         wx[1]=x;          wx[0]=1.0-wx[1];
         wy[1]=y;          wy[0]=1.0-wy[1];
//         WX[1]=xcc;        WX[0]=1.0-WX[1];
//         WY[1]=ycc;        WY[0]=1.0-WY[1];
//         i1=(int)(xc)-1;
//         j1=(int)(yc)-1;

         vx=px/gamma; vy=py/gamma; vz=pz/gamma;
//			for(jj=0; jj<2; jj++)
//           D->Jx[i][j+jj][k]+=coeff*charge*weight*wy[jj]*vx;
//			for(ii=0; ii<2; ii++)
//           D->Jy[i+ii][j][k]+=coeff*charge*weight*wx[ii]*vy;
//         for(ii=0; ii<2; ii++)
//           for(jj=0; jj<2; jj++) {
//             D->CurX[i1+ii][intYc+jj][k]-=coeff*charge*WX[ii]*wy[jj]*2.0*M_PI*vx;
//             D->CurY[intXc+ii][j1+jj][k]-=coeff*charge*wx[ii]*WY[jj]*2.0*M_PI*vy;
//             D->CurZ[intXc+ii][intYc+jj][k]-=coeff*charge*wx[ii]*wy[jj]*2.0*M_PI*vz;
//             D->JxOld[i1-1+ii][j+jj][k]+=coeff*charge*weight*WX[ii]*wy[jj]*vx;
//             D->JyOld[i+ii][j1-1+jj][k]+=coeff*charge*weight*wx[ii]*WY[jj]*vy;
//             D->JzOld[i+ii][j+jj][k]+=coeff*charge*weight*wx[ii]*wy[jj]*vz;
//             D->Jz[i+ii][j+jj][k]+=coeff*charge*weight*wx[ii]*wy[jj]*vz;
//           }
       } else ;
       cnt++; 	   
     }		//end of while(cnt)
//     gsl_qrng_free(q2);
   }			//End of for(i,j)
   gsl_qrng_free(q1);
   gsl_rng_free(ran);   

   if(fabs(lowGam)==1e10) lowGam=0.0; else ;
   if(fabs(upGam)==1e10) upGam=0.0; else ;
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Gather(&lowGam,1,MPI_DOUBLE,tmpMin,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//if(myrank==0) {
//  for(rank=0; rank<nTasks; rank++)
//    printf("tmpMin[%d]=%g\n",rank,tmpMin[rank]);
//
//}

   MPI_Bcast(tmpMin,nTasks,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Gather(&upGam,1,MPI_DOUBLE,tmpMax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(tmpMax,nTasks,MPI_DOUBLE,0,MPI_COMM_WORLD);

   lowGam=1e10;      upGam=-1e10;
   for(rank=0; rank<nTasks; rank++) {
     if(lowGam>=tmpMin[rank] && tmpMin[rank]>0) lowGam=tmpMin[rank]; else ;
     if(upGam<tmpMax[rank])   upGam=tmpMax[rank]; else ;
   }
   intervalGam=(upGam-lowGam)/(LL->numGam*1.0);

   if(myrank==0) {
     printf("minGamma=%g, maxGamma=%g\n",lowGam,upGam);
   } else ;




   // memory creation
   nxSub=D->nxSub;
   ny=D->ny+2;
   N=ny;
   den=(double **)malloc((nxSub+5)*sizeof(double *));
   phi=(double **)malloc((nxSub+5)*sizeof(double *));
   fftDen=(double **)malloc((nxSub+5)*sizeof(double *));
   fftPhi=(double **)malloc((nxSub+5)*sizeof(double *));
   for(i=0; i<nxSub+5; i++) {
     phi[i]=(double *)malloc(N*sizeof(double ));
     den[i]=(double *)malloc(N*sizeof(double ));
     fftDen[i]=(double *)malloc(N*sizeof(double ));
     fftPhi[i]=(double *)malloc(N*sizeof(double ));
   }
   for(i=0; i<nxSub+5; i++)
     for(j=0; j<N; j++) {
       den[i][j]=0.0;
       phi[i][j]=0.0;
       fftDen[i][j]=0.0;
       fftPhi[i][j]=0.0;
     }  
   for(i=0; i<nTasks; i++) minmax[i]=0;
   MPI_Allgather(&D->minXSub,1,MPI_INT,minmax,1,MPI_INT,MPI_COMM_WORLD);
   minmax[nTasks]=nx;
   for(i=0; i<nTasks; i++) {
     coreCntZ[i]=(minmax[i+1]-minmax[i])*nz;
     corePositionZ[i]=minmax[i]*nz;
     coreCntY[i]=(minmax[i+1]-minmax[i])*ny;
     corePositionY[i]=minmax[i]*ny;
   }
   
   shareXY=(double *)malloc(nx*ny*sizeof(double ));
   sendXY=(double *)malloc(nxSub*ny*sizeof(double ));
   shareXZ=(double *)malloc(nx*nz*sizeof(double ));
   sendXZ=(double *)malloc(nxSub*nz*sizeof(double ));

   for(n=0; n<LL->numGam; n++) {
     minGam=lowGam+n*intervalGam;
     maxGam=minGam+intervalGam;
     gamma0=(minGam+maxGam)*0.5;

     for(i=0; i<nxSub+5; i++)
       for(j=0; j<N; j++) {
         den[i][j]=0.0;
         phi[i][j]=0.0;
         fftDen[i][j]=0.0;
         fftPhi[i][j]=0.0;
       }  
     for(i=istart; i<iend; i++)
       for(j=jstart; j<jend; j++) {
         for(s=0; s<D->nSpecies; s++) {
           p=particle[i][j][k].head[s]->pt;
           while(p) {
             px=p->p1; py=p->p2; pz=p->p3; weight=p->weight;
             positionX=p->x; positionY=p->y;

             gamma=sqrt(1.0+px*px+py*py+pz*pz);
				 if(gamma>=minGam && gamma<maxGam) {
               wx[0]=1.0-positionX; wx[1]=positionX;
               wy[0]=1.0-positionY; wy[1]=positionY;
               for(ii=0; ii<2; ii++)
                 for(jj=0; jj<2; jj++)
                   den[i+ii][j-jstart+1+jj]+=density*weight*wx[ii]*wy[jj]*eCharge;	//real scale
				 } else ;

             p=p->next;
			  }
			}
		 }

     if(D->L>1)  {
       MPI_TransferBeamDen_Xplus(den,N,istart,iend);
       MPI_TransferBeamDen_Xminus(den,N,istart,iend);
     }  else     ;


     // save den file
//     for(i=istart; i<iend; i++)
//       for(j=0; j<ny; j++)
//         sendXY[j+(i-istart)*ny]=den[i][j];   
//     MPI_Allgatherv(sendXY,nxSub*ny,MPI_DOUBLE,shareXY,coreCntY,corePositionY,MPI_DOUBLE,MPI_COMM_WORLD);
//     MPI_Barrier(MPI_COMM_WORLD);

//     if(myrank==0) {   
//       sprintf(name,"den_XY");
//       out=fopen(name,"w");
//       for(i=0; i<nx; i++) {
//         x=(D->minXDomain+i)*dX;
//         for(j=0; j<ny; j++) {
//           y=(D->minYDomain+j)*dY;
//           fprintf(out,"%g %g %g\n",x,y,shareXY[j+i*ny]);
//         }
//         fprintf(out,"\n");
//       }
//       fclose(out);
//       printf("%s is made\n",name);
//     }  else ;



     for(i=istart; i<iend; i++) {
       fftp = fftw_plan_r2r_1d(N,den[i],fftDen[i],FFTW_RODFT00,FFTW_ESTIMATE);
       fftw_execute(fftp);
       fftw_destroy_plan(fftp);
     }   

     // fftPhi calculation   
     for(j=0; j<ny; j++) {
       alpha=(j+1.0)*M_PI/(N*1.0)/dY;

       for(i=istart; i<iend; i++)
         sendXZ[(i-istart)*nz]=fftDen[i][j*nz];
       MPI_Allgatherv(sendXZ,nxSub*nz,MPI_DOUBLE,shareXZ,coreCntZ,corePositionZ,MPI_DOUBLE,MPI_COMM_WORLD);
       MPI_Barrier(MPI_COMM_WORLD);

       for(i=istart; i<iend; i++) {
         x=(D->minXDomain+i-istart+minXSub)*dX;

         gam=alpha;  
         coef=0.5*gamma0/eps0/gam;
         sum=0.0;
         for(ii=0; ii<nx; ii++) {
           xPrime=(D->minXDomain+ii)*dX;
           absX=fabs(x-xPrime);
           if(i-istart+minXSub==ii) {
             aa=2.0/gamma0/gam;
             bb=exp(-0.5*gamma0*gam*dX);
             integralG=aa*(1.0-bb);
           } else  {
             aa=exp(gamma0*gam*(dX*0.5-absX));
             bb=exp(-1.0*gamma0*gam*(dX*0.5+absX));
             integralG=(aa-bb)/gamma0/gam;
           }
           tmp=shareXZ[ii*nz]/(1.0+ny)*integralG;
           sum+=tmp;
         }
         fftPhi[i][j*nz]=sum*coef;
       }               //End of for(i<nx)
     }

     for(i=istart; i<iend; i++) {
       fftp = fftw_plan_r2r_1d(N,fftPhi[i],phi[i],FFTW_RODFT00,FFTW_ESTIMATE);
       fftw_execute(fftp);
       fftw_destroy_plan(fftp);
     }
     if(D->L>1)  {
       MPI_TransferBeamPhi_Xplus(phi,N,istart,iend);
       MPI_TransferBeamPhi_Xminus(phi,N,istart,iend);
     }  else     ;

     norm=0.5*eCharge/eMass/velocityC/velocityC;	//0.5 is by fftw3

     // save phi file
//     for(i=istart; i<iend; i++)
//       for(j=0; j<ny; j++)
//         sendXY[j+(i-istart)*ny]=phi[i][j];   
//     MPI_Allgatherv(sendXY,nxSub*ny,MPI_DOUBLE,shareXY,coreCntY,corePositionY,MPI_DOUBLE,MPI_COMM_WORLD);
//     MPI_Barrier(MPI_COMM_WORLD);
   
//     if(myrank==0) {   
//       sprintf(name,"phi_XY");
//       out=fopen(name,"w");
//       for(i=0; i<nx; i++) {
//         x=(D->minXDomain+i)*dX;
//         for(j=0; j<ny; j++) {
//           y=(D->minYDomain+j)*dY;
//           fprintf(out,"%g %g %g\n",x,y,shareXY[j+i*ny]*norm);
//         }
//         fprintf(out,"\n");
//       }
//       fclose(out);
//       printf("%s is made\n",name);
//     } else ;

    //field assign
    beta2=1.0-1.0/gamma0/gamma0;
    beta=sqrt(beta2);
    if(D->fieldType==Yee || D->fieldType==Pukhov) {
      for(i=istart; i<iend; i++)
        for(j=1; j<ny-1; j++) {
          D->Ex[i][j-1+jstart][k]+=0.5/M_PI/gamma0/gamma0*(phi[i][j]-phi[i-1][j])*norm/dx;
          D->Ey[i][j-1+jstart][k]+=0.5/M_PI*(phi[i-1][j+1]-phi[i-1][j])*norm/dy;
//          D->Bz[i][j-1+jstart][k]+=0.25/M_PI*(phi[i][j+1]+phi[i-1][j+1]-phi[i][j]-phi[i-1][j])*norm/dy*beta;
          D->Bz[i][j-1+jstart][k]+=0.5/M_PI*(phi[i-1][j+1]-phi[i-1][j])*norm/dy*beta;
//          D->Bz[i][j-1+jstart][k]+=0.25/M_PI*(phi[i+1][j+1]+phi[i][j+1]-phi[i+1][j]-phi[i][j])*norm/dy*beta;
        }

      if(D->L>1)  {
        D->shareF[0]=D->Ex;
        D->shareF[1]=D->Ey;
        D->shareF[2]=D->Bz;
        MPI_TransferF_Xminus(D,3);
        MPI_TransferF_Xplus(D,3);
        D->Ex=D->shareF[0];
        D->Ey=D->shareF[1];
        D->Ez=D->shareF[2];
      } else      ;
    } else if(D->fieldType==Split) {
      for(i=istart; i<iend; i++)
        for(j=1; j<ny-1; j++) {
          D->Ex[i][j-1+jstart][k]+=0.5/M_PI/gamma0/gamma0*(phi[i][j]-phi[i-1][j])*norm/dx;
          D->Pr[i][j-1+jstart][k]+=0.25/M_PI*(phi[i][j+1]+phi[i-1][j+1]-phi[i][j]-phi[i-1][j])*norm/dy;
        }
      if(D->L>1)  {
        D->shareF[0]=D->Ex;
        D->shareF[1]=D->Pr;
        MPI_TransferF_Xminus(D,2);
        MPI_TransferF_Xplus(D,2);
        D->Ex=D->shareF[0];
        D->Pr=D->shareF[1];
      } else      ;
    } else ;

    if(myrank==0) printf("During loading beam, step/numGam=%d/%d\n",n,LL->numGam); else ;
	}		//End of for(n)


   for(i=0; i<nxSub+5; i++) {
     free(den[i]);
     free(phi[i]);
     free(fftDen[i]);
     free(fftPhi[i]);
   }
   free(den);
   free(phi);
   free(fftDen);
   free(fftPhi);
   free(shareXZ);
   free(sendXZ);
   free(shareXY);
   free(sendXY);
}


void loadBeamPlasma3D(Domain *D,LoadList *LL,int s,int iteration)
{
   int ii,jj,kk,i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,jC,kC,n,prog;
   int minXSub,minYSub,minZSub,maxYSub,maxZSub,totalCnt,i1,j1,intXc,intYc,flag,nxSub,nySub,nx,ny,nz,N;
   double posX,posY,posZ,tmp,weight,charge,weightCoef,ne,nenergy,testY,testZ;
   double px,py,pz,positionX,positionY,positionZ,yPrime,zPrime,sigYPrime,sigZPrime;
   double dt,dx,dy,dz,emitY,emitZ,sigY,sigZ,gammaY,gammaZ,gamma0,dGam,gamma,sigGam,delGam;
   double distanceY,distanceZ,delTY,delTZ,wx[2],wy[2],wz[2],density;
   double alpha,beta,aa,bb,x,y,z,xPrime,absX,phase,vx,integralG,gam,dX,dY,dZ,beta2,norm;
   double coeff,WX[2],WY[2],xc,yc,xcc,ycc,vy,vz;
   double **den,**phi,**fftDen,**fftPhi,*shareXZ,*sendXZ,*shareXY,*sendXY;
   FILE *out;
   char name[100];

   MPI_Status status;
   int myrank,nTasks,rankX,rank;
   double recv;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   int minmax[nTasks+1],coreCntZ[nTasks],corePositionZ[nTasks],coreCntY[nTasks],corePositionY[nTasks];

   Particle ***particle;
   particle=D->particle;
   ptclList *New,*p;   
   fftw_plan fftp;

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   minXSub=D->minXSub;
   minYSub=D->minYSub; maxYSub=D->maxYSub;
   minZSub=D->minZSub; maxZSub=D->maxZSub;
   nx=D->nx; ny=D->ny; nz=D->nz;
   dx=D->dx; dy=D->dy; dz=D->dz; dt=D->dt;
   dX=D->dx*D->lambda; dY=D->dy*D->lambda; dZ=D->dz*D->lambda;
   jC=0;
   kC=0;

   coeff=LL->density/LL->criticalDensity;

   gamma0=LL->energy/mc2;
   dGam=LL->spread*gamma0;
   emitY=LL->emitY/gamma0;  
   emitZ=LL->emitZ/gamma0;
   gammaY=(1+LL->alphaY*LL->alphaY)/LL->betaY;   
   gammaZ=(1+LL->alphaZ*LL->alphaZ)/LL->betaZ;   
   
   sigY=sqrt(emitY/gammaY)/D->lambda;
   sigZ=sqrt(emitZ/gammaZ)/D->lambda;
   sigYPrime=sqrt(emitY*gammaY);
   sigZPrime=sqrt(emitZ*gammaZ);
   
   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;

   charge=LL->charge;
//   distanceY=sqrt((LL->betaY-1.0/gammaY)/gammaY)/D->lambda;   
//   distanceZ=sqrt((LL->betaZ-1.0/gammaZ)/gammaZ)/D->lambda;   
   distanceY=LL->alphaY/gammaY/D->lambda+LL->focalLy/D->lambda;   
   distanceZ=LL->alphaZ/gammaZ/D->lambda+LL->focalLz/D->lambda;   
   vx=sqrt(gamma0*gamma0-1.0)/gamma0;				//normalized
   if(vx==0.0) { delTY=delTZ=0.0; }
   else  {
     delTY=distanceY/vx;						//normalized
     delTZ=distanceZ/vx;						//normalized
   }
   density=LL->density;


   srand(myrank+1);
   intNum=(int)LL->numberInCell;

   //position define     
   double sum,coef,v[4],tmpMin[nTasks],tmpMax[nTasks],lowGam,upGam,intervalGam;
	double minGam,maxGam;
	const gsl_rng_type * T;

   gsl_rng *ran;
   gsl_rng_env_setup();
   T = gsl_rng_default;
   ran = gsl_rng_alloc(T);
   gsl_qrng *q1=gsl_qrng_alloc(gsl_qrng_reversehalton,4);
	lowGam=1e10; upGam=-1e10;
   for(i=istart; i<iend; i++)
   {
     x=(i-istart+minXSub+D->minXDomain)*dX;
     ne=0.0;
     if(LL->gaussMode==OFF) {
       for(l=0; l<LL->xnodes-1; l++) {
         posX=(double)(i+D->minXSub-istart+D->minXDomain);
         if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1]) {
           ne=(LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l];
         } else ;
       }
     } else if(LL->gaussMode==ON) {
       posX=LL->posX;
       phase=pow((x-posX)/LL->sigX,LL->gaussPower);
       ne=exp(-phase*0.5);
     } else ;

     //energy chirp
     nenergy=1.0+(x-LL->posX)*LL->eChirp/LL->energy;       
    
     totalCnt=ne*sigY/dy*sigZ/dz*intNum*2.0*M_PI;
     weight=1.0/LL->numberInCell*weightCoef;

     cnt=0;
     while(cnt<totalCnt)  {      
       flag=0;
       while(flag==0) {
         gsl_qrng_get(q1,v);
         sum=0.0;
         for(ii=0; ii<4; ii++) sum+=v[ii];
         if(sum==0.0) flag=0; else flag=1 ;
       }

       coef=sqrt(-2.0*log(v[0]));
       y=coef*cos(2*M_PI*v[1]);
       yPrime=coef*sin(2*M_PI*v[1]);
       y*=sigY;
       yPrime*=sigYPrime;

       coef=sqrt(-2.0*log(v[2]));
       z=coef*cos(2*M_PI*v[3]);
       zPrime=coef*sin(2*M_PI*v[3]);
       z*=sigZ;
       zPrime*=sigZPrime;	   

       positionX=randomValue(1.0);
       sigGam=dGam*0.5*1.201122;		// 1.201122=1/sqrt(ln(2))
       delGam=gaussian_dist(sigGam);
       gamma=gamma0*nenergy+delGam;

       px=sqrt((gamma*gamma-1.0)/(1.0+yPrime*yPrime+zPrime*zPrime));
       py=yPrime*px;
       pz=zPrime*px;
       
       y-=delTY*py/gamma;
       z-=delTZ*pz/gamma;

       testY=y/dy; testZ=z/dz;
       if(testY>=minYSub && testY<maxYSub && testZ>=minZSub && testZ<maxZSub) {
         j=(int)(y/dy-minYSub+jstart);
         k=(int)(z/dz-minZSub+kstart);
         positionY=y/dy-minYSub-(int)(y/dy-minYSub);
         positionZ=z/dz-minZSub-(int)(z/dz-minZSub);

         New = (ptclList *)malloc(sizeof(ptclList)); 
         New->next = particle[i][j][k].head[s]->pt;
         particle[i][j][k].head[s]->pt = New;
 
         New->x = positionX; New->y = positionY; New->z = positionZ;
         New->oldX=i+positionX-px/gamma*dt/dx;
         New->oldY=j+positionY-py/gamma*dt/dy;
         New->oldZ=k+positionZ-py/gamma*dt/dy;
  
         New->E1=New->E2=New->E3=0.0;
         New->B1=New->B2=New->B3=0.0;
         New->p1=px;      New->p2=py;       New->p3=pz;
         New->p1Old1=0.0; New->p2Old1=0.0; New->p3Old1=0.0;
         New->weight=weight;
         New->charge=charge;
         D->index+=1;
         New->index=D->index;            
         New->core=myrank;            

         gamma=sqrt(1.0+px*px+py*py+pz*pz);
         if(gamma<=lowGam) lowGam=gamma; else ;
         if(gamma>upGam) upGam=gamma; else ;
         //current
//         xc=i+positionX+0.5*px/gamma*dt/dx;
//         yc=j+positionY+0.5*py/gamma*dt/dy;
//         intXc=(int)xc;   intYc=(int)yc;
//         xcc=xc-intXc;    ycc=yc-intYc;
//         xc=i+positionX+0.5; yc=j+positionY+0.5;
//         intXc=(int)xc;      intYc=(int)yc;
//         xcc=xc-intXc;       ycc=yc-intYc;
//         wx[1]=x;          wx[0]=1.0-wx[1];
//         wy[1]=y;          wy[0]=1.0-wy[1];
//         WX[1]=xcc;        WX[0]=1.0-WX[1];
//         WY[1]=ycc;        WY[0]=1.0-WY[1];
//         i1=(int)(xc)-1;
//         j1=(int)(yc)-1;

//         vx=px/gamma;
//         vy=py/gamma;
//         vz=pz/gamma;
//         for(ii=0; ii<2; ii++)
//           for(jj=0; jj<2; jj++) {
//             D->CurX[i1+ii][intYc+jj][k]-=coeff*charge*WX[ii]*wy[jj]*2.0*M_PI*vx;
//             D->CurY[intXc+ii][j1+jj][k]-=coeff*charge*wx[ii]*WY[jj]*2.0*M_PI*vy;
//             D->CurZ[intXc+ii][intYc+jj][k]-=coeff*charge*wx[ii]*wy[jj]*2.0*M_PI*vz;
//             D->JxOld[i1-1+ii][j+jj][k]+=coeff*charge*weight*WX[ii]*wy[jj]*vx;
//             D->JyOld[i+ii][j1-1+jj][k]+=coeff*charge*weight*wx[ii]*WY[jj]*vy;
//             D->JzOld[i+ii][j+jj][k]+=coeff*charge*weight*wx[ii]*wy[jj]*vz;
//             D->Jx[i1-1+ii][j+jj][k]+=coeff*charge*weight*WX[ii]*wy[jj]*vx;
//             D->Jy[i+ii][j1-1+jj][k]+=coeff*charge*weight*wx[ii]*WY[jj]*vy;
//             D->Jz[i+ii][j+jj][k]+=coeff*charge*weight*wx[ii]*wy[jj]*vz;
//           }
       } else ;
       cnt++; 	   
     }		//end of while(cnt)
//     gsl_qrng_free(q2);
   }			//End of for(i,j)
   gsl_qrng_free(q1);
   gsl_rng_free(ran);   

   if(fabs(lowGam)==1e10) lowGam=0.0; else ;
	if(fabs(upGam)==1e10) upGam=0.0; else ;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(&lowGam,1,MPI_DOUBLE,tmpMin,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(tmpMin,nTasks,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(&upGam,1,MPI_DOUBLE,tmpMax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(tmpMax,nTasks,MPI_DOUBLE,0,MPI_COMM_WORLD);

	lowGam=tmpMin[0]; upGam=tmpMax[0];
	for(rank=1; rank<nTasks; rank++) {
	  if(lowGam>=tmpMin[rank] && tmpMin[rank]!=0) lowGam=tmpMin[rank]; else ;
	  if(upGam<tmpMax[rank])   upGam=tmpMax[rank]; else ;
	}
	intervalGam=(upGam-lowGam)/(LL->numGam*1.0);

	if(myrank==0) {
	  printf("minGamma=%g, maxGamma=%g\n",lowGam,upGam);
	} else ;



   // memory creation
   nxSub=D->nxSub;
   ny=D->ny+2; nz=D->nz+2;
   N=ny*nz;
   den=(double **)malloc((nxSub+5)*sizeof(double *));
   phi=(double **)malloc((nxSub+5)*sizeof(double *));
   fftDen=(double **)malloc((nxSub+5)*sizeof(double *));
   fftPhi=(double **)malloc((nxSub+5)*sizeof(double *));
   for(i=0; i<nxSub+5; i++) {
     phi[i]=(double *)malloc(N*sizeof(double ));
     den[i]=(double *)malloc(N*sizeof(double ));
     fftDen[i]=(double *)malloc(N*sizeof(double ));
     fftPhi[i]=(double *)malloc(N*sizeof(double ));
   }
   for(i=0; i<nxSub+5; i++)
     for(j=0; j<N; j++) {
       den[i][j]=0.0;
       phi[i][j]=0.0;
       fftDen[i][j]=0.0;
       fftPhi[i][j]=0.0;
     }  
   for(i=0; i<nTasks; i++) minmax[i]=0;
   MPI_Allgather(&D->minXSub,1,MPI_INT,minmax,1,MPI_INT,MPI_COMM_WORLD);
   minmax[nTasks]=nx;
   for(i=0; i<nTasks; i++) {
     coreCntZ[i]=(minmax[i+1]-minmax[i])*nz;
     corePositionZ[i]=minmax[i]*nz;
     coreCntY[i]=(minmax[i+1]-minmax[i])*ny;
     corePositionY[i]=minmax[i]*ny;
   }
   
   shareXY=(double *)malloc(nx*ny*sizeof(double ));
   sendXY=(double *)malloc(nxSub*ny*sizeof(double ));
   shareXZ=(double *)malloc(nx*nz*sizeof(double ));
   sendXZ=(double *)malloc(nxSub*nz*sizeof(double ));

   for(n=0; n<LL->numGam; n++) {
     minGam=lowGam+n*intervalGam;
     maxGam=minGam+intervalGam;
     gamma0=(minGam+maxGam)*0.5;

     for(i=0; i<nxSub+5; i++)
       for(j=0; j<N; j++) {
         den[i][j]=0.0;
         phi[i][j]=0.0;
         fftDen[i][j]=0.0;
         fftPhi[i][j]=0.0;
       }  
     for(i=istart; i<iend; i++)
       for(j=jstart; j<jend; j++)
         for(k=kstart; k<kend; k++) 
           for(s=0; s<D->nSpecies; s++) {
             p=particle[i][j][k].head[s]->pt;
             while(p) {
               px=p->p1; py=p->p2; pz=p->p3; weight=p->weight;
               positionX=p->x; positionY=p->y; positionZ=p->z;

               gamma=sqrt(1.0+px*px+py*py+pz*pz);
     	  	  	   if(gamma>=minGam && gamma<maxGam) {
                 wx[0]=1.0-positionX; wx[1]=positionX;
                 wy[0]=1.0-positionY; wy[1]=positionY;
                 wz[0]=1.0-positionZ; wz[1]=positionZ;
                 for(ii=0; ii<2; ii++)
                   for(jj=0; jj<2; jj++)
                     for(kk=0; kk<2; kk++)
                       den[i+ii][k-kstart+1+kk+(j-jstart+1+jj)*nz]+=density*weight*wx[ii]*wy[jj]*wz[kk]*eCharge;	//real scale
				   } else ;

               p=p->next;
			    }
			  }

     if(D->L>1)  {
       MPI_TransferBeamDen_Xplus(den,N,istart,iend);
       MPI_TransferBeamDen_Xminus(den,N,istart,iend);
     }  else     ;

     // save den file
//	  k=nz/2;
//     for(i=istart; i<iend; i++)
//       for(j=0; j<ny; j++)
//         sendXY[j+(i-istart)*ny]=den[i][k+j*nz];   
//     MPI_Allgatherv(sendXY,nxSub*ny,MPI_DOUBLE,shareXY,coreCntY,corePositionY,MPI_DOUBLE,MPI_COMM_WORLD);
//     MPI_Barrier(MPI_COMM_WORLD);

//     if(myrank==0) {   
//       sprintf(name,"den_XY");
//       out=fopen(name,"w");
//       for(i=0; i<nx; i++) {
//         x=(D->minXDomain+i)*dX;
//         for(j=0; j<ny; j++) {
//           y=(D->minYDomain+j)*dY;
//           fprintf(out,"%g %g %g\n",x,y,shareXY[j+i*ny]);
//         }
//         fprintf(out,"\n");
//       }
//       fclose(out);
//       printf("%s is made\n",name);
//     }  else ;

     for(i=istart; i<iend; i++) {
       fftp = fftw_plan_r2r_2d(ny,nz,den[i],fftDen[i],FFTW_RODFT00,FFTW_RODFT00,FFTW_ESTIMATE);
       fftw_execute(fftp);
       fftw_destroy_plan(fftp);
     }   

     // fftPhi calculation   
     for(j=0; j<ny; j++) {
       alpha=(j+1.0)*M_PI/(ny*1.0)/dY;

       for(i=istart; i<iend; i++)
         for(k=0; k<nz; k++)
           sendXZ[k+(i-istart)*nz]=fftDen[i][k+j*nz];
       MPI_Allgatherv(sendXZ,nxSub*nz,MPI_DOUBLE,shareXZ,coreCntZ,corePositionZ,MPI_DOUBLE,MPI_COMM_WORLD);
       MPI_Barrier(MPI_COMM_WORLD);

       for(i=istart; i<iend; i++) {
         x=(D->minXDomain+i-istart+minXSub)*dX;
         for(k=0; k<nz; k++) {
           beta=(k+1.0)*M_PI/(nz*1.0)/dZ;

           gam=sqrt(alpha*alpha+beta*beta);  
           coef=0.5*gamma0/eps0/gam;
           sum=0.0;
           for(ii=0; ii<nx; ii++) {
             xPrime=(D->minXDomain+ii)*dX;
             absX=fabs(x-xPrime);
             if(i-istart+minXSub==ii) {
               aa=2.0/gamma0/gam;
               bb=exp(-0.5*gamma0*gam*dX);
               integralG=aa*(1.0-bb);
             } else  {
               aa=exp(gamma0*gam*(dX*0.5-absX));
               bb=exp(-1.0*gamma0*gam*(dX*0.5+absX));
               integralG=(aa-bb)/gamma0/gam;
             }
             tmp=shareXZ[k+ii*nz]/(1.0+ny)/(1.0+nz)*integralG;
             sum+=tmp;
           }
           fftPhi[i][k+j*nz]=sum*coef;
			}		//End of for(k)
       }       //End of for(i<nx)
		 prog=(int)(100.0*j/(ny*1.0))+1;
		 if(myrank==0 && prog%5==0) printf("progress of fftPhi is %d%% \n",prog); else ;
     }			//End of for(j)

     for(i=istart; i<iend; i++) {
       fftp = fftw_plan_r2r_2d(ny,nz,fftPhi[i],phi[i],FFTW_RODFT00,FFTW_RODFT00,FFTW_ESTIMATE);
       fftw_execute(fftp);
       fftw_destroy_plan(fftp);
     }
     if(D->L>1)  {
       MPI_TransferBeamPhi_Xplus(phi,N,istart,iend);
       MPI_TransferBeamPhi_Xminus(phi,N,istart,iend);
     }  else     ;

     norm=0.25*eCharge/eMass/velocityC/velocityC;	//0.25 is by fftw3

     // save phi file
//	  k=nz/2;
//     for(i=istart; i<iend; i++)
//       for(j=0; j<ny; j++)
//         sendXY[j+(i-istart)*ny]=phi[i][k+j*nz];   
//     MPI_Allgatherv(sendXY,nxSub*ny,MPI_DOUBLE,shareXY,coreCntY,corePositionY,MPI_DOUBLE,MPI_COMM_WORLD);
//     MPI_Barrier(MPI_COMM_WORLD);
   
//     if(myrank==0) {   
//       sprintf(name,"phi_XY");
//       out=fopen(name,"w");
//       for(i=0; i<nx; i++) {
//         x=(D->minXDomain+i)*dX;
//         for(j=0; j<ny; j++) {
//           y=(D->minYDomain+j)*dY;
//           fprintf(out,"%g %g %g\n",x,y,shareXY[j+i*ny]*norm);
//         }
//         fprintf(out,"\n");
//       }
//       fclose(out);
//       printf("%s is made\n",name);
//     } else ;

     //field assign
     beta2=1.0-1.0/gamma0/gamma0;
     beta=sqrt(beta2);
     if(D->fieldType==Yee || D->fieldType==Pukhov) {
       for(i=istart; i<iend; i++)
         for(j=1; j<ny-1; j++) 
           for(k=1; k<nz-1; k++) {
             D->Ex[i][j-1+jstart][k-1+kstart]+=0.5/M_PI/gamma0/gamma0*(phi[i][k+j*nz]-phi[i-1][k+j*nz])*norm/dx;
             D->Ey[i][j-1+jstart][k-1+kstart]+=0.5/M_PI*(phi[i-1][k+(j+1)*nz]-phi[i-1][k+j*nz])*norm/dy;
             D->Ez[i][j-1+jstart][k-1+kstart]+=0.5/M_PI*(phi[i-1][(k+1)+j*nz]-phi[i-1][k+j*nz])*norm/dz;
             D->Bz[i][j-1+jstart][k-1+kstart]+=0.5/M_PI*(phi[i-1][k+(j+1)*nz]-phi[i-1][k+j*nz])*norm/dy*beta;
             D->By[i][j-1+jstart][k-1+kstart]-=0.5/M_PI*(phi[i-1][(k+1)+j*nz]-phi[i-1][k+j*nz])*norm/dz*beta;
           }
       if(D->L>1)  {
        D->shareF[0]=D->Ex;
        D->shareF[1]=D->Ey;
        D->shareF[2]=D->Ez;
        D->shareF[3]=D->Bz;
        D->shareF[4]=D->By;
        MPI_TransferF_Xminus(D,5);
        MPI_TransferF_Xplus(D,5);
        D->Ex=D->shareF[0];
        D->Ey=D->shareF[1];
        D->Ez=D->shareF[2];
        D->Bz=D->shareF[3];
        D->By=D->shareF[4];
       } else      ;
     } else if(D->fieldType==Split) {
       for(i=istart; i<iend; i++)
         for(j=1; j<ny-1; j++) 
           for(k=1; k<nz-1; k++) {
             D->Ex[i][j-1+jstart][k-1+kstart]+=0.5/M_PI/gamma0/gamma0*(phi[i][k+j*nz]-phi[i-1][k+j*nz])*norm/dx;
             D->Pr[i][j-1+jstart][k-1+kstart]+=0.25/M_PI*(phi[i][k+(j+1)*nz]+phi[i-1][k+(j+1)*nz]-phi[i][k+j*nz]-phi[i-1][k+j*nz])*norm/dy;
             D->Sr[i][j-1+jstart][k-1+kstart]+=0.25/M_PI*(phi[i][k+1+j*nz]+phi[i-1][k+1+j*nz]-phi[i][k+j*nz]-phi[i-1][k+j*nz])*norm/dz;
           }
       if(D->L>1)  {
        D->shareF[0]=D->Ex;
        D->shareF[1]=D->Pr;
        D->shareF[2]=D->Sr;
        MPI_TransferF_Xminus(D,3);
        MPI_TransferF_Xplus(D,3);
        D->Ex=D->shareF[0];
        D->Pr=D->shareF[1];
        D->Sr=D->shareF[2];
       } else      ;
     } else ;

     if(myrank==0) printf("During loading beam, step/numGam=%d/%d\n",n,LL->numGam); else ;
   }		//End of for(n)


   for(i=0; i<nxSub+5; i++) {
     free(den[i]);
     free(phi[i]);
     free(fftDen[i]);
     free(fftPhi[i]);
   }
   free(den);
   free(phi);
   free(fftDen);
   free(fftPhi);
   free(shareXZ);
   free(sendXZ);
   free(shareXY);
   free(sendXY);
}


void MPI_TransferBeamDen_Xminus(double **f1,int N,int istart,int iend)
{
    int i,j,num,start;
    int myrank, nTasks;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    num=N*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores
    start=0;
    for(i=0; i<2; i++) {
      for(j=0; j<N; j++) data[j+start]=f1[istart-1-i][j]; start+=N;
    }

    if(myrank%2==0 && myrank!=nTasks-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,myrank+1,myrank+1, MPI_COMM_WORLD,&status);
      start=0;
      for(i=0; i<2; i++) {
        for(j=0; j<N; j++)  f1[iend-1-i][j]+=data[j+start]; start+=N;
      }
    }
    else if(myrank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,myrank-1,myrank,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores
    start=0;
    for(i=0; i<2; i++) {
      for(j=0; j<N; j++) data[j+start]=f1[istart-1-i][j]; start+=N;
    }

    if(myrank%2==1 && myrank!=nTasks-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,myrank+1,myrank+1, MPI_COMM_WORLD,&status);
      start=0;
      for(i=0; i<2; i++) {
        for(j=0; j<N; j++)  f1[iend-1-i][j]+=data[j+start]; start+=N;
      }
    }
    else if(myrank%2==0 && myrank!=0)
       MPI_Send(data,num,MPI_DOUBLE,myrank-1,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}


void MPI_TransferBeamDen_Xplus(double **f1,int N,int istart,int iend)
{
    int i,j,num,start;
    int myrank, nTasks;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    num=N*3;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores
    start=0;
    for(i=0; i<3; i++) {
      for(j=0; j<N; j++) data[j+start]=f1[iend+i][j]; start+=N;
    }

    if(myrank%2==1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,myrank-1,myrank-1, MPI_COMM_WORLD,&status);
      start=0;
      for(i=0; i<2; i++) {
        for(j=0; j<N; j++)  f1[istart+i][j]+=data[j+start];  start+=N;
      }
    }
    else if(myrank%2==0 && myrank!=nTasks-1)
       MPI_Send(data,num,MPI_DOUBLE,myrank+1,myrank,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores
    start=0;
    for(i=0; i<2; i++) {
      for(j=0; j<N; j++) data[j+start]=f1[iend+i][j]; start+=N;
    }

    if(myrank%2==0 && myrank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,myrank-1,myrank-1, MPI_COMM_WORLD,&status);
      start=0;
      for(i=0; i<2; i++) {
        for(j=0; j<N; j++) f1[istart+i][j]+=data[j+start]; start+=N;
      }
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
       MPI_Send(data,num,MPI_DOUBLE,myrank+1,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}




void MPI_TransferBeamPhi_Xminus(double **f1,int N,int istart,int iend)
{
    int i,j,num,start;
    int myrank, nTasks;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    num=N*3;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores
    start=0;
    for(i=0; i<3; i++) { 
      for(j=0; j<N; j++) data[j+start]=f1[istart+i][j]; start+=N;
    }

    if(myrank%2==0 && myrank!=nTasks-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,myrank+1,myrank+1, MPI_COMM_WORLD,&status);
      start=0;
      for(i=0; i<3; i++) { 
        for(j=0; j<N; j++)  f1[iend+i][j]=data[j+start]; start+=N;
      }
    }
    else if(myrank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,myrank-1,myrank,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores
    start=0;
    for(i=0; i<3; i++) { 
      for(j=0; j<N; j++) data[j+start]=f1[istart+i][j]; start+=N; 
    }

    if(myrank%2==1 && myrank!=nTasks-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,myrank+1,myrank+1, MPI_COMM_WORLD,&status);
      start=0;
      for(i=0; i<3; i++) { 
        for(j=0; j<N; j++) f1[iend+i][j]=data[j+start]; start+=N;
      }
    }
    else if(myrank%2==0 && myrank!=0)
       MPI_Send(data,num,MPI_DOUBLE,myrank-1,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_TransferBeamPhi_Xplus(double **f1,int N,int istart,int iend)
{
    int i,j,num,start;
    int myrank, nTasks;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    num=N*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores
    start=0;
    for(i=1; i<3; i++) { 
      for(j=0; j<N; j++) data[j+start]=f1[iend-i][j]; start+=N;
    }

    if(myrank%2==1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,myrank-1,myrank-1, MPI_COMM_WORLD,&status);
      start=0;
      for(i=1; i<3; i++) { 
        for(j=0; j<N; j++) f1[istart-i][j]=data[j+start]; start+=N;
      }
    }
    else if(myrank%2==0 && myrank!=nTasks-1)
       MPI_Send(data,num,MPI_DOUBLE,myrank+1,myrank,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores
    start=0;
    for(i=1; i<3; i++) { 
      for(j=0; j<N; j++)  data[j+start]=f1[iend-i][j]; start+=N;
    }

    if(myrank%2==0 && myrank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,myrank-1,myrank-1, MPI_COMM_WORLD,&status);
      start=0;
      for(i=1; i<3; i++) { 
        for(j=0; j<N; j++) f1[istart-i][j]=data[j+start]; start+=N;
      }
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
       MPI_Send(data,num,MPI_DOUBLE,myrank+1,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}   

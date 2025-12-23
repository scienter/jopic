#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

double maximum(double x1,double x2)  {
   double result;
   if(x1>=x2)      result=x1;
   else           result=x2;
   return result;
}

double minimum(double x1,double x2)  {
   double result;
   if(x1>=x2)      result=x2;
   else            result=x1; 
   return result;
}


void updateCurrent_villasenor(Domain *D,int nSpecies,int iteration);
void updateCurrent_Split_umeda(Domain *D,int nSpecies,int iteration);
void updateCurrent_Yee_umeda(Domain *D,int nSpecies,int iteration);
void calculaionRally(double *xr,double *yr,double rr,double x1,double x2,double y1,double y2,double r1,double r2,int iteration);
void compensateCurrent_Split_umeda(Domain *D,int iteration);
void MPI_Transfer6F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int share);

void updateCurrent(Domain D,int iteration)
{
  	int nxSub,nySub,i,j,n,m,iend,numMode;
  	int nSpecies;

  	nxSub=D.nxSub+5;  nySub=D.nySub+5;
  	nSpecies=D.nSpecies;
  	if(D.boostOn==ON && D.boostIon==OFF)
  	  	nSpecies=1;    

  	iend=D.iend;
  	numMode=D.numMode;    

	switch(D.currentType)  {
	case 1 :    

		if(D.fieldType==Split) updateCurrent_Split_umeda(&D,nSpecies,iteration);
    	else if(D.fieldType==Yee) updateCurrent_Yee_umeda(&D,nSpecies,iteration);
	 	else ;
		break;
	}



   D.shareF[0]=D.JzR;
   D.shareF[1]=D.JzI;
   D.shareF[2]=D.JrR;
   D.shareF[3]=D.JrI;
   D.shareF[4]=D.JpR;
   D.shareF[5]=D.JpI;	
   MPI_TransferJNew_Xplus(&D,6,D.nySub+5,3);
   MPI_TransferJNew_Xminus(&D,6,D.nySub+5,3);
	if(D.Period==ON) { MPI_TransferJNew_Period_X(&D,6,D.nySub+5,3); } else ;  
   D.JzR=D.shareF[0];
   D.JzI=D.shareF[1];
   D.JrR=D.shareF[2];
   D.JrI=D.shareF[3];
   D.JpR=D.shareF[4];
   D.JpI=D.shareF[5];

   D.shareF[0]=D.JzR;
   D.shareF[1]=D.JzI;
   D.shareF[2]=D.JrR;
   D.shareF[3]=D.JrI;
   D.shareF[4]=D.JpR;
   D.shareF[5]=D.JpI;	
   MPI_TransferFNew_Xplus(&D,6,D.nySub+5,3);
   MPI_TransferFNew_Xminus(&D,6,D.nySub+5,3);
	if(D.Period==ON) { MPI_TransferFNew_Period_X(&D,6,D.nySub+5,3); } else ;  
   D.JzR=D.shareF[0];
   D.JzI=D.shareF[1];
   D.JrR=D.shareF[2];
   D.JrI=D.shareF[3];
   D.JpR=D.shareF[4];
   D.JpI=D.shareF[5];
    
}



void updateCurrent_Split_umeda(Domain *D,int nSpecies,int iteration)
{
   int i,j,m,s,i1,i2,j1,j2,ii,jj,numMode,n,JJ;
   int istart,iend,jstart,jend,minRSub,dataI[2],dataJ[2],gridI,gridJ;
   int nxSub,nySub,jjstart;
   double oldZ,oldR,weight,gamma;
   double x1,x2,y1,y2,z1,z2,r1,r2,alpha;
   double Fz,Fr,factor[2],factM;
   double Wz[2],Wr[2],Wz1[2],Wz2[2],Wr1[2],Wr2[2];
   double tmpZ[2],tmpR[2],tmpP[2][2];
   double vp,xr,yr,zr,rr,xc,yc,zc,rc,vx,vy,v;
   double coss[D->numMode],sins[D->numMode];
   double cos1[D->numMode],sin1[D->numMode];
   double cos2[D->numMode],sin2[D->numMode];
   double dz,dr,dt,inverDt,drBydt,dzBydt,dtOverdr,dtOverdz;
	double posX1,posY1,posZ1,posR1;
	double posX2,posY2,posZ2,posR2;
	double dataX[3],dataY[3],dataZ[3],dataR[3];

   ptclList *p;
   LoadList *LL;
   Particle **particle;
   particle=D->particle;

   double coeff[nSpecies];

   istart=D->istart;    iend=D->iend;
   jstart=D->jstart;    jend=D->jend;
   nxSub=D->nxSub;      nySub=D->nySub;
   minRSub=D->minYSub;

   dt=D->dt; dz=D->dz; dr=D->dr; inverDt=1.0/dt;
   numMode=D->numMode;
   drBydt=dr/dt;  dzBydt=dz/dt;
	dtOverdz=dt/dz; dtOverdr=dt/dr;   

   s=0;
   LL=D->loadList;
   while(LL->next) {
      coeff[s]=LL->density/LL->criticalDensity;
      LL=LL->next;
      s++;
   }
   int myrank,nTasks;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   //initialize J
	//if(D->fieldIonization==OFF) {
		for(m=0; m<numMode; m++)
         for(i=0; i<nxSub+5; i++)
       	   for(j=0; j<nySub+5; j++)   {
               D->JzR[m][i][j]=0.0;
				   D->JrR[m][i][j]=0.0;
				   D->JpR[m][i][j]=0.0;
				   D->JzI[m][i][j]=0.0;
				   D->JrI[m][i][j]=0.0;
				   D->JpI[m][i][j]=0.0;
			  }
	//} else ;

   alpha=2.0; 
   for(i=istart; i<iend; i++)
      for(j=jstart+2; j<jend; j++)
         for(s=0; s<nSpecies; s++)
         {
            p=particle[i][j].head[s]->pt;
            while(p)
            {
               weight=p->weight*p->charge;
               gamma=sqrt(1.0+p->pz*p->pz+p->px*p->px+p->py*p->py);

               //------------------ Jr, Jp -------------------
               x1=p->oldX; y1=p->oldY; r1=sqrt(x1*x1+y1*y1);
               x2=p->x;    y2=p->y;    r2=sqrt(x2*x2+y2*y2);
               z2=p->z+i;  z1=p->oldZ;
               i1=(int)z1;           j1=(int)r1;
               i2=(int)z2;           j2=(int)r2;
               //rally calculation
               if(i1==i2)       zr=0.5*(z1+z2);
               else             zr=maximum(i1*1.0,i2*1.0);
               if(j1==j2)       rr=0.5*(r1+r2);
               else             rr=maximum(j1*1.0,j2*1.0);
               xr=x1; yr=y1;
               calculaionRally(&xr,&yr,rr,x1,x2,y1,y2,r1,r2,iteration);
		         dataX[0]=x1; dataY[0]=y1; dataR[0]=r1; dataZ[0]=z1;
		         dataX[1]=xr; dataY[1]=yr; dataR[1]=rr; dataZ[1]=zr;
		         dataX[2]=x2; dataY[2]=y2; dataR[2]=r2; dataZ[2]=z2;
		         dataI[0]=i1; dataJ[0]=j1;
		         dataI[1]=i2; dataJ[1]=j2;
		         for(n=0; n<2; n++)  {
		   	      gridI=dataI[n]; gridJ=dataJ[n];
		   	      posR1=dataR[n]; posR2=dataR[n+1];
		   	      posZ1=dataZ[n]; posZ2=dataZ[n+1];
		   	      posX1=dataX[n]; posX2=dataX[n+1];
		   	      posY1=dataY[n]; posY2=dataY[n+1];
		   	      rc=0.5*(posR1+posR2);     zc=0.5*(posZ1+posZ2);
		   	      //Wr[0]=((gridJ+1.0)*(gridJ+1.0)-rc*rc)/(2.0*gridJ+1.0);
		   	      Wr[0]=(gridJ+rc)*(1.0-(rc-gridJ))/(2.0*rc); 						
						Wr[1]=1.0-Wr[0];
		   	      Wz[1]=zc-gridI;        Wz[0]=1.0-Wz[1];
		   	      xc=posX1; yc=posY1;
		   	      calculaionRally(&xc,&yc,rc,posX1,posX2,posY1,posY2,posR1,posR2,iteration);
		   	      coss[1]=xc/rc; sins[1]=yc/rc;				
		   	      for(m=2; m<numMode; m++) {
		   	        coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
		   	        sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
		   	      }
		   	      factor[0]=weight*coeff[s]/(2.0*gridJ+1.0);
		   	      Fr=(posR2-posR1)*drBydt;
		   	      tmpR[0]=Fr*Wz[0]*factor[0];
		   	      tmpR[1]=Fr*Wz[1]*factor[0];
		   	      for(ii=0; ii<2; ii++)
		   	        D->JrR[0][ii+gridI][gridJ+jstart]+=tmpR[ii];
		   	      for(m=1; m<numMode; m++)  
		   	        for(ii=0; ii<2; ii++) {
		   	      		D->JrR[m][ii+gridI][gridJ+jstart]+=tmpR[ii]*coss[m]*alpha;
		   	      		D->JrI[m][ii+gridI][gridJ+jstart]-=tmpR[ii]*sins[m]*alpha;
		   	        }
		   	      //Jp
		   	      factor[0]=weight*coeff[s]/(2.0*gridJ);
		   	      factor[1]=weight*coeff[s]/(2.0*(gridJ+1));
		   	      vp=coss[1]*(posY2-posY1)*drBydt-sins[1]*(posX2-posX1)*drBydt;
		   	      tmpP[0][0]=Wz[0]*Wr[0]*factor[0]*vp;
		   	      tmpP[1][0]=Wz[1]*Wr[0]*factor[0]*vp;
		   	      tmpP[0][1]=Wz[0]*Wr[1]*factor[1]*vp;
		   	      tmpP[1][1]=Wz[1]*Wr[1]*factor[1]*vp;
		   	      for(ii=0; ii<2; ii++)
		   	        for(jj=0; jj<2; jj++)
		   	      		D->JpR[0][gridI+ii][gridJ+jj+jstart]+=tmpP[ii][jj];
   	   	   	    if(D->currentCons==Lifschitz) {
		   	            for(m=1; m<numMode; m++)  
		   	      		    for(ii=0; ii<2; ii++) 
		   	      			    for(jj=0; jj<2; jj++) {
		   	      				    D->JpR[m][gridI+ii][jj+gridJ+jstart]+=tmpP[ii][jj]*coss[m]*alpha;
		   	      				    D->JpI[m][gridI+ii][jj+gridJ+jstart]-=tmpP[ii][jj]*sins[m]*alpha;
		   	      			    }
		   	          } else if(D->currentCons==Davidson) {
		   	            cos1[1]=posX1/posR1; sin1[1]=posY1/posR1;
		   	            for(m=2; m<numMode; m++) {
		   	      		    cos1[m]=cos1[m-1]*cos1[1]-sin1[m-1]*sin1[1];
		   	      		    sin1[m]=sin1[m-1]*cos1[1]+cos1[m-1]*sin1[1];
		   	            }
		   	            cos2[1]=posX2/posR2; sin2[1]=posY2/posR2;
		   	            for(m=2; m<numMode; m++) {
		   	      	    	cos2[m]=cos2[m-1]*cos2[1]-sin2[m-1]*sin2[1];
		   	      	    	sin2[m]=sin2[m-1]*cos2[1]+cos2[m-1]*sin2[1];
		   	            }
		   	            Wz1[1]=posZ1-gridI; Wz1[0]=1.0-Wz1[1];
		   	            Wz2[1]=posZ2-gridI; Wz2[0]=1.0-Wz2[1];
		   	            Wr1[0]=(gridJ+posR1)*(1.0-(posR1-gridJ))/(2.0*posR1); Wr1[1]=1.0-Wr1[0];
		   	            Wr2[0]=(gridJ+posR2)*(1.0-(posR2-gridJ))/(2.0*posR2); Wr2[1]=1.0-Wr2[0];
   	   	   	      for(m=1; m<numMode; m++)  
		   	      	    	for(ii=0; ii<2; ii++) 
		   	      	    		for(jj=0; jj<2; jj++)  {
		   	      	    		  factM=alpha*(gridJ+jj)/(m*1.0)*drBydt*factor[jj];
		   	      	    		  D->JpR[m][ii+gridI][jj+gridJ+jstart]+=factM*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
		   	      	    		  D->JpI[m][ii+gridI][jj+gridJ+jstart]+=factM*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
		   	      	    		}
		   	          }
		         }

               //------------------- Jz ----------------------           
               x1=0.5*(p->oldX+p->x);
               y1=0.5*(p->oldY+p->y);
               r1=sqrt(x1*x1+y1*y1);
               z1=0.5*(p->oldZ+p->z+i);
               x2=x1+p->px/gamma*dtOverdr;
               y2=y1+p->py/gamma*dtOverdr;
		         r2=sqrt(x2*x2+y2*y2);
		         z2=z1+p->pz/gamma*dtOverdz;
               i1=(int)z1;           j1=(int)r1;
               i2=(int)z2;           j2=(int)r2;
               //rally calculation
               if(i1==i2)       zr=0.5*(z1+z2);
               else             zr=maximum(i1*1.0,i2*1.0);
               if(j1==j2)       rr=0.5*(r1+r2);
               else             rr=maximum(j1*1.0,j2*1.0);
               xr=x1; yr=y1;
               calculaionRally(&xr,&yr,rr,x1,x2,y1,y2,r1,r2,iteration);
		         dataX[0]=x1; dataY[0]=y1; dataR[0]=r1; dataZ[0]=z1;
		         dataX[1]=xr; dataY[1]=yr; dataR[1]=rr; dataZ[1]=zr;
		         dataX[2]=x2; dataY[2]=y2; dataR[2]=r2; dataZ[2]=z2;
		         dataI[0]=i1; dataJ[0]=j1;
		         dataI[1]=i2; dataJ[1]=j2;
		         for(n=0; n<2; n++)  {
		         	gridI=dataI[n]; gridJ=dataJ[n];
		         	posR1=dataR[n]; posR2=dataR[n+1];
		         	posZ1=dataZ[n]; posZ2=dataZ[n+1];
		         	posX1=dataX[n]; posX2=dataX[n+1];
		         	posY1=dataY[n]; posY2=dataY[n+1];
		         	rc=0.5*(posR1+posR2);     zc=0.5*(posZ1+posZ2);
   	   	      //Wr[0]=((gridJ+1.0)*(gridJ+1.0)-rc*rc)/(2.0*gridJ+1.0);
   	   	      Wr[0]=(gridJ+rc)*(1-(rc-gridJ))/(2.0*rc); 
						Wr[1]=1.0-Wr[0];
		         	Wz[1]=zc-gridI;        Wz[0]=1.0-Wz[1];
   	   	      xc=posX1; yc=posY1;
		         	calculaionRally(&xc,&yc,rc,posX1,posX2,posY1,posY2,posR1,posR2,iteration);
		         	coss[1]=xc/rc; sins[1]=yc/rc;				
		         	for(m=2; m<numMode; m++) {
		         	  coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
		         	  sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
		         	}
   	   	      factor[0]=weight*coeff[s]/(2.0*gridJ);
   	   	      factor[1]=weight*coeff[s]/(2.0*(gridJ+1));
   	            Fz=(posZ2-posZ1)*dzBydt;
		         	tmpZ[0]=Fz*Wr[0]*factor[0];
		         	tmpZ[1]=Fz*Wr[1]*factor[1];
		         	for(jj=0; jj<2; jj++)
		         	   D->JzR[0][gridI][jj+gridJ+jstart]+=tmpZ[jj];
		         	for(m=1; m<numMode; m++) 
		         	   for(jj=0; jj<2; jj++) {
		         			D->JzR[m][gridI][jj+gridJ+jstart]+=tmpZ[jj]*coss[m]*alpha;
		         			D->JzI[m][gridI][jj+gridJ+jstart]-=tmpZ[jj]*sins[m]*alpha;
		         	   }
		         }            
               p=p->next;
            }    //End of while(p)
         }    //End of for(s)     

	 
//-------------- for Axis -------------------------
   for(i=istart; i<iend; i++)
      for(j=jstart; j<jstart+2; j++)
         for(s=0; s<nSpecies; s++)
         {
            p=particle[i][j].head[s]->pt;
            while(p)
            {
               weight=p->weight*p->charge;
               gamma=sqrt(1.0+p->pz*p->pz+p->px*p->px+p->py*p->py);

               //------------------ Jr, Jp -------------------
               x1=p->oldX; y1=p->oldY; r1=sqrt(x1*x1+y1*y1);
               x2=p->x;    y2=p->y;    r2=sqrt(x2*x2+y2*y2);
               z2=p->z+i;  z1=p->oldZ;

               i1=(int)z1;           j1=(int)r1;
               i2=(int)z2;           j2=(int)r2;

               //rally calculation
               if(i1==i2)       zr=0.5*(z1+z2);
               else             zr=maximum(i1*1.0,i2*1.0);
               if(j1==j2)       rr=0.5*(r1+r2);
               else             rr=maximum(j1*1.0,j2*1.0);

				   xr=x1; yr=y1;
               calculaionRally(&xr,&yr,rr,x1,x2,y1,y2,r1,r2,iteration);

				   dataX[0]=x1; dataY[0]=y1; dataR[0]=r1; dataZ[0]=z1;
				   dataX[1]=xr; dataY[1]=yr; dataR[1]=rr; dataZ[1]=zr;
				   dataX[2]=x2; dataY[2]=y2; dataR[2]=r2; dataZ[2]=z2;
				   dataI[0]=i1; dataJ[0]=j1;
				   dataI[1]=i2; dataJ[1]=j2;

				   for(n=0; n<2; n++)  {
				   	gridI=dataI[n]; gridJ=dataJ[n];
				   	posR1=dataR[n]; posR2=dataR[n+1];
				   	posZ1=dataZ[n]; posZ2=dataZ[n+1];
				   	posX1=dataX[n]; posX2=dataX[n+1];
				   	posY1=dataY[n]; posY2=dataY[n+1];
				   	rc=0.5*(posR1+posR2);     zc=0.5*(posZ1+posZ2);
                  
						if(rc<0.5)      Wr[0]=(rc*rc-rc+0.5)/(2.0*rc*rc+0.5);
						else if(rc>=1)  Wr[0]=(gridJ+rc)*(1-(rc-gridJ))/(2.0*rc); 
						else            Wr[0]=(1.0-rc)*0.5;
						Wr[1]=1.0-Wr[0];
				   	Wz[1]=zc-gridI;        Wz[0]=1.0-Wz[1];

				   	xc=posX1; yc=posY1;
				   	calculaionRally(&xc,&yc,rc,posX1,posX2,posY1,posY2,posR1,posR2,iteration);

				   	coss[1]=xc/rc; sins[1]=yc/rc;				
				   	for(m=2; m<numMode; m++) {
				   	  coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
				   	  sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
				   	}


				   	factor[0]=weight*coeff[s]/(2.0*gridJ+1);
				   	Fr=(posR2-posR1)*drBydt;
				   	tmpR[0]=Fr*Wz[0]*factor[0];
				   	tmpR[1]=Fr*Wz[1]*factor[0];
				   	for(ii=0; ii<2; ii++)
				   	  D->JrR[0][ii+gridI][gridJ+jstart]+=tmpR[ii];
				   	for(m=1; m<numMode; m++)  
				   	  for(ii=0; ii<2; ii++) {
				   			D->JrR[m][ii+gridI][gridJ+jstart]+=tmpR[ii]*coss[m]*alpha;
				   			D->JrI[m][ii+gridI][gridJ+jstart]-=tmpR[ii]*sins[m]*alpha;
				   	  }

				   	//Jp
				   	if(gridJ==0) {
						  factor[0]=weight*coeff[s]/(0.25);
						  factor[1]=weight*coeff[s]/(2.0);
				   	  jjstart=1;
				   	} else {
				   	  factor[0]=weight*coeff[s]/(2.0*gridJ);
				   	  factor[1]=weight*coeff[s]/(2.0*(gridJ+1));
				   	  jjstart=0;
				   	}
				   	vp=coss[1]*(posY2-posY1)*drBydt-sins[1]*(posX2-posX1)*drBydt;
				   	tmpP[0][0]=Wz[0]*Wr[0]*factor[0]*vp;
				   	tmpP[1][0]=Wz[1]*Wr[0]*factor[0]*vp;
				   	tmpP[0][1]=Wz[0]*Wr[1]*factor[1]*vp;
				   	tmpP[1][1]=Wz[1]*Wr[1]*factor[1]*vp;
				   	for(ii=0; ii<2; ii++)
				   	   for(jj=jjstart; jj<2; jj++)
				   			D->JpR[0][gridI+ii][gridJ+jj+jstart]+=tmpP[ii][jj];
				   	if(D->currentCons==Lifschitz) {
				   	   m=1;  
				   			for(ii=0; ii<2; ii++) 
				   				for(jj=0; jj<2; jj++) {
				   					D->JpR[m][gridI+ii][jj+gridJ+jstart]+=tmpP[ii][jj]*coss[m]*alpha;
				   					D->JpI[m][gridI+ii][jj+gridJ+jstart]-=tmpP[ii][jj]*sins[m]*alpha;
				   				}
				   	  for(m=2; m<numMode; m++)  
				   			for(ii=0; ii<2; ii++) 
				   				for(jj=jjstart; jj<2; jj++) {
				   					D->JpR[m][gridI+ii][jj+gridJ+jstart]+=tmpP[ii][jj]*coss[m]*alpha;
				   					D->JpI[m][gridI+ii][jj+gridJ+jstart]-=tmpP[ii][jj]*sins[m]*alpha;
				   				}
				   	} else if(D->currentCons==Davidson) {
				   	   cos1[1]=posX1/posR1; sin1[1]=posY1/posR1;
				   	   for(m=2; m<numMode; m++) {
				   			cos1[m]=cos1[m-1]*cos1[1]-sin1[m-1]*sin1[1];
				   			sin1[m]=sin1[m-1]*cos1[1]+cos1[m-1]*sin1[1];
				   	   }
				   	   cos2[1]=posX2/posR2; sin2[1]=posY2/posR2;
				   	   for(m=2; m<numMode; m++) {
				   			cos2[m]=cos2[m-1]*cos2[1]-sin2[m-1]*sin2[1];
				   			sin2[m]=sin2[m-1]*cos2[1]+cos2[m-1]*sin2[1];
				   	   }
				   	   Wz1[1]=posZ1-gridI; Wz1[0]=1.0-Wz1[1];
				   	   Wz2[1]=posZ2-gridI; Wz2[0]=1.0-Wz2[1];
						   if(posR1<0.5)      Wr1[0]=(posR1*posR1-posR1+0.5)/(2.0*posR1*posR1+0.5);
						   else if(posR1>=1)  Wr1[0]=(gridJ+posR1)*(1-(posR1-gridJ))/(2.0*posR1); 
						   else               Wr1[0]=(1.0-posR1)*0.5;
                     Wr1[1]=1.0-Wr1[0];
						   
                     if(posR2<0.5)      Wr2[0]=(posR2*posR2-posR2+0.5)/(2.0*posR2*posR2+0.5);
						   else if(posR2>=1)  Wr2[0]=(gridJ+posR2)*(1-(posR2-gridJ))/(2.0*posR2); 
						   else               Wr2[0]=(1.0-posR2)*0.5;
                     Wr2[1]=1.0-Wr2[0];

                     m=1; 
                     for(ii=0; ii<2; ii++) 
                       for(jj=0; jj<2; jj++)  {
                     	  factM=alpha*(gridJ+jj)/(m*1.0)*drBydt;
                          D->JpR[m][ii+gridI][jj+gridJ+jstart]+=factM*factor[jj]*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
                     	  D->JpI[m][ii+gridI][jj+gridJ+jstart]+=factM*factor[jj]*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
                       }
						   for(m=2; m<numMode; m++)  
								for(ii=0; ii<2; ii++) 
									for(jj=jjstart; jj<2; jj++)  {
										factM=alpha*(gridJ+jj)/(m*1.0)*drBydt*factor[jj];
										D->JpR[m][ii+gridI][jj+gridJ+jstart]+=factM*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
										D->JpI[m][ii+gridI][jj+gridJ+jstart]+=factM*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
									}
					   }
				   }

               //------------------ Jz -------------------
               x1=0.5*(p->oldX+p->x);
               y1=0.5*(p->oldY+p->y);
               r1=sqrt(x1*x1+y1*y1);
               z1=0.5*(p->oldZ+p->z+i);
               x2=x1+p->px/gamma*dtOverdr;
               y2=y1+p->py/gamma*dtOverdr;
		         r2=sqrt(x2*x2+y2*y2);
		         z2=z1+p->pz/gamma*dtOverdz;

               i1=(int)z1;           j1=(int)r1;
               i2=(int)z2;           j2=(int)r2;

               //rally calculation
               if(i1==i2)       zr=0.5*(z1+z2);
               else             zr=maximum(i1*1.0,i2*1.0);
               if(j1==j2)       rr=0.5*(r1+r2);
               else             rr=maximum(j1*1.0,j2*1.0);

				   xr=x1; yr=y1;
               calculaionRally(&xr,&yr,rr,x1,x2,y1,y2,r1,r2,iteration);

				   dataX[0]=x1; dataY[0]=y1; dataR[0]=r1; dataZ[0]=z1;
				   dataX[1]=xr; dataY[1]=yr; dataR[1]=rr; dataZ[1]=zr;
				   dataX[2]=x2; dataY[2]=y2; dataR[2]=r2; dataZ[2]=z2;
				   dataI[0]=i1; dataJ[0]=j1;
				   dataI[1]=i2; dataJ[1]=j2;

				   for(n=0; n<2; n++)  {
					  gridI=dataI[n]; gridJ=dataJ[n];
					  posR1=dataR[n]; posR2=dataR[n+1];
					  posZ1=dataZ[n]; posZ2=dataZ[n+1];
					  posX1=dataX[n]; posX2=dataX[n+1];
					  posY1=dataY[n]; posY2=dataY[n+1];
					  rc=0.5*(posR1+posR2);     zc=0.5*(posZ1+posZ2);
					
                 //Wr[0]=((gridJ+1.0)*(gridJ+1.0)-rc*rc)/(2.0*gridJ+1.0); Wr[1]=1.0-Wr[0];
					  if(rc<0.5)      Wr[0]=(rc*rc-rc+0.5)/(2.0*rc*rc+0.5);
					  else if(rc>=1)  Wr[0]=(gridJ+rc)*(1-(rc-gridJ))/(2.0*rc); 
					  else            Wr[0]=(1.0-rc)*0.5;
					  Wr[1]=1.0-Wr[0];
					  Wz[1]=zc-gridI;        Wz[0]=1.0-Wz[1];
					
                 xc=posX1; yc=posY1;
					  calculaionRally(&xc,&yc,rc,posX1,posX2,posY1,posY2,posR1,posR2,iteration);
					
                 coss[1]=xc/rc; sins[1]=yc/rc;				
					  for(m=2; m<numMode; m++) {
					    coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
					    sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
					  }
				     if(gridJ==0) {
					    factor[0]=weight*coeff[s]/(0.25);
					    factor[1]=weight*coeff[s]/(2.0);
				       jjstart=1;
				     } else {
				       factor[0]=weight*coeff[s]/(2.0*gridJ);
				       factor[1]=weight*coeff[s]/(2.0*(gridJ+1));
				       jjstart=0;
				     }
	            
                 Fz=(posZ2-posZ1)*dzBydt;
					  tmpZ[0]=Fz*Wr[0]*factor[0];
					  tmpZ[1]=Fz*Wr[1]*factor[1];
					  for(jj=0; jj<2; jj++)
					    D->JzR[0][gridI][jj+gridJ+jstart]+=tmpZ[jj];
					  for(m=1; m<numMode; m++) 
					    for(jj=jjstart; jj<2; jj++) {
							D->JzR[m][gridI][jj+gridJ+jstart]+=tmpZ[jj]*coss[m]*alpha;
							D->JzI[m][gridI][jj+gridJ+jstart]-=tmpZ[jj]*sins[m]*alpha;
					    }
				   }            
            
               p=p->next;
            }    //End of while(p)

      }    //End of for(s)     

}


void updateCurrent_Yee_umeda(Domain *D,int nSpecies,int iteration)
{
  int i,j,m,s,i1,i2,j1,j2,ii,jj,numMode,n,JJ;
  int istart,iend,jstart,jend,minRSub,dataI[2],dataJ[2],gridI,gridJ;
  int nxSub,nySub,jjstart;
  double oldZ,oldR,weight,gamma;
  double x1,x2,y1,y2,z1,z2,r1,r2,alpha;
  double Fz,Fr,factor[2],factM;
  double Wz[2],Wr[2],Wz1[2],Wz2[2],Wr1[2],Wr2[2];
  double tmpZ[2],tmpR[2],tmpP[2][2];
  double vp,xr,yr,zr,rr,xc,yc,zc,rc,vx,vy,v;
  double coss[D->numMode],sins[D->numMode];
  double cos1[D->numMode],sin1[D->numMode];
  double cos2[D->numMode],sin2[D->numMode];
  double dz,dr,dt,inverDt,drBydt,dzBydt;
	double posX1,posY1,posZ1,posR1;
	double posX2,posY2,posZ2,posR2;
	double dataX[3],dataY[3],dataZ[3],dataR[3];

  ptclList *p;
  LoadList *LL;
  Particle **particle;
  particle=D->particle;

  double coeff[nSpecies];

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  nxSub=D->nxSub;      nySub=D->nySub;
  minRSub=D->minYSub;

  dt=D->dt; dz=D->dz; dr=D->dr; inverDt=1.0/dt;
  numMode=D->numMode;
  drBydt=dr/dt;  dzBydt=dz/dt;

  s=0;
  LL=D->loadList;
  while(LL->next) {
     coeff[s]=LL->density/LL->criticalDensity;
     LL=LL->next;
     s++;
  }
  int myrank,nTasks;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  //initialize J
	if(D->fieldIonization==OFF) {
		for(m=0; m<numMode; m++)
         for(i=0; i<nxSub+5; i++)
       	   for(j=0; j<nySub+5; j++)   {
               D->JzR[m][i][j]=0.0;
				   D->JrR[m][i][j]=0.0;
				   D->JpR[m][i][j]=0.0;
				   D->JzI[m][i][j]=0.0;
				   D->JrI[m][i][j]=0.0;
				   D->JpI[m][i][j]=0.0;
			   }
	} else ;

  alpha=2.0; 
  for(i=istart; i<iend; i++)
    for(j=jstart+2; j<jend; j++)
      for(s=0; s<nSpecies; s++)
      {
        p=particle[i][j].head[s]->pt;
        while(p)
        {
          weight=p->weight*p->charge;
          gamma=sqrt(1.0+p->pz*p->pz+p->px*p->px+p->py*p->py);

          x1=p->oldX; y1=p->oldY; r1=sqrt(x1*x1+y1*y1);
          x2=p->x;    y2=p->y;    r2=sqrt(x2*x2+y2*y2);
          z2=p->z+i;  z1=p->oldZ;

          i1=(int)z1;           j1=(int)r1;
          i2=(int)z2;           j2=(int)r2;

          //rally calculation
          if(i1==i2)       zr=0.5*(z1+z2);
          else             zr=maximum(i1*1.0,i2*1.0);
          if(j1==j2)       rr=0.5*(r1+r2);
          else             rr=maximum(j1*1.0,j2*1.0);
          xr=x1; yr=y1;
          calculaionRally(&xr,&yr,rr,x1,x2,y1,y2,r1,r2,iteration);

		    dataX[0]=x1; dataY[0]=y1; dataR[0]=r1; dataZ[0]=z1;
		    dataX[1]=xr; dataY[1]=yr; dataR[1]=rr; dataZ[1]=zr;
		    dataX[2]=x2; dataY[2]=y2; dataR[2]=r2; dataZ[2]=z2;
		    dataI[0]=i1; dataJ[0]=j1;
		    dataI[1]=i2; dataJ[1]=j2;

		    for(n=0; n<2; n++)  {
			    gridI=dataI[n]; gridJ=dataJ[n];
			    posR1=dataR[n]; posR2=dataR[n+1];
			    posZ1=dataZ[n]; posZ2=dataZ[n+1];
			    posX1=dataX[n]; posX2=dataX[n+1];
			    posY1=dataY[n]; posY2=dataY[n+1];
			    rc=0.5*(posR1+posR2);     zc=0.5*(posZ1+posZ2);
             Wr[0]=(gridJ+rc)*(1.0-(rc-gridJ))/(2.0*rc);
             Wr[1]=1.0-Wr[0];
   		    Wz[1]=zc-gridI;        Wz[0]=1.0-Wz[1];

			    xc=posX1; yc=posY1;
			    calculaionRally(&xc,&yc,rc,posX1,posX2,posY1,posY2,posR1,posR2,iteration);
			    coss[1]=xc/rc; sins[1]=yc/rc;				
			    for(m=2; m<numMode; m++) {
			       coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
			       sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
	          }
             factor[0]=weight*coeff[s]/(2.0*gridJ+1.0);
			    Fr=(posR2-posR1)*drBydt;
			    tmpR[0]=Fr*Wz[0]*factor[0];
			    tmpR[1]=Fr*Wz[1]*factor[0];
			    for(ii=0; ii<2; ii++)
				    D->JrR[0][ii+gridI][gridJ+jstart]+=tmpR[ii];
				 for(m=1; m<numMode; m++)  
				    for(ii=0; ii<2; ii++) {
				       D->JrR[m][ii+gridI][gridJ+jstart]+=tmpR[ii]*coss[m]*alpha;
				    	 D->JrI[m][ii+gridI][gridJ+jstart]-=tmpR[ii]*sins[m]*alpha;
				    }

             factor[0]=weight*coeff[s]/(2.0*gridJ);
             factor[1]=weight*coeff[s]/(2.0*(gridJ+1));
				 //Jp
				 vp=coss[1]*(posY2-posY1)*drBydt-sins[1]*(posX2-posX1)*drBydt;
				 tmpP[0][0]=Wz[0]*Wr[0]*factor[0]*vp;
			    tmpP[1][0]=Wz[1]*Wr[0]*factor[0]*vp;
			    tmpP[0][1]=Wz[0]*Wr[1]*factor[1]*vp;
			    tmpP[1][1]=Wz[1]*Wr[1]*factor[1]*vp;
			    for(ii=0; ii<2; ii++)
                for(jj=0; jj<2; jj++)
			          D->JpR[0][gridI+ii][gridJ+jj+jstart]+=tmpP[ii][jj];
  				 if(D->currentCons==Lifschitz) {
				    for(m=1; m<numMode; m++)  
				       for(ii=0; ii<2; ii++) 
				          for(jj=0; jj<2; jj++) {
				    	       D->JpR[m][gridI+ii][jj+gridJ+jstart]+=tmpP[ii][jj]*coss[m]*alpha;
				    		    D->JpI[m][gridI+ii][jj+gridJ+jstart]-=tmpP[ii][jj]*sins[m]*alpha;
				    		 }
				 } else if(D->currentCons==Davidson) {
				    cos1[1]=posX1/posR1; sin1[1]=posY1/posR1;
				    for(m=2; m<numMode; m++) {
				       cos1[m]=cos1[m-1]*cos1[1]-sin1[m-1]*sin1[1];
				    	 sin1[m]=sin1[m-1]*cos1[1]+cos1[m-1]*sin1[1];
				    }
				    cos2[1]=posX2/posR2; sin2[1]=posY2/posR2;
				    for(m=2; m<numMode; m++) {
				       cos2[m]=cos2[m-1]*cos2[1]-sin2[m-1]*sin2[1];
				    	 sin2[m]=sin2[m-1]*cos2[1]+cos2[m-1]*sin2[1];
				    }
				    Wz1[1]=posZ1-gridI; Wz1[0]=1.0-Wz1[1];
				    Wz2[1]=posZ2-gridI; Wz2[0]=1.0-Wz2[1];
                Wr1[0]=(gridJ+posR1)*(1.0-(posR1-gridJ))/(2.0*posR1); Wr1[1]=1.0-Wr1[0];
                Wr2[0]=(gridJ+posR2)*(1.0-(posR2-gridJ))/(2.0*posR2); Wr2[1]=1.0-Wr2[0]; 
                for(m=1; m<numMode; m++)  
			          for(ii=0; ii<2; ii++) 
				    	    for(jj=0; jj<2; jj++)  {
				    	       factM=alpha*(gridJ+jj)/(m*1.0)*drBydt*factor[jj];
				    	    	 D->JpR[m][ii+gridI][jj+gridJ+jstart]+=factM*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
				    	    	 D->JpI[m][ii+gridI][jj+gridJ+jstart]+=factM*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
				    	    }
				 }

             //Jz
  	          Fz=(posZ2-posZ1)*dzBydt;
			    tmpZ[0]=Fz*Wr[0]*factor[0];
			    tmpZ[1]=Fz*Wr[1]*factor[1];
             for(jj=0; jj<2; jj++)
			       D->JzR[0][gridI][jj+gridJ+jstart]+=tmpZ[jj];
			    for(m=1; m<numMode; m++) 
			       for(jj=0; jj<2; jj++) {
			    	    D->JzR[m][gridI][jj+gridJ+jstart]+=tmpZ[jj]*coss[m]*alpha;
			    		 D->JzI[m][gridI][jj+gridJ+jstart]-=tmpZ[jj]*sins[m]*alpha;
			    	 }

			 }     

          p=p->next;
        }    //End of while(p)

      }    //End of for(s)     

	 
  //-------------- for Axis -------------------------
  for(i=istart; i<iend; i++)
     for(j=jstart; j<jstart+2; j++)
        for(s=0; s<nSpecies; s++)
        {
           p=particle[i][j].head[s]->pt;
           while(p)
           {
              weight=p->weight*p->charge;
              gamma=sqrt(1.0+p->pz*p->pz+p->px*p->px+p->py*p->py);

              //------------------ Jr, Jp -------------------
              x1=p->oldX; y1=p->oldY; r1=sqrt(x1*x1+y1*y1);
              x2=p->x;    y2=p->y;    r2=sqrt(x2*x2+y2*y2);
              z2=p->z+i;  z1=p->oldZ;

              i1=(int)z1;           j1=(int)r1;
              i2=(int)z2;           j2=(int)r2;

              //rally calculation
              if(i1==i2)       zr=0.5*(z1+z2);
              else             zr=maximum(i1*1.0,i2*1.0);
              if(j1==j2)       rr=0.5*(r1+r2);
              else             rr=maximum(j1*1.0,j2*1.0);

			     xr=x1; yr=y1;
              calculaionRally(&xr,&yr,rr,x1,x2,y1,y2,r1,r2,iteration);

			     dataX[0]=x1; dataY[0]=y1; dataR[0]=r1; dataZ[0]=z1;
			     dataX[1]=xr; dataY[1]=yr; dataR[1]=rr; dataZ[1]=zr;
			     dataX[2]=x2; dataY[2]=y2; dataR[2]=r2; dataZ[2]=z2;
			     dataI[0]=i1; dataJ[0]=j1;
			     dataI[1]=i2; dataJ[1]=j2;

			     for(n=0; n<2; n++)  {
			        gridI=dataI[n]; gridJ=dataJ[n];
			        posR1=dataR[n]; posR2=dataR[n+1];
			        posZ1=dataZ[n]; posZ2=dataZ[n+1];
			        posX1=dataX[n]; posX2=dataX[n+1];
			        posY1=dataY[n]; posY2=dataY[n+1];
			        rc=0.5*(posR1+posR2);     zc=0.5*(posZ1+posZ2);

                 if(rc<0.5)      Wr[0]=(rc*rc-rc+0.5)/(2.0*rc*rc+0.5);
                 else if(rc>=1)  Wr[0]=(gridJ+rc)*(1-(rc-gridJ))/(2.0*rc);
                 else            Wr[0]=(1.0-rc)*0.5;
                 Wr[1]=1.0-Wr[0];
			        Wz[1]=zc-gridI;        Wz[0]=1.0-Wz[1];
  			    
                 xc=posX1; yc=posY1;
				     calculaionRally(&xc,&yc,rc,posX1,posX2,posY1,posY2,posR1,posR2,iteration);
				    
                 coss[1]=xc/rc; sins[1]=yc/rc;				
				     for(m=2; m<numMode; m++) {
					     coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
					     sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
				     }

                 //Jr
                 factor[0]=weight*coeff[s]/(2.0*gridJ+1);
				     Fr=(posR2-posR1)*drBydt;
				     tmpR[0]=Fr*Wz[0]*factor[0];
				     tmpR[1]=Fr*Wz[1]*factor[0];
				     for(ii=0; ii<2; ii++)
					     D->JrR[0][ii+gridI][gridJ+jstart]+=tmpR[ii];
				     for(m=1; m<numMode; m++)  
					     for(ii=0; ii<2; ii++) {
						     D->JrR[m][ii+gridI][gridJ+jstart]+=tmpR[ii]*coss[m]*alpha;
							  D->JrI[m][ii+gridI][gridJ+jstart]-=tmpR[ii]*sins[m]*alpha;
					     }

                if(gridJ==0) {
                   factor[0]=weight*coeff[s]/(0.25);
                   factor[1]=weight*coeff[s]/(2.0);
                   jjstart=1;
                } else {
                   factor[0]=weight*coeff[s]/(2.0*gridJ);
                   factor[1]=weight*coeff[s]/(2.0*(gridJ+1));
                   jjstart=0;
                }
				    //Jp
				    vp=coss[1]*(posY2-posY1)*drBydt-sins[1]*(posX2-posX1)*drBydt;
				    tmpP[0][0]=Wz[0]*Wr[0]*factor[0]*vp;
				    tmpP[1][0]=Wz[1]*Wr[0]*factor[0]*vp;
				    tmpP[0][1]=Wz[0]*Wr[1]*factor[1]*vp;
				    tmpP[1][1]=Wz[1]*Wr[1]*factor[1]*vp;
				    for(ii=0; ii<2; ii++)
				       for(jj=jjstart; jj<2; jj++)
				    		 D->JpR[0][gridI+ii][gridJ+jj+jstart]+=tmpP[ii][jj];
				    if(D->currentCons==Lifschitz) {
				       m=1;  
				    	 for(ii=0; ii<2; ii++) 
				    	    for(jj=0; jj<2; jj++) {
				    		    D->JpR[m][gridI+ii][jj+gridJ+jstart]+=tmpP[ii][jj]*coss[m]*alpha;
				    			 D->JpI[m][gridI+ii][jj+gridJ+jstart]-=tmpP[ii][jj]*sins[m]*alpha;
				    		 }
				       for(m=2; m<numMode; m++)  
				    	    for(ii=0; ii<2; ii++) 
				    		    for(jj=jjstart; jj<2; jj++) {
				    			    D->JpR[m][gridI+ii][jj+gridJ+jstart]+=tmpP[ii][jj]*coss[m]*alpha;
				    				 D->JpI[m][gridI+ii][jj+gridJ+jstart]-=tmpP[ii][jj]*sins[m]*alpha;
				    			 }
				    } else if(D->currentCons==Davidson) {
				       cos1[1]=posX1/posR1; sin1[1]=posY1/posR1;
				       for(m=2; m<numMode; m++) {
				    	    cos1[m]=cos1[m-1]*cos1[1]-sin1[m-1]*sin1[1];
				    		 sin1[m]=sin1[m-1]*cos1[1]+cos1[m-1]*sin1[1];
				       }
				       cos2[1]=posX2/posR2; sin2[1]=posY2/posR2;
				       for(m=2; m<numMode; m++) {
				    		 cos2[m]=cos2[m-1]*cos2[1]-sin2[m-1]*sin2[1];
				    		 sin2[m]=sin2[m-1]*cos2[1]+cos2[m-1]*sin2[1];
				       }
				       Wz1[1]=posZ1-gridI; Wz1[0]=1.0-Wz1[1];
				       Wz2[1]=posZ2-gridI; Wz2[0]=1.0-Wz2[1];
                   if(posR1<0.5)      Wr1[0]=(posR1*posR1-posR1+0.5)/(2.0*posR1*posR1+0.5);
                   else if(posR1>=1)  Wr1[0]=(gridJ+posR1)*(1-(posR1-gridJ))/(2.0*posR1);
                   else               Wr1[0]=(1.0-posR1)*0.5;
                   Wr1[1]=1.0-Wr1[0];

                   if(posR2<0.5)      Wr2[0]=(posR2*posR2-posR2+0.5)/(2.0*posR2*posR2+0.5);
                   else if(posR2>=1)  Wr2[0]=(gridJ+posR2)*(1-(posR2-gridJ))/(2.0*posR2);
                   else               Wr2[0]=(1.0-posR2)*0.5;
                   Wr2[1]=1.0-Wr2[0];
						 
                   m=1;
						    for(ii=0; ii<2; ii++) 
							    for(jj=jjstart; jj<2; jj++)  {
								    factM=alpha*(gridJ+jj)/(m*1.0)*drBydt*factor[jj];
									 D->JpR[m][ii+gridI][jj+gridJ+jstart]+=factM*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
									 D->JpI[m][ii+gridI][jj+gridJ+jstart]+=factM*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
								 }
						 for(m=2; m<numMode; m++)  
						    for(ii=0; ii<2; ii++) 
							    for(jj=jjstart; jj<2; jj++)  {
								    factM=alpha*(gridJ+jj)/(m*1.0)*drBydt*factor[jj];
									 D->JpR[m][ii+gridI][jj+gridJ+jstart]+=factM*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
									 D->JpI[m][ii+gridI][jj+gridJ+jstart]+=factM*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
								 }
                }      
					  
                // Jz
	             Fz=(posZ2-posZ1)*dzBydt;
					 tmpZ[0]=Fz*Wr[0]*factor[0];
					 tmpZ[1]=Fz*Wr[1]*factor[1];
					 for(jj=0; jj<2; jj++)
					    D->JzR[0][gridI][jj+gridJ+jstart]+=tmpZ[jj];
					 for(m=1; m<numMode; m++) 
					    for(jj=jjstart; jj<2; jj++) {
						    D->JzR[m][gridI][jj+gridJ+jstart]+=tmpZ[jj]*coss[m]*alpha;
							 D->JzI[m][gridI][jj+gridJ+jstart]-=tmpZ[jj]*sins[m]*alpha;
						 }

				 }

             p=p->next;
          }    //End of while(p)

       }    //End of for(s)     

}



void calculaionRally(double *xr,double *yr,double rr,double x1,double x2,double y1,double y2,double r1,double r2,int iteration)
{
  int mode;
  double xc,yc,rc,cosTh,sinTh,alpha,signX,tmp;

  xc=0.5*(x1+x2);
  yc=0.5*(y1+y2);
  rc=sqrt(xc*xc+yc*yc);
  cosTh=xc/rc;
  sinTh=yc/rc;
  
  if(rc>0)  {
    *xr=cosTh*rr;
    *yr=sinTh*rr;
  } else {
		*xr=x1;
		*yr=y1;
	}

/*
   if(x1!=x2 && xc!=0.0)  {
     alpha=yc/xc;
//     sign=xc/fabs(xc);
     signX=copysign(1.0,xc);
	  tmp=sqrt(rr*rr/(1.0+alpha*alpha));
     *xr=signX*tmp;
     *yr=alpha*(*xr);
     mode=1;
   } else {
		*xr=x1;
		*yr=y1;
	}
*/
}



void findPhaseDif(int m,double  cosP1,double cosP2,double sinP1,double sinP2,double s1,double s2,double *A,double *B,double *X,double *Y,double *coss,double *sins)
{
  switch (m) {
    case 1 :
      *A=0.5*(sqrt((1.0+cosP2)*(1.0+cosP1))+s1*s2*sqrt((1.0-cosP2)*(1.0-cosP1)));
      *B=0.5*(s2*sqrt((1.0-cosP2)*(1.0+cosP1))-s1*sqrt((1.0+cosP2)*(1.0-cosP1)));
      *coss=*A;
      *sins=*B;
      break;
    case 2 :
      *X=cosP2*cosP1+sinP2*sinP1;
      *Y=sinP2*cosP1-cosP2*sinP1;
      *coss=*X;
      *sins=*Y;
      break;
    case 3 :
      *coss=(*A)*(*X)-(*B)*(*Y);
      *sins=(*B)*(*X)+(*A)*(*Y);
      break;
    case 4 :
      *coss=(*X)*(*X)-(*Y)*(*Y);
      *sins=2.0*(*X)*(*Y);
      break;
    default :
      printf("out of mode. Check number_mode.\n");
  }

}

void findPhaseAdd(int m,double  cosP1,double cosP2,double sinP1,double sinP2,double s1,double s2,double *A,double *B,double *X,double *Y,double *coss,double *sins)
{
  switch (m) {
    case 1 :
      *A=0.5*(sqrt((1.0+cosP2)*(1.0+cosP1))-s1*s2*sqrt((1.0-cosP2)*(1.0-cosP1)));
      *B=0.5*(s2*sqrt((1.0-cosP2)*(1.0+cosP1))+s1*sqrt((1.0+cosP2)*(1.0-cosP1)));
      *coss=*A;
      *sins=*B;
      break;
    case 2 :
      *X=cosP2*cosP1-sinP2*sinP1;
      *Y=sinP2*cosP1+cosP2*sinP1;
      *coss=*X;
      *sins=*Y;
      break;
    case 3 :
      *coss=(*A)*(*X)-(*B)*(*Y);
      *sins=(*B)*(*X)+(*A)*(*Y);
      break;
    case 4 :
      *coss=(*X)*(*X)-(*Y)*(*Y);
      *sins=2.0*(*X)*(*Y);
      break;
    default :
      printf("out of mode. Check number_mode.\n");
  }

}

double findR(double x1, double x2,double x3, double x4)
{
//  double minimum();
//  double maximum();
  double result,result1,result2,result3;

  result1=MIN(x1-0.5,x2-0.5);
  result2=MAX(x1-1.5,x2-1.5);
  result3=MAX(result2,(x3+x4)*0.5);
  result=MIN(result1,result3);

  return result;
}


int intmaximum(int x1,int x2)
{
   int result;

   if(x1>=x2)
      result=x1;
   else
      result=x2;
  
   return result;
}

int intminimum(int x1,int x2)
{
   int result;

   if(x1>=x2)
      result=x2;
   else
      result=x1;
  
   return result;
}

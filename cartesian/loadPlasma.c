#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_qrng.h>


double maxwellianVelocity(double temperature,double mass);
void random1D_sobol(double *x,gsl_qrng *q1);
void random2D_sobol(double *x,double *y,gsl_qrng *q2);
void random3D_sobol(double *x,double *y,double *z,gsl_qrng *q3);
void loadDefinedPlasma(Domain *D,LoadList *LL,int s);
double gaussian_dist(double sig);
void loadPolygonPlasma(Domain *D,LoadList *LL,int s,int iteration,int istart,int iend);

double applyFunctionX(int mode,double centerX,double x,double gaussCoefX,double polyCoefX)
{
  double result;

  switch (mode)  {
  case 0 :	//Costant
    result=1.0;
    break;
  case 1 :	//Gaussian
    result=exp(-(x-centerX)*(x-centerX)/gaussCoefX/gaussCoefX);
    break;
  case 2 :	//2nd polynomial
    result=1.0+polyCoefX*(x-centerX)*(x-centerX);
    break;
  }
  return result;
}

double applyFunctionYZ(int mode,double centerY,double y,double centerZ,double z,double gaussCoefYZ,double polyCoefYZ)
{
  double result;

  switch (mode)  {
  case 0 :	//Costant
    result=1.0;
    break;
  case 1 :	//Gaussian
    result=exp(-((y-centerY)*(y-centerY)+(z-centerZ)*(z-centerZ))/gaussCoefYZ/gaussCoefYZ);
    break;
  case 2 :	//2nd polynomial
    result=1.0+polyCoefYZ*((y-centerY)*(y-centerY)+(z-centerZ)*(z-centerZ));
    break;
  }
  return result;
}

void loadPlasma(Domain *D,LoadList *LL,int s,int istart,int iend,int iteration)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  switch(LL->type)  {
  case Polygon:
    loadPolygonPlasma(D,LL,s,istart,iend,iteration); 
    break;

  case (Defined-1)*3+3:
//    if(LL->minLoadTime<=iteration && iteration<=LL->maxLoadTime) {
    loadDefinedPlasma(D,LL,s); 
//    }
    MPI_Barrier(MPI_COMM_WORLD);
    break;

  case Beam:
    break;

  default:
    ;
  }
}


void loadDefinedPlasma(Domain *D,LoadList *LL,int s)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,cnt;
   int l,n,myrank;
   double dx,dy,dz,v1,v2,v3,charge,weightCoef;
   double x,y,z,xPos,yPos,zPos,coef2D,coef3D;
   Particle ***particle;
   particle=D->particle;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   ptclList *New;
   DefPtcl *p;

   dx=D->dx;   dy=D->dy;   dz=D->dz;

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   cnt=LL->numDefPtcls;
   charge=LL->charge;
   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;
   coef2D=0; coef3D=0;
   if(D->dimension>1) { coef2D=1; } else ;
   if(D->dimension>2) { coef3D=1; } else ;

   //position define      
   p=LL->def;
   while(p) {
     cnt=p->numDefPtcls;
     for(n=0; n<cnt; n++)
     {
       xPos=p->define[n];
       i=(int)(xPos/dx-D->minXSub+istart);
       if(i<iend && i>=istart)
       {
         yPos=p->define[n+cnt];
         j=(int)(yPos/dy-D->minYSub+jstart);
         if(jstart<=j && j<jend)
         {
           zPos=p->define[n+2*cnt];
           k=(int)(zPos/dz-D->minZSub+kstart);
           if(kstart<=k && k<kend)
           {
             x=xPos/dx-D->minXSub+istart-i;
             y=yPos/dy-D->minYSub+jstart-j;
             z=zPos/dz-D->minZSub+kstart-k;
             New = (ptclList *)malloc(sizeof(ptclList));
             New->next = particle[i][j][k].head[s]->pt;
             particle[i][j][k].head[s]->pt = New;

             New->x = x;          New->oldX=i+x;
             New->y = y*coef2D;   New->oldY=(j+y)*coef2D;
             New->z = z*coef3D;   New->oldZ=(k+z)*coef3D;

             New->E1=New->E2=New->E3=0.0;
             New->B1=New->B2=New->B3=0.0;
             v1=maxwellianVelocity(LL->temperature,LL->mass)/velocityC;
             v2=maxwellianVelocity(LL->temperature,LL->mass)/velocityC;
             v3=maxwellianVelocity(LL->temperature,LL->mass)/velocityC;
             New->p1=-D->gamma*D->beta+v1;
             New->p2=v2;
             New->p3=v3;
          	 New->p1Old1=New->p1; New->p2Old1=New->p2; New->p3Old1=New->p3;
             New->weight=1.0/LL->numberInCell*weightCoef;
             New->charge=charge;
             D->index+=1;
             New->index=D->index;
             New->core=myrank;
           }
         }
       }                //if(i,j,k<range)
     }                  //End of for(n<cnt)
     p=p->next;
   }                    //End of while(p)
}


void loadPolygonPlasma(Domain *D,LoadList *LL,int s,int istart,int iend,int iteration)
{
  	int nn,i,j,k,jstart,jend,kstart,kend,intNum,cnt,l,m,n;
  	int modeX,modeYZ;
  	double posX,posY,posZ,v1,v2,v3,weight,charge,centerX,centerY,centerZ,tmp,weightCoef,invNum,realDx,realDy,realDz,y0,z0;
  	double randTest,positionX[3],positionY[3],positionZ[3],gaussCoefX,polyCoefX,gaussCoefYZ,polyCoefYZ;
   double neX,neY,neZ,neFuncX,neFuncY,neFuncZ;
  	Particle ***particle;
  	particle=D->particle;
  	int myrank;
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  	double coef2D,coef3D;

  	ptclList *New,*p;   

  	jstart=D->jstart;   jend=D->jend;
  	kstart=D->kstart;   kend=D->kend;
  	centerX=LL->centerX;
  	centerY=LL->centerY;
  	centerZ=LL->centerZ;
  	gaussCoefX=LL->gaussCoefX;
  	polyCoefX=LL->polyCoefX;
  	gaussCoefYZ=LL->gaussCoefYZ;
  	polyCoefYZ=LL->polyCoefYZ;
  	modeX=LL->modeX;
  	modeYZ=LL->modeYZ;
  	charge=LL->charge;
  	if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;

  	srand(iteration+myrank);
  	intNum=(int)LL->numberInCell;

  	coef2D=coef3D=0.0;
  	if(D->dimension>1) coef2D=1; else;
  	if(D->dimension>2) coef3D=1; else;
  	if(LL->numberInCell==0) invNum=0.0;
  	else                    invNum=1.0/LL->numberInCell;
  	realDx=D->dx*D->lambda;
  	realDy=D->dy*D->lambda;
  	realDz=D->dz*D->lambda;

  	//position define   
	gsl_qrng *q3 = gsl_qrng_alloc (gsl_qrng_sobol,3);   
	for(i=istart; i<iend; i++)  
	{
 		neX=neFuncX=0.0;
		posX=(double)(i+D->minXSub-D->istart);
     	for(l=0; l<LL->xnodes-1; l++) {
      	if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1]) {
        		y0=evaluate_rpn(LL->rpn_x_y0[l],LL->rpn_size_x_y0[l],posX*realDx)*coef2D;
        		z0=evaluate_rpn(LL->rpn_x_z0[l],LL->rpn_size_x_z0[l],posX*realDx)*coef3D;
        		neX=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
        		neFuncX=evaluate_rpn(LL->rpn_x[l],LL->rpn_size_x[l],posX*realDx);
			} else ;
		}
		for(j=jstart; j<jend; j++)
		{
        	neY=neFuncY=0.0;
      	posY=(double)(j+D->minYSub-jstart+y0/realDy);
         for(m=0; m<LL->ynodes-1; m++) {
         	if(posY>=LL->ypoint[m] && posY<LL->ypoint[m+1]) {
         		neY=((LL->yn[m+1]-LL->yn[m])/(LL->ypoint[m+1]-LL->ypoint[m])*(posY-LL->ypoint[m])+LL->yn[m]);
         		neFuncY=evaluate_rpn(LL->rpn_y[m],LL->rpn_size_y[m],posY*realDy);
				} else ;
			}
      	for(k=kstart; k<kend; k++)
      	{
 				neZ=neFuncZ=0.0;
           	posZ=(double)(k+D->minZSub-kstart+z0/realDz);
           	for(n=0; n<LL->znodes-1; n++)	{
         		if(posZ>=LL->zpoint[n] && posZ<LL->zpoint[n+1])  {
		      		neZ=((LL->zn[n+1]-LL->zn[n])/(LL->zpoint[n+1]-LL->zpoint[n])*(posZ-LL->zpoint[n])+LL->zn[n]);
            		neFuncZ=evaluate_rpn(LL->rpn_z[n],LL->rpn_size_z[n],posZ*realDz);
            	} else ;
				}
            weight=neX*neY*neZ*neFuncX*neFuncY*neFuncZ*invNum*weightCoef;

		      cnt=0;
      		while(cnt<intNum/2)
            {               
              	positionX[0]=randomValue(1.0);
               positionY[0]=randomValue(1.0);
               positionZ[0]=randomValue(1.0);
               positionX[1]=1.0-positionX[0];
               positionY[1]=1.0-positionY[0];
               positionZ[1]=1.0-positionZ[0];
//                  random3D_sobol(&positionX,&positionY,&positionZ,q3);
               for(nn=0; nn<2; nn++) {
                	New = (ptclList *)malloc(sizeof(ptclList)); 
                 	New->next = particle[i][j][k].head[s]->pt;
                 	particle[i][j][k].head[s]->pt = New;
                 	New->x = positionX[nn];
                 	New->oldX=i+positionX[nn];
                 	New->y = positionY[nn]*coef2D;
                 	New->oldY=(j+positionY[nn])*coef2D;
                 	New->z = positionZ[nn]*coef3D;
                 	New->oldZ=(k +positionZ[nn])*coef3D;
                 	New->E1=New->E2=New->E3=0.0;
                	New->B1=New->B2=New->B3=0.0;
                	v1=maxwellianVelocity(LL->temperature,LL->mass)/velocityC;
                 	v2=maxwellianVelocity(LL->temperature,LL->mass)/velocityC;
                 	v3=maxwellianVelocity(LL->temperature,LL->mass)/velocityC;
                 	New->p1=-D->gamma*D->beta+v1;
                 	New->p2=v2;
                 	New->p3=v3;
                 	New->weight=weight;
                 	New->charge=charge;
                 	D->index+=1;
                 	New->index=D->index;            
                 	New->core=myrank;            
               } 
               cnt++;
            }		//end of while(cnt)
         }	   //End of for (k)
      }		//End of for(j)    
   }		//End of for(i)    
	gsl_qrng_free(q3);     
}

double maxwellianVelocity(double temperature,double mass)
{
   double vth,r,prob,v,random;
   int intRand,randRange=1e5;

   vth=sqrt(2.0*eCharge*temperature/(eMass*mass));
   
   r=1.0;
   prob=0.0;
   while (r>prob)  {
      intRand = rand() % randRange;
      r = ((double)intRand)/randRange;
      intRand = rand() % randRange;
      random = ((double)intRand)/randRange;
      v = 6.0*(random-0.5);
      prob=exp(-v*v);
   }
   return vth*v;
}

double gaussian_dist(double sig)
{
   double r,prob,v,random;
   int intRand,randRange=1e5;

   r=1.0;
   prob=0.0;
   while (r>prob)  {
      intRand = rand() % randRange;
      r = ((double)intRand)/randRange;
      intRand = rand() % randRange;
      random = ((double)intRand)/randRange;
      v = 4.0*(random-0.5);
      prob=exp(-v*v);
   }
   return sig*v;
}

void assignDefParticle(Domain *D)
{
  int n,cnt,flag;
  double tmp,tmpX,tmpY,tmpZ,tmpR,factY,factZ,randX,randY,randZ;
  LoadList *LL,*prevL;
  DefPtcl *p,*prevP;

  factY=factZ=0.0;
  if(D->dimension>1) factY=1.0; else;
  if(D->dimension>2) factZ=1.0; else;

  LL=D->loadList;
  while(LL->next)      {
    if(LL->type==Defined) {
      if(LL->pair==OFF) {
        p=LL->def;
        while(p) {
          if(p->flag==OFF) {
            if(p->xPos<=LL->createGuard) {
              cnt=p->numDefPtcls=LL->numDefPtcls;           
              p->define=(double *)malloc((cnt*3)*sizeof(double ));
              gsl_qrng *q3 = gsl_qrng_alloc (gsl_qrng_sobol,3);
              for(n=0; n<cnt; n++)  {
                flag=0;
                while(flag==0) {
//                  random2D_sobol(&randX,&randY,q3); randZ=0.0;
                  random3D_sobol(&randX,&randY,&randZ,q3);
                  tmpX=-LL->RDef+randX*2.0*LL->RDef;
                  tmpY=(-LL->RDef+randY*2.0*LL->RDef)*factY;
                  tmpZ=(-LL->RDef+randZ*2.0*LL->RDef)*factZ;
//                  tmp=(double)(randomValue(1.0));
//                  tmpX=-LL->RDef+tmp*2.0*LL->RDef;
//                  tmp=(double)(randomValue(1.0));
//                  tmpY=(-LL->RDef+tmp*2.0*LL->RDef)*factY;
//                  tmp=(double)(randomValue(1.0));
//                  tmpZ=(-LL->RDef+tmp*2.0*LL->RDef)*factZ;
                  tmpR=sqrt(tmpX*tmpX+tmpY*tmpY+tmpZ*tmpZ);
                  if(tmpR<=LL->RDef) {
                    flag=1;
                    tmpX=p->xPos+tmpX;
                    tmpY=p->yPos+tmpY;
                    tmpZ=p->zPos+tmpZ;
                  } else  flag=0;
                }
                p->define[n]=tmpX;
                p->define[n+cnt]=tmpY;
                p->define[n+2*cnt]=tmpZ;
              }	//End of for(n<cnt)
              p->flag=ON;  
            } else ;
          } else ;	//End of if(flag==Off)
          p=p->next;
        }
      } else { 	//if pair==ON
        p=LL->def; prevP=prevL->def;
        while(p) {
          if(p->flag==OFF) {
            if(p->xPos<=LL->createGuard) {
              cnt=p->numDefPtcls=prevL->numDefPtcls;           
              p->define=(double *)malloc((cnt*3)*sizeof(double ));
              for(n=0; n<cnt; n++)  {
                p->define[n]=prevP->define[n];
                p->define[n+cnt]=prevP->define[n+cnt];
                p->define[n+2*cnt]=prevP->define[n+2*cnt];
              }
              p->flag=ON;  
            } else ;
          } else ;	//End of if(flag==Off)
          p=p->next;
          prevP=prevP->next;	
        } 
      }
    } else ; 		//End of if(LL->type==Defined) 
    prevL=LL;
    LL=LL->next;
  }
}

void deleteDefParticle(Domain *D)
{
  int cnt;
  LoadList *LL;
  DefPtcl *p,*prev;

  LL=D->loadList;
  while(LL->next)      {
    if(LL->type==Defined) {
      p=LL->def;
      cnt=1;
      while(p) {
        if(cnt==1) prev=p; else ;
  
        if(p->xPos<=LL->deleteGuard) {
          if(cnt==1)  {
            LL->def=p->next;
            p->next=NULL;
            free(p->define);
            free(p);
            p=LL->def;
            cnt=1;
          } else {
            prev->next=p->next;
            p->next=NULL;
            free(p->define);
            free(p);
            p=prev->next;
          } 
        }
        else  {
          prev=p;
          p=p->next;
          cnt++;
        }
      }
    }	else ;	//End of if(LL->type==Defined)
    LL=LL->next;
  }
}

void saveDefParticle(Domain *D,int iteration)
{
  int n,cnt,num;
  double x,y,z;
  char name[100];
  FILE *out;

  LoadList *LL;
  DefPtcl *p;
  LL=D->loadList;
  num=0;
  while(LL->next)      {
    sprintf(name,"def%d_%d",iteration,num);
    out = fopen(name,"w");    

    p=LL->def;
    while(p) {
      cnt=p->numDefPtcls;
      for(n=0; n<cnt; n++)  {
        x=p->define[n]*D->lambda;
        y=p->define[n+cnt]*D->lambda;
        z=p->define[n+2*cnt]*D->lambda;
        fprintf(out,"%g %g %g\n",x,y,z);
      }	//End of for(n<cnt)
      p=p->next;
    }
    fclose(out);
    num++;
    LL=LL->next;
  }
}

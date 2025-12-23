#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"

void absorb_RL3(Domain *D);
void absorb_UD3(Domain *D);
void absorb_FB3(Domain *D);
void ExBx_solve1D_Split(Domain *D);
void PS_solve1D_Split(Domain *D);
void ExBx_solve2D_Split(Domain *D,int itertation);
void PS_solve2D_Split(Domain *D,int iteration);
void ExBx_solve3D_Split(Domain *D,int itertation);
void PS_solve3D_Split(Domain *D,int iteration);
void Bsolve2D_Pukhov(Domain *D,int iteration);
void Esolve2D_Pukhov(Domain *D,int iteration);
void Bsolve3D_Pukhov(Domain *D,int iteration);
void Esolve3D_Pukhov(Domain *D,int iteration);
void Bsolve3D_Yee(Domain *D,int iteration);
void Esolve3D_Yee(Domain *D,int iteration);
void Bsolve2D_Yee(Domain *D,int iteration);
void Esolve2D_Yee(Domain *D,int iteration);
void Bsolve1D_Yee_Pukhov(Domain *D);
void Esolve1D_Yee_Pukhov(Domain *D);
void Bfield_Treatment(Domain *D,int iteration);


void fieldSolve1(Domain D,double t,int iteration)
{
  	LaserList *L;
  	int myrank, nTasks,rank,rankM,rankN;

  	MPI_Status status;
  	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


   if(D.boostOn==OFF)    {
     L=D.laserList;
     while(L->next)  {
       if(L->loadMethod==Boundary) loadLaser(&D,L,t); else ;
		//         if(L->direction==1)     loadLaser2D(&D,L,t); 
		//         else if(L->direction==-1)     loadLaserOpp2D(&D,L,t); 
       L=L->next;
     }
   }

  	switch((D.fieldType-1)*3+D.dimension) {

  	// Split
  	case (Split-1)*3+1:
  		PS_solve1D_Split(&D);
  		break;

  	case (Split-1)*3+2:
  		PS_solve2D_Split(&D,iteration);
  		break;

  	case (Split-1)*3+3:
  		PS_solve3D_Split(&D,iteration);
  		break;

  	//Yee, Pukhov
  	case (Yee-1)*3+1:
  	case (Pukhov-1)*3+1:
  		Esolve1D_Yee_Pukhov(&D);
  		break;

  	case (Yee-1)*3+2:
  		Esolve2D_Yee(&D,iteration);
  		break;

  	case (Yee-1)*3+3:
  		Esolve3D_Yee(&D,iteration);
  		break;

  	case (Pukhov-1)*3+2:
  		Esolve2D_Pukhov(&D,iteration);
  		break;

  	case (Pukhov-1)*3+3:
  		Esolve3D_Pukhov(&D,iteration);
  		break;

  	default:
   	printf("what fieldType? and what dimension?\n");
	}

  	switch (D.fieldType) {
  	case Split:
  		D.shareF[0]=D.Pr;
  		D.shareF[1]=D.Pl;
  		D.shareF[2]=D.Sr;
  		D.shareF[3]=D.Sl;
  		MPI_TransferF_Xminus(&D,4);
  		MPI_TransferF_Xplus(&D,4);
  		MPI_TransferF_Yminus(&D,4);
  		MPI_TransferF_Yplus(&D,4);
  		MPI_TransferF_Zminus(&D,4);
  		MPI_TransferF_Zplus(&D,4);
  		if(D.Period==ON)  {
  			MPI_TransferF_Period_X(&D,4);
  			MPI_TransferF_Period_Y(&D,4);
  			MPI_TransferF_Period_Z(&D,4);
  		} else ;
  		D.Pr=D.shareF[0];
  		D.Pl=D.shareF[1];
  		D.Sr=D.shareF[2];
  		D.Sl=D.shareF[3];   
  		break;
  	case Yee:
  	case Pukhov:
  		D.shareF[0]=D.Ex;
  		D.shareF[1]=D.Ey;
  		D.shareF[2]=D.Ez;
  		MPI_TransferF_Xminus(&D,3);
  		MPI_TransferF_Xplus(&D,3);
  		MPI_TransferF_Yminus(&D,3);
  		MPI_TransferF_Yplus(&D,3);
  		MPI_TransferF_Zminus(&D,3);
  		MPI_TransferF_Zplus(&D,3);
  		if(D.Period==ON)  {
  			MPI_TransferF_Period_X(&D,3);
  			MPI_TransferF_Period_Y(&D,3);
  			MPI_TransferF_Period_Z(&D,3);
  		} else ;
  		D.Ex=D.shareF[0];
  		D.Ey=D.shareF[1];
  		D.Ez=D.shareF[2];
  		break;
  	default:
  		printf("what fieldType? and what dimension?\n");
  	}
}

void fieldSolve2(Domain D,double t,int iteration)
{
  LaserList *L;
  int myrank, nTasks,rank,rankM,rankN;

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);



  switch((D.fieldType-1)*3+D.dimension) {

  // Split
  	case (Split-1)*3+1:
  		ExBx_solve1D_Split(&D);
  		break;
	
  	case (Split-1)*3+2:
  		ExBx_solve2D_Split(&D,iteration);
  		break;
	
  	case (Split-1)*3+3:
  		ExBx_solve3D_Split(&D,iteration);  
  		break;

  	//Yee, Pukhov
  	case (Yee-1)*3+1:
  	case (Pukhov-1)*3+1:
  		Bsolve1D_Yee_Pukhov(&D);  
  		break;

  	case (Yee-1)*3+2:
  		Bsolve2D_Yee(&D,iteration);  
  		break;

  	case (Yee-1)*3+3:
  		Bsolve3D_Yee(&D,iteration);
  		break;

  	case (Pukhov-1)*3+2:
		Bsolve2D_Pukhov(&D,iteration);
		break;

  	case (Pukhov-1)*3+3:
		Bsolve3D_Pukhov(&D,iteration);
		break;

	default:
   	printf("what fieldType? and what dimension?\n");
	}

  	switch (D.fieldType) {
	  
  	case Split:
  		D.shareF[0]=D.Ex;
  		D.shareF[1]=D.Bx;
  		D.shareF[2]=D.ExNow;
  		D.shareF[3]=D.BxNow;
  		MPI_TransferF_Xminus(&D,4);
  		MPI_TransferF_Xplus(&D,4);
  		MPI_TransferF_Yminus(&D,4);
  		MPI_TransferF_Yplus(&D,4);
  		MPI_TransferF_Zminus(&D,4);
  		MPI_TransferF_Zplus(&D,4);
  		if(D.Period==ON)  {
  			MPI_TransferF_Period_X(&D,4);
  			MPI_TransferF_Period_Y(&D,4);
  			MPI_TransferF_Period_Z(&D,4);
  		} else ;
  		D.Ex=D.shareF[0];
  		D.Bx=D.shareF[1];
  		D.ExNow=D.shareF[2];
  		D.BxNow=D.shareF[3];   
  		break;
  	case Yee:
  	case Pukhov:
  		D.shareF[0]=D.Bx;
  		D.shareF[1]=D.By;
  		D.shareF[2]=D.Bz;
  		D.shareF[3]=D.BxNow;
  		D.shareF[4]=D.ByNow;
  		D.shareF[5]=D.BzNow;
  		MPI_TransferF_Xminus(&D,6);
  		MPI_TransferF_Xplus(&D,6);
  		MPI_TransferF_Yminus(&D,6);
  		MPI_TransferF_Yplus(&D,6);
  		MPI_TransferF_Zminus(&D,6);
  		MPI_TransferF_Zplus(&D,6);
  		if(D.Period==ON)  {
  			MPI_TransferF_Period_X(&D,6);
  			MPI_TransferF_Period_Y(&D,6);
  			MPI_TransferF_Period_Z(&D,6);
  		} else ;
  		D.Bx=D.shareF[0];
  		D.By=D.shareF[1];
  		D.Bz=D.shareF[2];
  		D.BxNow=D.shareF[3];
  		D.ByNow=D.shareF[4];  
  		D.BzNow=D.shareF[5];  
  		break;
  	default:
  		printf("what fieldType? and what dimension?\n");
  	}

}


void Esolve3D_Yee(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,rankY,rankZ;
    double dtOverdx,dtOverdy,dtOverdz,dbPiDt;
    double oldEx,oldEy,oldEz;
    double rtr,rtd,ltr,ltd,upr,upd,dnr,dnd;
    double tmp,tmpr,tmpd;

    int nTasks,myrank;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    istart=D->istart; iend=D->iend;
    jstart=D->jstart; jend=D->jend;
    kstart=D->kstart; kend=D->kend;

    for(k=kstart; k<kend; k++) {
      D->frr[k]=1.0;
      D->frd[k]=1.0;
      D->bkr[k]=1.0;
      D->bkd[k]=1.0;
    }
    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0;
      D->upd[j]=1.0;
      D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_FB3(D);
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    dtOverdy=D->dt/D->dy; 
	 dtOverdx=D->dt/D->dx; 
	 dtOverdz=D->dt/D->dz;
    dbPiDt=2.0*M_PI*D->dt;	 

    //Solving E field
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        upr=D->upr[j]; dnr=D->dnr[j];
        upd=D->upd[j]; dnd=D->dnd[j];
        for(k=kstart; k<kend; k++)
        {
          tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
          tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];	 

          oldEx=D->Ex[i][j][k];
          oldEy=D->Ey[i][j][k];
          oldEz=D->Ez[i][j][k];

          tmp=dtOverdy*(D->Bz[i][j][k]-D->Bz[i][j-1][k])-dtOverdz*(D->By[i][j][k]-D->By[i][j][k-1])-dbPiDt*D->Jx[i][j][k];
          D->Ex[i][j][k]=tmpd*(oldEx+tmp*tmpr);

          tmp=-dtOverdx*(D->Bz[i][j][k]-D->Bz[i-1][j][k])+dtOverdz*(D->Bx[i][j][k]-D->Bx[i][j][k-1])-dbPiDt*D->Jy[i][j][k];
          D->Ey[i][j][k]=tmpd*(oldEy+tmp*tmpr);

          tmp=dtOverdx*(D->By[i][j][k]-D->By[i-1][j][k])-dtOverdy*(D->Bx[i][j][k]-D->Bx[i][j-1][k])-dbPiDt*D->Jz[i][j][k];
          D->Ez[i][j][k]=tmpd*(oldEz+tmp*tmpr);
        }
      }
    }

   //  rankY=(myrank%(D->M*D->N))%D->M;
   //  rankZ=(myrank%(D->M*D->N))/D->M;
	//  if(rankY==0) {
	// 	j=jstart;
   //    for(i=istart; i<iend; i++)
   //      for(k=kstart; k<kend; k++)
	//        D->Ex[i][j][k]=0.0;
	//  } else ;
	//  if(rankZ==0) {
	// 	k=kstart;
   //    for(i=istart; i<iend; i++)
   //      for(j=jstart; j<jend; j++)
	//        D->Ex[i][j][k]=0.0;
	//  } else ;
}

void Bsolve3D_Yee(Domain *D,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend;
	double ax,ay,az,bx,by,bz,dtOverdx,dtOverdy,dtOverdz;
	double oldBx,oldBy,oldBz;
	double rtr,rtd,ltr,ltd,upr,upd,dnr,dnd;
   double tmp,tmpr,tmpd;

   int nTasks,myrank;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;    iend=D->iend;
   jstart=D->jstart;    jend=D->jend;
   kstart=D->kstart;    kend=D->kend;

    for(k=kstart; k<kend; k++) {
      D->frr[k]=1.0;
      D->frd[k]=1.0;
      D->bkr[k]=1.0;
      D->bkd[k]=1.0;
    }
    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0;
      D->upd[j]=1.0;
      D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_FB3(D);
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    dtOverdy=D->dt/D->dy; 
	 dtOverdx=D->dt/D->dx; 
	 dtOverdz=D->dt/D->dz;

    //Solving E field
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        upr=D->upr[j]; dnr=D->dnr[j];
        upd=D->upd[j]; dnd=D->dnd[j];
        for(k=kstart; k<kend; k++)
        {
          tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
          tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];	 

          oldBx=D->Bx[i][j][k];
          oldBy=D->By[i][j][k];
          oldBz=D->Bz[i][j][k];

          tmp=-dtOverdy*(D->Ez[i][j+1][k]-D->Ez[i][j][k])+dtOverdz*(D->Ey[i][j][k+1]-D->Ey[i][j][k]);
          D->Bx[i][j][k]=tmpd*(oldBx+tmp*tmpr);

          tmp=dtOverdx*(D->Ez[i+1][j][k]-D->Ez[i][j][k])-dtOverdz*(D->Ex[i][j][k+1]-D->Ex[i][j][k]);
          D->By[i][j][k]=tmpd*(oldBy+tmp*tmpr);

          tmp=dtOverdy*(D->Ex[i][j+1][k]-D->Ex[i][j][k])-dtOverdx*(D->Ey[i+1][j][k]-D->Ey[i][j][k]);
          D->Bz[i][j][k]=tmpd*(oldBz+tmp*tmpr);

          D->BxNow[i][j][k]=0.5*(D->Bx[i][j][k]+oldBx);
          D->ByNow[i][j][k]=0.5*(D->By[i][j][k]+oldBy);
          D->BzNow[i][j][k]=0.5*(D->Bz[i][j][k]+oldBz);
        }
      }
    }
}

void ExBx_solve1D_Split(Domain *D)
{
  int i,j,k,istart,iend,nxSub;
  double dx,dt;

  dx=D->dx;          dt=D->dt;
  istart=D->istart;  iend=D->iend;

  j=k=0;
  for(i=iend-1; i>=istart; i--)
    D->Ex[i][j][k]+=-2.0*M_PI*dt*D->Jx[i][j][k];  
}

void PS_solve1D_Split(Domain *D)
{
  int i,j,k,istart,iend,nxSub;
  double dx,dt;

  dx=D->dx;          dt=D->dt;
  istart=D->istart;  iend=D->iend;
  nxSub=D->nxSub;

  j=k=0;

  for(i=istart; i<iend; i++)
  {
    D->Pl[i][j][k]=D->Pl[i+1][j][k]-M_PI*dt*D->Jy[i][j][k];
    D->Sl[i][j][k]=D->Sl[i+1][j][k]-M_PI*dt*D->Jz[i][j][k];
  }  

  for(i=iend-1; i>=istart; i--)
  {
    D->Pr[i][j][k]=D->Pr[i-1][j][k]-M_PI*dt*D->Jy[i][j][k];
    D->Sr[i][j][k]=D->Sr[i-1][j][k]-M_PI*dt*D->Jz[i][j][k];
  }  

}

void Bsolve1D_Yee_Pukhov(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt,oldBx,oldBy,oldBz;

    dx=D->dx;    dy=D->dy;    dt=D->dt; 
    nxSub=D->nxSub;    nySub=D->nySub;
    
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
        {
          oldBx=D->Bx[i][j][k];
          oldBy=D->By[i][j][k];
          oldBz=D->Bz[i][j][k];
          D->Bx[i][j][k]=0.0;
          D->By[i][j][k]+=dt/dx*(D->Ez[i+1][j][k]-D->Ez[i][j][k]);
          D->Bz[i][j][k]+=-dt/dx*(D->Ey[i+1][j][k]-D->Ey[i][j][k]);
          D->BxNow[i][j][k]=0.5*(D->Bx[i][j][k]+oldBx);
          D->ByNow[i][j][k]=0.5*(D->By[i][j][k]+oldBy);
          D->BzNow[i][j][k]=0.5*(D->Bz[i][j][k]+oldBz);
        }
}

void Esolve1D_Yee_Pukhov(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt;

    dx=D->dx;    dy=D->dy;    dz=D->dz;    dt=D->dt;
    nxSub=D->nxSub;    nySub=D->nySub;    nzSub=D->nzSub;
    
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
        {
          D->Ex[i][j][k]+=-2*M_PI*dt*D->Jx[i][j][k];
          D->Ey[i][j][k]+=-dt/dx*(D->Bz[i][j][k]-D->Bz[i-1][j][k])-2*M_PI*dt*D->Jy[i][j][k];
          D->Ez[i][j][k]+=dt/dx*(D->By[i][j][k]-D->By[i-1][j][k])-2*M_PI*dt*D->Jz[i][j][k];
        }
}

void Bsolve3D_Pukhov(Domain *D,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend;
	double ax,ay,az,bx,by,bz,dtOverdx,dtOverdy,dtOverdz;
	double oldBx,oldBy,oldBz;
	double rtr,rtd,ltr,ltd,upr,upd,dnr,dnd;
	double tmp,tmpr,tmpd;

   int nTasks,myrank;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;    iend=D->iend;
   jstart=D->jstart;    jend=D->jend;
   kstart=D->kstart;    kend=D->kend;

   for(k=kstart; k<kend; k++) {
     D->frr[k]=1.0;
     D->frd[k]=1.0;
     D->bkr[k]=1.0;
     D->bkd[k]=1.0;
   }
   for(j=jstart; j<jend; j++) {
     D->upr[j]=1.0;
     D->upd[j]=1.0;
     D->dnr[j]=1.0;
     D->dnd[j]=1.0;
   }
   for(i=istart; i<iend; i++) {
     D->rtr[i]=1.0;
     D->rtd[i]=1.0;
     D->ltr[i]=1.0;
     D->ltd[i]=1.0;
   }

   ay=0.125*D->dx/D->dy;
   az=0.125*D->dx/D->dz;
   ax=ay+az;
   bx=1.0-2.0*ax;
   by=1.0-2.0*ay;
   bz=1.0-2.0*az;
   dtOverdx=D->dt/D->dx;
   dtOverdy=D->dt/D->dy;
   dtOverdz=D->dt/D->dz;

   for(i=istart; i<iend; i++)
   {
     rtr=D->rtr[i]; ltr=D->ltr[i];
     rtd=D->rtd[i]; ltd=D->ltd[i];
     for(j=jstart; j<jend; j++)
     {
       upr=D->upr[j]; dnr=D->dnr[j];
       upd=D->upd[j]; dnd=D->dnd[j];
       for(k=kstart; k<kend; k++)
       {
         tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
         tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];

         oldBx=D->Bx[i][j][k];
         oldBy=D->By[i][j][k];
         oldBz=D->Bz[i][j][k];

         tmp=-dtOverdy*(bz*(D->Ez[i][j+1][k]-D->Ez[i][j][k])+az*(D->Ez[i][j+1][k+1]+D->Ez[i][j+1][k-1]-D->Ez[i][j][k+1]-D->Ez[i][j][k-1]))+dtOverdz*(by*(D->Ey[i][j][k+1]-D->Ey[i][j][k])+ay*(D->Ey[i][j+1][k+1]+D->Ey[i][j-1][k+1]-D->Ey[i][j+1][k]-D->Ey[i][j-1][k]));
         D->Bx[i][j][k]=tmpd*(oldBx+tmp*tmpr);

         tmp=dtOverdx*(bz*(D->Ez[i+1][j][k]-D->Ez[i][j][k])+az*(D->Ez[i+1][j][k+1]+D->Ez[i+1][j][k-1]-D->Ez[i][j][k+1]-D->Ez[i][j][k-1]))-dtOverdz*(bx*(D->Ex[i][j][k+1]-D->Ex[i][j][k])+ax*(D->Ex[i+1][j][k+1]+D->Ex[i-1][j][k+1]-D->Ex[i+1][j][k]-D->Ex[i-1][j][k]));
         D->By[i][j][k]=tmpd*(oldBy+tmp*tmpr);

         tmp=dtOverdy*(bx*(D->Ex[i][j+1][k]-D->Ex[i][j][k])+ax*(D->Ex[i+1][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i+1][j][k]-D->Ex[i-1][j][k]))-dtOverdx*(by*(D->Ey[i+1][j][k]-D->Ey[i][j][k])+ay*(D->Ey[i+1][j+1][k]+D->Ey[i+1][j-1][k]-D->Ey[i][j+1][k]-D->Ey[i][j-1][k]));
         D->Bz[i][j][k]=tmpd*(oldBz+tmp*tmpr);

         D->BxNow[i][j][k]=0.5*(D->Bx[i][j][k]+oldBx);
         D->ByNow[i][j][k]=0.5*(D->By[i][j][k]+oldBy);
         D->BzNow[i][j][k]=0.5*(D->Bz[i][j][k]+oldBz);
       }
     }
   }
}

void Esolve3D_Pukhov(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;
    double ax,ay,az,bx,by,bz,dtOverdx,dtOverdy,dtOverdz,dbPiDt;
    double oldEx,oldEy,oldEz;
    double rtr,rtd,ltr,ltd,upr,upd,dnr,dnd;
    double tmp,tmpr,tmpd;

    int nTasks,myrank,rankY,rankZ;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

   for(k=kstart; k<kend; k++) {
     D->frr[k]=1.0;
     D->frd[k]=1.0;
     D->bkr[k]=1.0;
     D->bkd[k]=1.0;
   }
   for(j=jstart; j<jend; j++) {
     D->upr[j]=1.0;
     D->upd[j]=1.0;
     D->dnr[j]=1.0;
     D->dnd[j]=1.0;
   }
   for(i=istart; i<iend; i++) {
     D->rtr[i]=1.0;
     D->rtd[i]=1.0;
     D->ltr[i]=1.0;
     D->ltd[i]=1.0;
   }

    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_FB3(D);
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    ay=0.125*D->dx/D->dy;
    az=0.125*D->dx/D->dz;
    ax=ay+az;
    bx=1.0-2.0*ax;
    by=1.0-2.0*ay;
    bz=1.0-2.0*az;
    dtOverdx=D->dt/D->dx;
    dtOverdy=D->dt/D->dy;
    dtOverdz=D->dt/D->dz;
    dbPiDt=2.0*M_PI*D->dt;

    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_FB3(D);
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    ay=0.125*D->dx/D->dy;
    az=0.125*D->dx/D->dz;
    ax=ay+az;
    bx=1.0-2.0*ax;
    by=1.0-2.0*ay;
    bz=1.0-2.0*az;
    dtOverdx=D->dt/D->dx;
    dtOverdy=D->dt/D->dy;
    dtOverdz=D->dt/D->dz;
    dbPiDt=2.0*M_PI*D->dt;

    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        upr=D->upr[j]; dnr=D->dnr[j];
        upd=D->upd[j]; dnd=D->dnd[j];
        for(k=kstart; k<kend; k++)
        {
          tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
          tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];

          oldEx=D->Ex[i][j][k];
          oldEy=D->Ey[i][j][k];
          oldEz=D->Ez[i][j][k];

          tmp=dtOverdy*(bz*(D->Bz[i][j][k]-D->Bz[i][j-1][k])+az*(D->Bz[i][j][k+1]+D->Bz[i][j][k-1]-D->Bz[i][j-1][k+1]-D->Bz[i][j-1][k-1]))-dtOverdz*(by*(D->By[i][j][k]-D->By[i][j][k-1])+ay*(D->By[i][j+1][k]+D->By[i][j-1][k]-D->By[i][j+1][k-1]-D->By[i][j-1][k-1]))-dbPiDt*D->Jx[i][j][k];
          D->Ex[i][j][k]=tmpd*(oldEx+tmp*tmpr);

          tmp=-dtOverdx*(bz*(D->Bz[i][j][k]-D->Bz[i-1][j][k])+az*(D->Bz[i][j][k+1]+D->Bz[i][j][k-1]-D->Bz[i-1][j][k+1]-D->Bz[i-1][j][k-1]))+dtOverdz*(bx*(D->Bx[i][j][k]-D->Bx[i][j][k-1])+ax*(D->Bx[i+1][j][k]+D->Bx[i-1][j][k]-D->Bx[i+1][j][k-1]-D->Bx[i-1][j][k-1]))-dbPiDt*D->Jy[i][j][k];
          D->Ey[i][j][k]=tmpd*(oldEy+tmp*tmpr);

          tmp=dtOverdx*(by*(D->By[i][j][k]-D->By[i-1][j][k])+ay*(D->By[i][j+1][k]+D->By[i][j-1][k]-D->By[i-1][j+1][k]-D->By[i-1][j-1][k]))-dtOverdy*(bx*(D->Bx[i][j][k]-D->Bx[i][j-1][k])+ax*(D->Bx[i+1][j][k]+D->Bx[i-1][j][k]-D->Bx[i+1][j-1][k]-D->Bx[i-1][j-1][k]))-dbPiDt*D->Jz[i][j][k];
          D->Ez[i][j][k]=tmpd*(oldEz+tmp*tmpr);
        }
      }
    }

	//  rankY=(myrank%(D->M*D->N))%D->M;
   //  rankZ=(myrank%(D->M*D->N))/D->M;
	//  if(rankY==0) {
	// 	j=jstart;
   //    for(i=istart; i<iend; i++)
   //      for(k=kstart; k<kend; k++)
	//        D->Ex[i][j][k]=0.0;
	//  } else ;
	//  if(rankZ==0) {
	// 	k=kstart;
   //    for(i=istart; i<iend; i++)
   //      for(j=jstart; j<jend; j++)
	//        D->Ex[i][j][k]=0.0;
	//  } else ;
}

void Bsolve2D_Pukhov(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;  
    double oldBx,oldBy,oldBz,ax,ay,bx,by,x,y;
    double dtOverdx,dtOverdy,rtr,ltr,rtd,ltd;
    double tmp,tmpr,tmpd;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;

    k=0;
    ay=ax=0.125*D->dx/D->dy;
    bx=1.0-2.0*ax;
    by=1.0-2.0*ay;
    dtOverdy=D->dt/D->dy;
    dtOverdx=D->dt/D->dx;

    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0;                                                                                                       D->upd[j]=1.0;                                                                                                       D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    //Solving B field
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        oldBx=D->Bx[i][j][k];
        oldBy=D->By[i][j][k];
        oldBz=D->Bz[i][j][k];

        tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
		  tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

        tmp=-dtOverdy*(D->Ez[i][j+1][k]-D->Ez[i][j][k]);
        D->Bx[i][j][k]=tmpd*(oldBx+tmp*tmpr);

        tmp=dtOverdx*(D->Ez[i+1][j][k]-D->Ez[i][j][k]);
        D->By[i][j][k]=tmpd*(oldBy+tmp*tmpr);

        tmp=dtOverdy*(bx*(D->Ex[i][j+1][k]-D->Ex[i][j][k])+ax*(D->Ex[i+1][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i+1][j][k]-D->Ex[i-1][j][k]))-dtOverdx*(by*(D->Ey[i+1][j][k]-D->Ey[i][j][k])+ay*(D->Ey[i+1][j+1][k]+D->Ey[i+1][j-1][k]-D->Ey[i][j+1][k]-D->Ey[i][j-1][k]));
        D->Bz[i][j][k]=tmpd*(oldBz+tmp*tmpr);

        D->BxNow[i][j][k]=0.5*(D->Bx[i][j][k]+oldBx);
        D->ByNow[i][j][k]=0.5*(D->By[i][j][k]+oldBy);
        D->BzNow[i][j][k]=0.5*(D->Bz[i][j][k]+oldBz);
      }
    }
}



void Esolve2D_Pukhov(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;  
    double dx,dy,dz,dt,ax,ay,bx,by,x1,x2,x,y;
    double oldEx,oldEy,oldEz;
    double dtOverdx,dtOverdy,dbPiDt,rtr,ltr,rtd,ltd;
    double tmp,tmpr,tmpd;

	 int nTasks,myrank;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0;                                                                                                       D->upd[j]=1.0;                                                                                                       D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    ay=ax=0.125*D->dx/D->dy;
    bx=1.0-2.0*ax;
    by=1.0-2.0*ay;
    dtOverdy=D->dt/D->dy;
    dtOverdx=D->dt/D->dx;
	 dbPiDt=2.0*M_PI*D->dt;

    k=0;
    //Solving E field
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
		  tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

        oldEx=D->Ex[i][j][k];
        oldEy=D->Ey[i][j][k];
        oldEz=D->Ez[i][j][k];

        tmp=dtOverdy*(D->Bz[i][j][k]-D->Bz[i][j-1][k])-dbPiDt*D->Jx[i][j][k];
        D->Ex[i][j][k]=tmpd*(oldEx+tmp*tmpr);

        tmp=-dtOverdx*(D->Bz[i][j][k]-D->Bz[i-1][j][k])-dbPiDt*D->Jy[i][j][k];
        D->Ey[i][j][k]=tmpd*(oldEy+tmp*tmpr);

        tmp=dtOverdx*(by*(D->By[i][j][k]-D->By[i-1][j][k])+ay*(D->By[i][j+1][k]+D->By[i][j-1][k]-D->By[i-1][j+1][k]-D->By[i-1][j-1][k]))-dtOverdy*(bx*(D->Bx[i][j][k]-D->Bx[i][j-1][k])+ax*(D->Bx[i+1][j][k]+D->Bx[i-1][j][k]-D->Bx[i+1][j-1][k]-D->Bx[i-1][j-1][k]))-dbPiDt*D->Jz[i][j][k];
        D->Ez[i][j][k]=tmpd*(oldEz+tmp*tmpr);
      }

//		j=jstart-1;
//      D->Ey[i][j][k]-=dtOverdx*(D->Bz[i][j][k]-D->Bz[i-1][j][k]);
    }

 	//  if(myrank%D->M==0 && D->Period==OFF) {
   //   for(i=istart; i<iend; i++) D->Ex[i][jstart][k]=0.0;
   // } else ;
}

void Bsolve2D_Yee(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;  
    double oldBx,oldBy,oldBz,ax,ay,bx,by;
    double dtOverdx,dtOverdy,rtr,ltr,rtd,ltd;
    double tmp,tmpr,tmpd;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;

    k=0;
    ay=ax=0.125*D->dx/D->dy;
    bx=1.0-2.0*ax;
    by=1.0-2.0*ay;
    dtOverdy=D->dt/D->dy;
    dtOverdx=D->dt/D->dx;

    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0;
      D->upd[j]=1.0;
      D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    //Solving B field
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
        tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

        oldBx=D->Bx[i][j][k];
        oldBy=D->By[i][j][k];
        oldBz=D->Bz[i][j][k];

        tmp=-dtOverdy*(D->Ez[i][j+1][k]-D->Ez[i][j][k]);
        D->Bx[i][j][k]=tmpd*(oldBx+tmp*tmpr);

        tmp=dtOverdx*(D->Ez[i+1][j][k]-D->Ez[i][j][k]);
        D->By[i][j][k]=tmpd*(oldBy+tmp*tmpr);

        tmp=dtOverdy*(D->Ex[i][j+1][k]-D->Ex[i][j][k])-dtOverdx*(D->Ey[i+1][j][k]-D->Ey[i][j][k]);
        D->Bz[i][j][k]=tmpd*(oldBz+tmp*tmpr);

        D->BxNow[i][j][k]=0.5*(D->Bx[i][j][k]+oldBx);
        D->ByNow[i][j][k]=0.5*(D->By[i][j][k]+oldBy);
        D->BzNow[i][j][k]=0.5*(D->Bz[i][j][k]+oldBz);
      }
    }

}

void Esolve2D_Yee(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;  
    double dx,dy,dz,dt,ax,ay,bx,by;
    double oldEx,oldEy,oldEz;
    double dtOverdx,dtOverdy,dbPiDt,rtr,ltr,rtd,ltd;
    double tmp,tmpr,tmpd;
    int nTasks,myrank;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0;
      D->upd[j]=1.0;
      D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    ay=ax=0.125*D->dx/D->dy;
    bx=1.0-2.0*ax;
    by=1.0-2.0*ay;
    dtOverdy=D->dt/D->dy;
    dtOverdx=D->dt/D->dx;
    dbPiDt=2.0*M_PI*D->dt;

    k=0;
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
        tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

        oldEx=D->Ex[i][j][k];
        oldEy=D->Ey[i][j][k];
        oldEz=D->Ez[i][j][k];

        tmp=dtOverdy*(D->Bz[i][j][k]-D->Bz[i][j-1][k])-dbPiDt*D->Jx[i][j][k];
        D->Ex[i][j][k]=tmpd*(oldEx+tmp*tmpr);

        tmp=-dtOverdx*(D->Bz[i][j][k]-D->Bz[i-1][j][k])-dbPiDt*D->Jy[i][j][k];
        D->Ey[i][j][k]=tmpd*(oldEy+tmp*tmpr);

        tmp=dtOverdx*(D->By[i][j][k]-D->By[i-1][j][k])-dtOverdy*(D->Bx[i][j][k]-D->Bx[i][j-1][k])-dbPiDt*D->Jz[i][j][k];
        D->Ez[i][j][k]=tmpd*(oldEz+tmp*tmpr);
     }

   }

   // if(myrank%D->M==0 && D->Period==OFF) {
   //   for(i=istart; i<iend; i++) D->Ex[i][jstart][k]=0.0;
   // } else ;
}

void ExBx_solve2D_Split(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub;  
    double dtOverdx,dtOverdy,piDt,qtDtByDy,hfPiDt,dF;
    double tmp,tmpd,tmpr,ltr,rtr,ltd,rtd,ExOld,BxOld;

    int nTasks,myrank;
    MPI_Status status;          
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    
    
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend; k=0;
    dF=D->dF;

    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0; 
      D->upd[j]=1.0;
      D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    dtOverdy=D->dt/D->dy;
    dtOverdx=D->dt/D->dx;
    piDt=M_PI*D->dt;
    qtDtByDy=0.25*D->dt/D->dy;
    hfPiDt=0.5*M_PI*D->dt;

    // PrC,PlC,E1C,SrC,SlC,B1C
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        ExOld=D->Ex[i][j][k];
        BxOld=D->Bx[i][j][k];

        tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
        tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

        tmp=dtOverdy*(D->Pr[i][j][k]-D->Pr[i][j-1][k]-D->Pl[i][j][k]+D->Pl[i][j-1][k])-2*piDt*D->Jx[i][j][k]+dF*dtOverdx*(D->F[i+1][j][k]-D->F[i][j][k]);
        D->Ex[i][j][k]=tmpd*(D->Ex[i][j][k]+tmp*tmpr);

        tmp=-dtOverdy*(D->Sr[i][j+1][k]-D->Sr[i][j][k]+D->Sl[i][j+1][k]-D->Sl[i][j][k]);
        D->Bx[i][j][k]=tmpd*(D->Bx[i][j][k]+tmp*tmpr);

        D->ExNow[i][j][k]=0.5*(ExOld+D->Ex[i][j][k]);
        D->BxNow[i][j][k]=0.5*(BxOld+D->Bx[i][j][k]);
      }	
    }

     if(myrank%D->M==0 && D->Period==OFF) {
       for(i=istart; i<iend; i++) {
         D->Ex[i][jstart][k]=0.0;
         D->ExNow[i][jstart][k]=0.0;
       }
     } else ;
}

void PS_solve2D_Split(Domain *D,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend;  
   double dtOverdx,dtOverdy,piDt,qtDtByDy,dbPiDt;
   double tmp,tmpr,tmpd,ltr,ltd,rtr,rtd,dF,oldPr,oldPl,oldSr,oldSl;
   int nTasks,myrank;

   MPI_Status status;          
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    

   istart=D->istart;    iend=D->iend;
   jstart=D->jstart;    jend=D->jend;
   kstart=D->kstart;    kend=D->kend; k=0;
   dF=D->dF;

   for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0; 
      D->upd[j]=1.0;
      D->dnr[j]=1.0;
      D->dnd[j]=1.0;
   }
   for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
   }
   if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_UD3(D);
      absorb_RL3(D);
   } else ;

   dtOverdy=D->dt/D->dy;
   piDt=M_PI*D->dt;
   qtDtByDy=0.25*D->dt/D->dy;
   dbPiDt=2.0*M_PI*D->dt;

   for(i=istart; i<iend; i++)
   {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
         oldPl=D->Pl[i][j][k];
         oldSl=D->Sl[i][j][k];
 	 tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
	 tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

     	 tmp=-qtDtByDy*(D->Ex[i+1][j+1][k]+D->Ex[i][j+1][k]-D->Ex[i+1][j][k]-D->Ex[i][j][k])-piDt*D->Jy[i+1][j][k]+dF*dtOverdy*(D->F[i+1][j+1][k]-D->F[i+1][j][k]);
     	 D->Pl[i][j][k]=tmpd*(D->Pl[i+1][j][k]+tmp*tmpr);

     	 tmp=-qtDtByDy*(D->Bx[i+1][j][k]+D->Bx[i][j][k]-D->Bx[i+1][j-1][k]-D->Bx[i][j-1][k])-piDt*D->Jz[i+1][j][k];
     	 D->Sl[i][j][k]=tmpd*(D->Sl[i+1][j][k]+tmp*tmpr);

         D->PlC[i][j][k]=0.5*(oldPl+D->Pl[i][j][k]);
         D->SlC[i][j][k]=0.5*(oldSl+D->Sl[i][j][k]);				 
      }
   }

   for(i=iend-1; i>=istart; i--)
   {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
         oldPr=D->Pr[i][j][k];
         oldSr=D->Sr[i][j][k];			  
     	 tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
	 tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

     	 tmp=qtDtByDy*(D->Ex[i][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-piDt*D->Jy[i][j][k]+dF*dtOverdy*(D->F[i][j+1][k]-D->F[i][j][k]);
     	 D->Pr[i][j][k]=tmpd*(D->Pr[i-1][j][k]+tmp*tmpr);

     	 tmp=-qtDtByDy*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j-1][k]-D->Bx[i-1][j-1][k])-piDt*D->Jz[i][j][k];
     	 D->Sr[i][j][k]=tmpd*(D->Sr[i-1][j][k]+tmp*tmpr);

         D->PrC[i][j][k]=0.5*(oldPr+D->Pr[i][j][k]);
         D->SrC[i][j][k]=0.5*(oldSr+D->Sr[i][j][k]);					 
      }	
   }

}


void ExBx_solve3D_Split(Domain *D,int iteration)
{
  	int i,j,k,istart,iend,jstart,jend,kstart,kend;
  	double dtOverdx,dtOverdy,dtOverdz,qtDtByDy,qtDtByDz,hfPiDt,piDt;
  	double rtr,rtd,ltr,ltd,upr,upd,dnr,dnd,dF;
  	double tmp,tmpr,tmpd,ExOld,BxOld;

  	int nTasks,myrank,rankY,rankZ;
  	MPI_Status status;
  	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  	istart=D->istart;    iend=D->iend;
  	jstart=D->jstart;    jend=D->jend;
  	kstart=D->kstart;    kend=D->kend;
  	dF=D->dF;

  	for(k=kstart; k<kend; k++) {
  	  	D->frr[k]=1.0;
  	  	D->frd[k]=1.0;
  	  	D->bkr[k]=1.0;
  	  	D->bkd[k]=1.0;
  	}
  	for(j=jstart; j<jend; j++) {
  	  	D->upr[j]=1.0;
  	  	D->upd[j]=1.0;
  	  	D->dnr[j]=1.0;
  	  	D->dnd[j]=1.0;
  	}
  	for(i=istart; i<iend; i++) {
  	  	D->rtr[i]=1.0;
  	  	D->rtd[i]=1.0;
  	  	D->ltr[i]=1.0;
  	  	D->ltd[i]=1.0;
  	}
  	if(D->pmlOn==ON && iteration>D->pmlStart) {
  	  	absorb_FB3(D);
  	  	absorb_UD3(D);
  	  	absorb_RL3(D);
  	} else ;

  	dtOverdx=D->dt/D->dx;
  	dtOverdy=D->dt/D->dy;
  	dtOverdz=D->dt/D->dz;
  	piDt=M_PI*D->dt;
  	qtDtByDy=0.25*D->dt/D->dy;
  	qtDtByDz=0.25*D->dt/D->dz;
  	hfPiDt=0.5*M_PI*D->dt;

  	for(i=istart; i<iend; i++)
  	{
  	  	rtr=D->rtr[i]; ltr=D->ltr[i];
  	  	rtd=D->rtd[i]; ltd=D->ltd[i];
  	  	for(j=jstart; j<jend; j++)
  	  	{
  	    	upr=D->upr[j]; dnr=D->dnr[j];
  	    	upd=D->upd[j]; dnd=D->dnd[j];
  	    	for(k=kstart; k<kend; k++)
  	    	{
  	      	tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
  	      	tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];

  	      	ExOld=D->Ex[i][j][k];
  	      	BxOld=D->Bx[i][j][k];

  	      	tmp=dtOverdy*(D->Pr[i][j][k]-D->Pr[i][j-1][k]-D->Pl[i][j][k]+D->Pl[i][j-1][k])+dtOverdz*(D->Sr[i][j][k]-D->Sr[i][j][k-1]-D->Sl[i][j][k]+D->Sl[i][j][k-1])
  	      	  -2*piDt*D->Jx[i][j][k]+dF*dtOverdx*(D->F[i+1][j][k]-D->F[i][j][k]);;
  	      	D->Ex[i][j][k]=tmpd*(ExOld+tmp*tmpr);

  	      	tmp=dtOverdz*(D->Pr[i][j][k+1]-D->Pr[i][j][k]+D->Pl[i][j][k+1]-D->Pl[i][j][k])-dtOverdy*(D->Sr[i][j+1][k]-D->Sr[i][j][k]+D->Sl[i][j+1][k]-D->Sl[i][j][k]);
  	      	D->Bx[i][j][k]=tmpd*(BxOld+tmp*tmpr);

  	      	D->ExNow[i][j][k]=0.5*(ExOld+D->Ex[i][j][k]);
  	      	D->BxNow[i][j][k]=0.5*(BxOld+D->Bx[i][j][k]);

  	    	}
  	  	}
  	}
	 
}

void PS_solve3D_Split(Domain *D,int iteration)
{
  	int i,j,k,istart,iend,jstart,jend,kstart,kend;
  	double dtOverdy,dtOverdz,qtDtByDy,qtDtByDz,dbPiDt,piDt;
  	double rtr,rtd,ltr,ltd,upr,upd,dnr,dnd,dF;
  	double tmp,tmpr,tmpd,oldPr,oldPl,oldSr,oldSl;

  	int nTasks,myrank,rankY,rankZ;
  	MPI_Status status;
  	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  	istart=D->istart;    iend=D->iend;
  	jstart=D->jstart;    jend=D->jend;
  	kstart=D->kstart;    kend=D->kend;
  	dF=D->dF;

  	for(k=kstart; k<kend; k++) {
  	  D->frr[k]=1.0;
  	  D->frd[k]=1.0;
  	  D->bkr[k]=1.0;
  	  D->bkd[k]=1.0;
  	}
  	for(j=jstart; j<jend; j++) {
  	  D->upr[j]=1.0;
  	  D->upd[j]=1.0;
  	  D->dnr[j]=1.0;
  	  D->dnd[j]=1.0;
  	}
  	for(i=istart; i<iend; i++) {
  	  D->rtr[i]=1.0;
  	  D->rtd[i]=1.0;
  	  D->ltr[i]=1.0;
  	  D->ltd[i]=1.0;
  	}
  	if(D->pmlOn==ON && iteration>D->pmlStart) {
  	  absorb_FB3(D);
  	  absorb_UD3(D);
  	  absorb_RL3(D);
  	} else ;


  	dtOverdy=D->dt/D->dy;
  	dtOverdz=D->dt/D->dz;
  	piDt=M_PI*D->dt;
  	qtDtByDy=0.25*D->dt/D->dy;
  	qtDtByDz=0.25*D->dt/D->dz;
  	dbPiDt=2.0*M_PI*D->dt;
							
  	for(i=istart; i<iend; i++)
  	{
    	rtr=D->rtr[i]; ltr=D->ltr[i];
    	rtd=D->rtd[i]; ltd=D->ltd[i];
    	for(j=jstart; j<jend; j++)
    	{
      	upr=D->upr[j]; dnr=D->dnr[j];
      	upd=D->upd[j]; dnd=D->dnd[j];
      	for(k=kstart; k<kend; k++)
      	{
				oldPl=D->Pl[i][j][k];
				oldSl=D->Sl[i][j][k];
        		tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
        		tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];	

        		tmp=qtDtByDz*(D->Bx[i+1][j][k]+D->Bx[i][j][k]-D->Bx[i+1][j][k-1]-D->Bx[i][j][k-1])-qtDtByDy*(D->Ex[i+1][j+1][k]+D->Ex[i][j+1][k]-D->Ex[i+1][j][k]-D->Ex[i][j][k])
        		  -piDt*D->Jy[i+1][j][k]+dF*dtOverdy*(D->F[i+1][j+1][k]-D->F[i+1][j][k]);
        		D->Pl[i][j][k]=tmpd*(D->Pl[i+1][j][k]+tmp*tmpr);

        		tmp=-qtDtByDy*(D->Bx[i+1][j][k]+D->Bx[i][j][k]-D->Bx[i+1][j-1][k]-D->Bx[i][j-1][k])-qtDtByDz*(D->Ex[i+1][j][k+1]+D->Ex[i][j][k+1]-D->Ex[i+1][j][k]-D->Ex[i][j][k])
        		  -piDt*D->Jz[i+1][j][k]+dF*dtOverdz*(D->F[i+1][j][k+1]-D->F[i+1][j][k]);
        		D->Sl[i][j][k]=tmpd*(D->Sl[i+1][j][k]+tmp*tmpr);

				D->PlC[i][j][k]=0.5*(oldPl+D->Pl[i][j][k]);
				D->SlC[i][j][k]=0.5*(oldSl+D->Sl[i][j][k]);
      	}
    	}
  	}		
	
  	for(i=iend-1; i>=istart; i--)
  	{
    	rtr=D->rtr[i]; ltr=D->ltr[i];
    	rtd=D->rtd[i]; ltd=D->ltd[i];
    	for(j=jstart; j<jend; j++)
    	{
    	  	upr=D->upr[j]; dnr=D->dnr[j];
    	  	upd=D->upd[j]; dnd=D->dnd[j];
    	  	for(k=kstart; k<kend; k++)
    	  	{
				oldPr=D->Pr[i][j][k];
				oldSr=D->Sr[i][j][k];					
    	  	  	tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
    	  	  	tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];

    	  	  	tmp=qtDtByDz*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j][k-1]-D->Bx[i-1][j][k-1])+qtDtByDy*(D->Ex[i][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i][j][k]-D->Ex[i-1][j][k])
    	  	  	  -piDt*D->Jy[i][j][k]+dF*dtOverdy*(D->F[i][j+1][k]-D->F[i][j][k]);
    	  	  	D->Pr[i][j][k]=tmpd*(D->Pr[i-1][j][k]+tmp*tmpr);

    	  	  	tmp=-qtDtByDy*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j-1][k]-D->Bx[i-1][j-1][k])+qtDtByDz*(D->Ex[i][j][k+1]+D->Ex[i-1][j][k+1]-D->Ex[i][j][k]-D->Ex[i-1][j][k])
    	  	  	  -piDt*D->Jz[i][j][k]+dF*dtOverdz*(D->F[i][j][k+1]-D->F[i][j][k]);
    	  	  	D->Sr[i][j][k]=tmpd*(D->Sr[i-1][j][k]+tmp*tmpr);

				D->PrC[i][j][k]=0.5*(oldPr+D->Pr[i][j][k]);
				D->SrC[i][j][k]=0.5*(oldSr+D->Sr[i][j][k]);					  
    	  	}
    	}
  	}	

}

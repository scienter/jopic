#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>

void cleanMemory(Domain *D)
{
    Particle ***particle;
    particle=D->particle;
    LoadList *LL,*tmpLL;
    PlasmaLens *PL, *tmpPL;
    LaserList *L, *tmpL;
    int i,j,k,istart,iend,jstart,jend,kstart,kend;
    int s,nxSub,nySub,nzSub,start1D,end1D,start2D,end2D,start3D,end3D;
    double *plusX;
    void deleteField();
    ptclList *p,*tmp;
    DefPtcl *defP;
    IDCore *tr,*trTmp;	 
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    //remove particles
    start1D=0; end1D=iend+1;
    start2D=0; end2D=1;
    start3D=0; end3D=1;
    if(D->dimension>1) { start2D=jstart-1; end2D=jend+1; } else ;
    if(D->dimension>2) { start3D=kstart-1; end3D=kend+1; } else ;

    for(i=start1D; i<end1D; i++)
      for(j=start2D; j<end2D; j++)
        for(k=start3D; k<end3D; k++) {
          for(s=0; s<D->nSpecies; s++) {
            p=particle[i][j][k].head[s]->pt;
            while(p) {	
              tmp=p->next;
              particle[i][j][k].head[s]->pt=tmp; 
              p->next=NULL;
              free(p);
              p=particle[i][j][k].head[s]->pt;
            }
            free(particle[i][j][k].head[s]);
          }
          free(particle[i][j][k].head);
        }

    for(i=start1D; i<end1D; i++) {
      for(j=start2D; j<end2D; j++) {
        free(D->particle[i][j]);
      }
      free(D->particle[i]);
    }
    free(particle);

      
    LL=D->loadList;
    while(LL->next)
    {
      switch (LL->type)  {
      case Polygon :
        if(LL->xnodes>0)  {
          free(LL->xpoint);	
          free(LL->xn);	
          for(i=0; i<LL->xnodes; i++) free(LL->rpn_x[i]); 
          free(LL->expr_x);
          free(LL->rpn_size_x); 
          free(LL->rpn_x); 
          free(LL->expr_x_y0);
          free(LL->rpn_size_x_y0); 
          free(LL->rpn_x_y0); 
          free(LL->expr_x_z0);
          free(LL->rpn_size_x_z0); 
          free(LL->rpn_x_z0); 
        } else ;
        if(D->dimension>1)  {
          if(LL->ynodes>0)   {
            free(LL->ypoint);
            free(LL->yn);
          	for(i=0; i<LL->ynodes; i++) free(LL->rpn_y[i]); 
          	free(LL->expr_y);
          	free(LL->rpn_size_y); 
          	free(LL->rpn_y); 
          } else ;
        } else ;
        if(D->dimension>2) {
          if(LL->znodes>0)  {
            free(LL->zpoint);
            free(LL->zn);
          	for(i=0; i<LL->znodes; i++) free(LL->rpn_z[i]); 
          	free(LL->expr_z);
          	free(LL->rpn_size_z); 
          	free(LL->rpn_z); 
          } else ;
        } else ;
        break;
      case Defined :
//        if(LL->defineMode==byDensity) {
          defP=LL->def;
          while(defP) {
            LL->def=defP->next;
            defP->next=NULL;
            free(defP->define);
            free(defP);
            defP=LL->def;
          }
//        }
        break;
      case Beam :
        if(LL->xnodes>0 && LL->gaussMode==OFF) {
          free(LL->xpoint);	
          free(LL->xn);	
        } else ;
        break;
      }

      tmpLL=LL->next;
      D->loadList=tmpLL; 
      LL->next=NULL;
      free(LL);
      LL=D->loadList;
    }
    free(D->loadList);

    L=D->laserList;
    while(L)
    {	
      tmpL=L->next;
      D->laserList=tmpL; 
      L->next=NULL;
      free(L);
      L=D->laserList;
    }
    free(D->laserList);

	 
    // Plamsa Lens
    PL=D->lensList;
//	 if(D->nPlasmaLens>0) {
      while(PL) {
		  if(PL->xnodes>0) {
		    free(PL->xpoint);
          free(PL->xn);
		  } else ;
        tmpPL=PL->next;
        D->lensList=tmpPL;
        PL->next=NULL;
        free(PL);
        PL=D->lensList;
      }
//	 } else ;
    free(D->lensList);


    //remove field
    nxSub=D->nxSub+5;
    nySub=1;
    nzSub=1;
    if(D->dimension>1) nySub=D->nySub+5; else;
    if(D->dimension>2) nzSub=D->nzSub+5; else;

	 free(D->shareF);
   deleteField(D->F,nxSub,nySub,nzSub);
   deleteField(D->Rho,nxSub,nySub,nzSub);
    deleteField(D->RhoPair,nxSub,nySub,nzSub);
    deleteField(D->RhoNoPair,nxSub,nySub,nzSub);
    deleteField(D->Phi,nxSub,nySub,nzSub);
    deleteField(D->PhiOld,nxSub,nySub,nzSub);
    deleteField(D->Den,nxSub,nySub,nzSub);
    deleteField(D->DenOld,nxSub,nySub,nzSub);
    deleteField(D->CurX,nxSub,nySub,nzSub);
    deleteField(D->CurY,nxSub,nySub,nzSub);
    deleteField(D->CurZ,nxSub,nySub,nzSub);
    deleteField(D->Ax,nxSub,nySub,nzSub);
    deleteField(D->Ay,nxSub,nySub,nzSub);
    deleteField(D->Az,nxSub,nySub,nzSub);

    if(D->fieldType==Pukhov)
    {
      deleteField(D->Ex,nxSub,nySub,nzSub);
      deleteField(D->Ey,nxSub,nySub,nzSub);
      deleteField(D->Ez,nxSub,nySub,nzSub);
      deleteField(D->Bx,nxSub,nySub,nzSub);
      deleteField(D->By,nxSub,nySub,nzSub);
      deleteField(D->Bz,nxSub,nySub,nzSub);
      deleteField(D->BxNow,nxSub,nySub,nzSub);
      deleteField(D->ByNow,nxSub,nySub,nzSub);
      deleteField(D->BzNow,nxSub,nySub,nzSub);
      deleteField(D->Jx,nxSub,nySub,nzSub);
      deleteField(D->Jy,nxSub,nySub,nzSub);
      deleteField(D->Jz,nxSub,nySub,nzSub);
    }
    else if(D->fieldType==Split)
    {
      deleteField(D->Ex,nxSub,nySub,nzSub);
      deleteField(D->Bx,nxSub,nySub,nzSub);
      deleteField(D->ExNow,nxSub,nySub,nzSub);
      deleteField(D->BxNow,nxSub,nySub,nzSub);
      deleteField(D->Pr,nxSub,nySub,nzSub);
      deleteField(D->Pl,nxSub,nySub,nzSub);
      deleteField(D->Sr,nxSub,nySub,nzSub);
      deleteField(D->Sl,nxSub,nySub,nzSub);
      deleteField(D->PrC,nxSub,nySub,nzSub);
      deleteField(D->PlC,nxSub,nySub,nzSub);
      deleteField(D->SrC,nxSub,nySub,nzSub);
      deleteField(D->SlC,nxSub,nySub,nzSub);		
      deleteField(D->Jx,nxSub,nySub,nzSub);
      deleteField(D->Jy,nxSub,nySub,nzSub);
      deleteField(D->Jz,nxSub,nySub,nzSub);
    }

    //PML
    free(D->frr); free(D->bkr);
    free(D->upr); free(D->dnr);
    free(D->rtr); free(D->ltr);
    free(D->frd); free(D->bkd);
    free(D->upd); free(D->dnd);
    free(D->rtd); free(D->ltd);



    //remove track	
	free(D->trackParticle);
	free(D->trackID);
	free(D->trackCore);
	free(D->trackWe);
	if(D->idNums>0) {

      tr=D->trhead->pt;
		while(tr) {
			trTmp=tr->next;
			D->trhead->pt=trTmp;
			tr->next=NULL;
			free(tr);
			tr=D->trhead->pt;
		}
		free(D->trhead);
	} else ;

      
	for(i=start1D; i<end1D; i++)
		for(j=start2D; j<end2D; j++)
			for(k=start3D; k<end3D; k++) {
				for(s=0; s<1; s++) {
					p=D->track[i][j][k].head[s]->pt;
					while(p) {
						tmp=p->next;
						D->track[i][j][k].head[s]->pt=tmp;
						p->next=NULL;
						free(p);
						p=D->track[i][j][k].head[s]->pt;
					}
					free(D->track[i][j][k].head[s]);
				}
				free(D->track[i][j][k].head);
			}
	for(i=start1D; i<end1D; i++) {
		for(j=start2D; j<end2D; j++) {
			free(D->track[i][j]);
		}
		free(D->track[i]);
	}
	free(D->track);



/*
    if(D->tracking==ON)
    {
      if(D->idNums>0)
      {
        for(i=0; i<D->idNums; i++)
         free(D->track[i]);
        free(D->track);
        free(D->trackID);
        free(D->trackCore);
        free(D->trackS);
      }
    }

    //remove probe
    if(D->probeNum>0)
    {
      for(i=0; i<D->probeNum; i++)
       free(D->probe[i]);
      free(D->probe);
      free(D->probeX);
      free(D->probeY);
    }


    //remove boost field
    Boost **boost;
    boost=D->boost;

    for(i=0; i<D->nxSub+5; i++)
      free(boost[i]);
    free(boost);
*/
}

void deleteField(double ***field,int nx,int ny,int nz)
{
   int i,j,k;
   for(i=0; i<nx; i++)  
   {
     for(j=0; j<ny; j++)  
       free(field[i][j]);
     free(field[i]);
   }
   free(field);
}

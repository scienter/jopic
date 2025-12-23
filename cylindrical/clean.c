#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>

void deleteField(double ***field,int nx,int ny,int nz);

void cleanMemory(Domain *D)
{
   Particle **particle;
   particle=D->particle;
   LoadList *LL,*tmpLL;
   LaserList *L, *tmpL;
   int i,j,istart,iend,jstart,jend,numMode;
   int s,nxSub,nySub;
   double *plusX;
   void deleteField();
   ptclList *p,*tmp;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 
   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   numMode=D->numMode;

   //remove particles
   for(i=0; i<iend+1; i++)
      for(j=jstart-1; j<jend+1; j++)
      {
        for(s=0; s<D->nSpecies; s++)
        {
          p=particle[i][j].head[s]->pt;
          while(p) {	
            tmp=p->next;
            particle[i][j].head[s]->pt=tmp; 
            p->next=NULL;
            free(p);
            p=particle[i][j].head[s]->pt;
          }
          free(particle[i][j].head[s]);
        }
        free(particle[i][j].head);
      }

   for(i=0; i<iend+1; i++) 
      free(D->particle[i]);
   free(particle);
      
   LL=D->loadList;
   while(LL->next)
   {
      switch (LL->type)  {
      case Polygon :
        	if(LL->xnodes>0)
        	{
          	free(LL->xpoint);	
          	free(LL->xn);
          	for(i=0; i<LL->xnodes; i++) free(LL->rpn_x[i]);
          	free(LL->expr_x);
          	free(LL->rpn_size_x);
          	free(LL->rpn_x);	
          	free(LL->infix_x);	
        	}
        	if(LL->ynodes>0)
        	{
          	free(LL->ypoint);
          	free(LL->yn);
            for(i=0; i<LL->ynodes; i++) free(LL->rpn_y[i]);
            free(LL->expr_y);
            free(LL->rpn_size_y);
            free(LL->rpn_y);
            free(LL->infix_y);
        	}
        	break;
      case Defined :
        	if(LL->numDefined>0)
        	{
         	free(LL->xPosition);
            free(LL->yPosition);
            for(i=0; i<LL->numDefined; i++)
              	free(LL->define[i]);
        	}
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
    free(D->laserI); free(D->laserPhase);

    //remove field
    numMode=D->numMode;
    nxSub=D->nxSub+5;
    nySub=D->nySub+5;
	  free(D->shareF);
    switch (D->fieldType) {
    case Yee :
    case NoCherenkov :
      deleteField(D->EzR,numMode,nxSub,nySub);
      deleteField(D->ErR,numMode,nxSub,nySub);
      deleteField(D->EpR,numMode,nxSub,nySub);
      deleteField(D->EzI,numMode,nxSub,nySub);
      deleteField(D->ErI,numMode,nxSub,nySub);
      deleteField(D->EpI,numMode,nxSub,nySub);
      deleteField(D->BzR,numMode,nxSub,nySub);
      deleteField(D->BrR,numMode,nxSub,nySub);
      deleteField(D->BpR,numMode,nxSub,nySub);
      deleteField(D->BzI,numMode,nxSub,nySub);
      deleteField(D->BrI,numMode,nxSub,nySub);
      deleteField(D->BpI,numMode,nxSub,nySub);
      deleteField(D->BzNowR,numMode,nxSub,nySub);
      deleteField(D->BrNowR,numMode,nxSub,nySub);
      deleteField(D->BpNowR,numMode,nxSub,nySub);
      deleteField(D->BzNowI,numMode,nxSub,nySub);
      deleteField(D->BrNowI,numMode,nxSub,nySub);
      deleteField(D->BpNowI,numMode,nxSub,nySub);
      deleteField(D->JzR,numMode,nxSub,nySub);
      deleteField(D->JrR,numMode,nxSub,nySub);
      deleteField(D->JpR,numMode,nxSub,nySub);
      deleteField(D->JzI,numMode,nxSub,nySub);
      deleteField(D->JrI,numMode,nxSub,nySub);
      deleteField(D->JpI,numMode,nxSub,nySub);
      deleteField(D->RhoNoPairR,numMode,nxSub,nySub);
      deleteField(D->RhoNoPairI,numMode,nxSub,nySub);
      deleteField(D->RhoPairR,numMode,nxSub,nySub);
      deleteField(D->RhoPairI,numMode,nxSub,nySub);
      deleteField(D->FR,numMode,nxSub,nySub);
      deleteField(D->FI,numMode,nxSub,nySub);
      deleteField(D->CnR,numMode,nxSub,nySub);
      deleteField(D->CnI,numMode,nxSub,nySub);
      break;
    case Split :
      deleteField(D->PrR,numMode,nxSub,nySub);
      deleteField(D->PlR,numMode,nxSub,nySub);
      deleteField(D->SrR,numMode,nxSub,nySub);
      deleteField(D->SlR,numMode,nxSub,nySub);
      deleteField(D->EzR,numMode,nxSub,nySub);
      deleteField(D->BzR,numMode,nxSub,nySub);
      deleteField(D->PrI,numMode,nxSub,nySub);
      deleteField(D->PlI,numMode,nxSub,nySub);
      deleteField(D->SrI,numMode,nxSub,nySub);
      deleteField(D->SlI,numMode,nxSub,nySub);
      deleteField(D->EzI,numMode,nxSub,nySub);
      deleteField(D->BzI,numMode,nxSub,nySub);
      deleteField(D->EzNowR,numMode,nxSub,nySub);
      deleteField(D->EzNowI,numMode,nxSub,nySub);
      deleteField(D->BzNowR,numMode,nxSub,nySub);
      deleteField(D->BzNowI,numMode,nxSub,nySub);
      deleteField(D->JzR,numMode,nxSub,nySub);
      deleteField(D->JrR,numMode,nxSub,nySub);
      deleteField(D->JpR,numMode,nxSub,nySub);
      deleteField(D->JzI,numMode,nxSub,nySub);
      deleteField(D->JrI,numMode,nxSub,nySub);
      deleteField(D->JpI,numMode,nxSub,nySub);
      deleteField(D->RhoNoPairR,numMode,nxSub,nySub);
      deleteField(D->RhoNoPairI,numMode,nxSub,nySub);
      deleteField(D->RhoPairR,numMode,nxSub,nySub);
      deleteField(D->RhoPairI,numMode,nxSub,nySub);
      deleteField(D->FR,numMode,nxSub,nySub);
      deleteField(D->FI,numMode,nxSub,nySub);
      deleteField(D->CnR,numMode,nxSub,nySub);
      deleteField(D->CnI,numMode,nxSub,nySub);
      break;
    }

/*  
    //remove track
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
//lala
void tmpCleanMemory(Domain *D)
{
    Particle **particle;
    particle=D->particle;
    LoadList *LL,*tmpLL;
    LaserList *L, *tmpL;
    int i,j,istart,iend,jstart,jend,numMode;
    int s,nxSub,nySub;
    double *plusX;
    ptclList *p,*tmp;
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;

    //remove particles
    for(i=0; i<iend+1; i++)
      for(j=jstart-1; j<jend+1; j++)
      {
        for(s=0; s<D->nSpecies; s++)
        {
          p=particle[i][j].head[s]->pt;
          while(p) {	
            tmp=p->next;
            particle[i][j].head[s]->pt=tmp; 
            p->next=NULL;
            free(p);
            p=particle[i][j].head[s]->pt;
          }
          free(particle[i][j].head[s]);
        }
        free(particle[i][j].head);
      }

    for(i=0; i<iend+1; i++) 
      free(D->particle[i]);
    free(particle);

    //remove field
    numMode=D->numMode;
    nxSub=D->nxSub+5;
    nySub=D->nySub+5;
    switch (D->fieldType) {
    case Yee :
    case NoCherenkov :
      deleteField(D->EzR,numMode,nxSub,nySub);
      deleteField(D->ErR,numMode,nxSub,nySub);
      deleteField(D->EpR,numMode,nxSub,nySub);
      deleteField(D->EzI,numMode,nxSub,nySub);
      deleteField(D->ErI,numMode,nxSub,nySub);
      deleteField(D->EpI,numMode,nxSub,nySub);
      deleteField(D->BzR,numMode,nxSub,nySub);
      deleteField(D->BrR,numMode,nxSub,nySub);
      deleteField(D->BpR,numMode,nxSub,nySub);
      deleteField(D->BzI,numMode,nxSub,nySub);
      deleteField(D->BrI,numMode,nxSub,nySub);
      deleteField(D->BpI,numMode,nxSub,nySub);
      deleteField(D->BzNowR,numMode,nxSub,nySub);
      deleteField(D->BrNowR,numMode,nxSub,nySub);
      deleteField(D->BpNowR,numMode,nxSub,nySub);
      deleteField(D->BzNowI,numMode,nxSub,nySub);
      deleteField(D->BrNowI,numMode,nxSub,nySub);
      deleteField(D->BpNowI,numMode,nxSub,nySub);
      deleteField(D->JzR,numMode,nxSub,nySub);
      deleteField(D->JrR,numMode,nxSub,nySub);
      deleteField(D->JpR,numMode,nxSub,nySub);
      deleteField(D->JzI,numMode,nxSub,nySub);
      deleteField(D->JrI,numMode,nxSub,nySub);
      deleteField(D->JpI,numMode,nxSub,nySub);
      deleteField(D->RhoNoPairR,numMode,nxSub,nySub);
      deleteField(D->RhoNoPairI,numMode,nxSub,nySub);
      deleteField(D->RhoPairR,numMode,nxSub,nySub);
      deleteField(D->RhoPairI,numMode,nxSub,nySub);
      deleteField(D->FR,numMode,nxSub,nySub);
      deleteField(D->FI,numMode,nxSub,nySub);
      break;
    case Split :
      deleteField(D->PrR,numMode,nxSub,nySub);
      deleteField(D->PlR,numMode,nxSub,nySub);
      deleteField(D->SrR,numMode,nxSub,nySub);
      deleteField(D->SlR,numMode,nxSub,nySub);
      deleteField(D->EzR,numMode,nxSub,nySub);
      deleteField(D->BzR,numMode,nxSub,nySub);
      deleteField(D->PrI,numMode,nxSub,nySub);
      deleteField(D->PlI,numMode,nxSub,nySub);
      deleteField(D->SrI,numMode,nxSub,nySub);
      deleteField(D->SlI,numMode,nxSub,nySub);
      deleteField(D->EzI,numMode,nxSub,nySub);
      deleteField(D->BzI,numMode,nxSub,nySub);
      deleteField(D->EzNowR,numMode,nxSub,nySub);
      deleteField(D->BzNowR,numMode,nxSub,nySub);
      deleteField(D->EzNowI,numMode,nxSub,nySub);
      deleteField(D->BzNowI,numMode,nxSub,nySub);
      deleteField(D->JzR,numMode,nxSub,nySub);
      deleteField(D->JrR,numMode,nxSub,nySub);
      deleteField(D->JpR,numMode,nxSub,nySub);
      deleteField(D->JzI,numMode,nxSub,nySub);
      deleteField(D->JrI,numMode,nxSub,nySub);
      deleteField(D->JpI,numMode,nxSub,nySub);
      deleteField(D->RhoNoPairR,numMode,nxSub,nySub);
      deleteField(D->RhoNoPairI,numMode,nxSub,nySub);
      deleteField(D->RhoPairR,numMode,nxSub,nySub);
      deleteField(D->RhoPairI,numMode,nxSub,nySub);
      deleteField(D->FR,numMode,nxSub,nySub);
      deleteField(D->FI,numMode,nxSub,nySub);
      break;
    }

    //PML
    free(D->upr); free(D->upd);
    free(D->rtr); free(D->rtd);
    free(D->ltr); free(D->ltd);
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


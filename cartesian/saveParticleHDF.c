#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);

void saveParticleHDF(Domain *D,int iteration,int s,double minPx,double density)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;
    int cnt,totalCnt,index,start,dataCnt=9;
    int minXSub,minYSub,minZSub,nxSub,nySub,nzSub;
    char name[100],dataName[100];
    double *data,dx,dy,dz,lambda;
    int *recv,*offSetRank;
    Particle ***particle;
    particle=D->particle;
    ptclList *p;
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,plist_id,dset_id,attr_id,as_id;
    hid_t filespace,memspace;
    hsize_t dimsf[2],count[2],offset[2],block[2],stride[2],a_dims;
    herr_t ierr;

    nxSub=D->nxSub; nySub=D->nySub; nzSub=D->nzSub;
    lambda=D->lambda;
    dx=D->dx*lambda; dy=1; dz=1;
    if(D->dimension>1) dy=D->dy*lambda; else ;
    if(D->dimension>2) dz=D->dz*lambda; else ;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;
    minXSub=D->minXSub; minYSub=D->minYSub; minZSub=D->minZSub;

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    ierr=H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
    //create file
    sprintf(name,"%dParticle%d.h5",s,iteration);
    file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
    ierr=H5Pclose(plist_id);

    recv = (int *)malloc(nTasks*sizeof(int ));

    cnt=0;
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)          {
          p=particle[i][j][k].head[s]->pt;
          while(p)   {
            if(p->p1>=minPx)  cnt++; else;
            p=p->next;
          }
        }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&cnt,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);

    start=0;
    for(i=0; i<myrank; i++)        start+=recv[i];
    totalCnt=0;
    for(i=0; i<nTasks; i++)        totalCnt+=recv[i];
    //file space
    dimsf[0]=totalCnt;
    dimsf[1]=dataCnt;
    filespace=H5Screate_simple(2,dimsf,NULL);
    sprintf(dataName,"%d",s);
    dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

    if(totalCnt>0)
    {
      data = (double *)malloc(cnt*dataCnt*sizeof(double ));

      index=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)  {
            p=particle[i][j][k].head[s]->pt;
            while(p)    {
              if(p->p1>=minPx)  {
                data[index*dataCnt+0]=(p->x+i-istart+minXSub)*dx;
                data[index*dataCnt+1]=(p->y+j-jstart+minYSub)*dy;
                data[index*dataCnt+2]=(p->z+k-kstart+minZSub)*dz;
                data[index*dataCnt+3]=p->p1;
                data[index*dataCnt+4]=p->p2;
                data[index*dataCnt+5]=p->p3;
                data[index*dataCnt+6]=p->index;
                data[index*dataCnt+7]=p->core;
                data[index*dataCnt+8]=p->weight*density*dx*dy*dz;
//                data[index*dataCnt+9]=p->charge;
                index++;
              }	else; 
              p=p->next;
            }
          }

      //memory space
      dimsf[0]=cnt;
      dimsf[1]=dataCnt;
      memspace=H5Screate_simple(2,dimsf,NULL);
      stride[0]=1;
      stride[1]=1;
      count[0]=1;
      count[1]=1;

      //hyperslab in file space
      block[0]=cnt;
      block[1]=dataCnt;
      offset[0]=start;
      offset[1]=0;
      H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,stride,count,block);

      //hyperslab in memory space
      offset[0]=0;
      H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset,stride,count,block);
      
      plist_id=H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
      H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,data);
      H5Pclose(plist_id);
      H5Sclose(memspace);
    
      free(data);
    }	else	; 	//End of totalCnt>0
    MPI_Barrier(MPI_COMM_WORLD);
/*
    //write meta
    a_dims=1;
    as_id=H5Screate_simple(1,&a_dims,NULL);
    sprintf(dataName,"totalCnt");
//lala
    attr_id=H5Acreate2(dset_id,dataName,H5T_NATIVE_INT,as_id,H5P_DEFAULT,H5P_DEFAULT);
    H5Awrite(attr_id,H5T_NATIVE_INT,&totalCnt);
    H5Aclose(attr_id);
    H5Sclose(as_id);
*/
    H5Dclose(dset_id);

    H5Sclose(filespace);

    free(recv);
    H5Fclose(file_id);

    if(myrank==0)
      saveIntMeta(name,"totalCnt",&totalCnt,1);
    else        ;

}


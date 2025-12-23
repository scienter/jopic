#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void calParameter(int nx,int *istart,int *iend,int *nxSub,int rankX,int *biasX,int L);
void saveParticleComp_Double(double *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet);
void saveParticleComp_Int(int *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet);
void saveDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
double calHighInter(double y1,double y2,double y3,double resolX);
double calLowInter(double y1,double y2,double y3,int resolX);

void saveDump(Domain D,int iteration)
{
  void saveDumpFieldHDF();
  void saveDumpParticleHDF();
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  if(D.saveDumpMode==HDF)  {
    saveDumpFieldHDF(&D,iteration);
    saveDumpParticleHDF(&D,iteration);
    if(myrank==0)  {
      printf("dumpField%d.h5\n",iteration);
      printf("dumpParticle%d.h5\n",iteration);
    }   else	;
  }

}

void saveDumpParticleHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend;
    int cnt,totalCnt,index,start,cntList[D->nSpecies],dataCnt=13;
    int minXSub,minYSub,minZSub,nxSub,nySub,nzSub;
    char name[100],dataName[100];
    double *data;
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
    hsize_t dimsf[2],count[2],offset[2],block[2],stride[2],a_dims,metaDim[1];
    herr_t ierr;

    nxSub=D->nxSub; nySub=D->nySub; nzSub=D->nzSub;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;
    minXSub=D->minXSub; minYSub=D->minYSub; minZSub=D->minZSub;

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    ierr=H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
    //create file
    sprintf(name,"dumpParticle%d.h5",iteration);
    file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
    ierr=H5Pclose(plist_id);

    recv = (int *)malloc(nTasks*sizeof(int ));

    for(s=0; s<D->nSpecies; s++)
    {
      cnt=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)          {
            p=particle[i][j][k].head[s]->pt;
            while(p)   {
              cnt++;
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
      cntList[s]=totalCnt;

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
                data[index*dataCnt+0]=p->x+i-istart+minXSub;
                data[index*dataCnt+1]=p->y+j-jstart+minYSub;
                data[index*dataCnt+2]=p->z+k-kstart+minZSub;
                data[index*dataCnt+3]=p->p1;
                data[index*dataCnt+4]=p->p2;
                data[index*dataCnt+5]=p->p3;
                data[index*dataCnt+6]=p->index;
                data[index*dataCnt+7]=p->core;
                data[index*dataCnt+8]=p->weight;
                data[index*dataCnt+9]=p->charge;
                data[index*dataCnt+10]=p->oldX-istart+minXSub;
                data[index*dataCnt+11]=p->oldY-jstart+minYSub;
                data[index*dataCnt+12]=p->oldZ-kstart+minZSub;
                index++;
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
        H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
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
      sprintf(dataName,"%dtotalCnt",s);
      attr_id=H5Acreate2(dset_id,dataName,H5T_NATIVE_INT,as_id,H5P_DEFAULT,H5P_DEFAULT);
      H5Awrite(attr_id,H5T_NATIVE_INT,&totalCnt);
      H5Aclose(attr_id);
      H5Sclose(as_id);
*/
      H5Dclose(dset_id);

      H5Sclose(filespace);
    }	//End of nSpecies
    free(recv);
    H5Fclose(file_id);

    if(myrank==0)  {
      sprintf(dataName,"nSpecies");
      saveIntMeta(name,dataName,&D->nSpecies,1);
      for(s=0; s<D->nSpecies; s++)  {
        sprintf(dataName,"%dtotalCnt",s);
        saveIntMeta(name,dataName,&cntList[s],1);
      }
      saveIntMeta(name,"/nx",&D->nx,1);
      saveIntMeta(name,"/ny",&D->ny,1);
      saveIntMeta(name,"/nz",&D->nz,1);
      saveIntMeta(name,"/minXDomain",&D->minXDomain,1);
      saveIntMeta(name,"/minYDomain",&D->minYDomain,1);
      saveIntMeta(name,"/minZDomain",&D->minZDomain,1);
    }  else	;
}

void saveJDump(Domain D,int iteration)
{
  void saveJDumpHDF();
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  if(D.saveDumpMode==HDF)  {
    saveJDumpHDF(&D,iteration);
    if(myrank==0) 
      printf("resolJ%d.h5\n",iteration);
  } else;
}


void saveBDump(Domain D,int iteration)
{
  void saveBDumpHDF();
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  if(D.saveDumpMode==HDF)  {
    saveBDumpHDF(&D,iteration);
    if(myrank==0) 
      printf("resolB%d.h5\n",iteration);
  } else;
}

void saveEDump(Domain D,int iteration)
{
  void saveEDumpHDF();
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  if(D.saveDumpMode==HDF)  {
    saveEDumpHDF(&D,iteration);
    if(myrank==0) 
      printf("resolE%d.h5\n",iteration);
  } else;
}

void saveDumpFieldHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int biasX,biasY,biasZ,offSetY,rankX,rankY,rankZ;
    int nxSub,nySub,nzSub;
    int offset[3];
    char name[100],name2[100];

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;
    void saveFieldComp();

    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    sprintf(name,"dumpField%d.h5",iteration);
    if(myrank==0)
    {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else ;

    rankX=myrank/(D->M*D->N);
    rankZ=(myrank%(D->M*D->N))/D->M;
    rankY=(myrank%(D->M*D->N))%D->M;

    if(myrank==0)  {
      saveIntMeta(name,"/maxXDomain",&D->maxXDomain,1);
      saveIntMeta(name,"/minXDomain",&D->minXDomain,1);
      saveIntMeta(name,"/minYDomain",&D->minYDomain,1);
      saveIntMeta(name,"/minZDomain",&D->minZDomain,1);
      saveIntMeta(name,"/nSpecies",&D->nSpecies,1);
      saveIntMeta(name,"/nx",&D->nx,1);
      saveIntMeta(name,"/ny",&D->ny,1);
      saveIntMeta(name,"/nz",&D->nz,1);
    } else	;
    MPI_Barrier(MPI_COMM_WORLD);
    
    switch((D->fieldType-1)*3+D->dimension) {
    //1D
    case (Split-1)*3+1:
      nx=D->nx+5; 
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);
      ny=1;  nySub=1; jstart=0; jend=1;
      nz=1;  nzSub=1; kstart=0; kend=1;

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=0;
      offset[2]=0;
      
      saveFieldComp(D->Ex,name,"/Ex",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->Pr,name,"/Pr",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->Pl,name,"/Pl",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->Bx,name,"/Bx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->Sr,name,"/Sr",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->Sl,name,"/Sl",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->Jx,name,"/Jx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->Jy,name,"/Jy",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->Jz,name,"/Jz",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->F,name,"/F",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->RhoPair,name,"/RhoPair",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      break;

    //1D
    case (Yee-1)*3+1:
    case (Pukhov-1)*3+1:
      nx=D->nx+5;  ny=1;
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);
      nySub=1;  nzSub=1;

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=0;
      offset[2]=0;
      
      saveFieldComp(D->Ex,name,"/Ex",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Ey,name,"/Ey",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Ez,name,"/Ez",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Bx,name,"/Bx",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->By,name,"/By",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Bz,name,"/Bz",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Jx,name,"/Jx",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Jy,name,"/Jy",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Jz,name,"/Jz",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->RhoPair,name,"/RhoPair",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      break;

    //2D
    case (Yee-1)*3+2:
    case (Pukhov-1)*3+2:
      nx=D->nx+5;
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);
      ny=D->ny+5;
      calParameter(ny,&jstart,&jend,&nySub,rankY,&biasY,D->M);

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=(D->minYSub-D->minYDomain)+biasY;
      offset[2]=0;
      
      saveFieldComp(D->Ex,name,"/Ex",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Ey,name,"/Ey",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Ez,name,"/Ez",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Bx,name,"/Bx",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->By,name,"/By",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Bz,name,"/Bz",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Jx,name,"/Jx",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Jy,name,"/Jy",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Jz,name,"/Jz",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->RhoPair,name,"/RhoPair",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      break;
    case (Split-1)*3+2:
      nx=D->nx+5;
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);
      ny=D->ny+5;
      calParameter(ny,&jstart,&jend,&nySub,rankY,&biasY,D->M);

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=(D->minYSub-D->minYDomain)+biasY;
      offset[2]=0;
      
      saveFieldComp(D->Ex,name,"/Ex",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Bx,name,"/Bx",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->ExNow,name,"/ExNow",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->BxNow,name,"/BxNow",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Pr,name,"/Pr",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Pl,name,"/Pl",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Sr,name,"/Sr",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Sl,name,"/Sl",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Jx,name,"/Jx",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Jy,name,"/Jy",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Jz,name,"/Jz",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->F,name,"/F",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->RhoPair,name,"/RhoPair",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      break;

    //3D
    case (Yee-1)*3+3:
    case (Pukhov-1)*3+3:
      nx=D->nx+5;
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);
      ny=D->ny+5;
      calParameter(ny,&jstart,&jend,&nySub,rankY,&biasY,D->M);
      nz=D->nz+5;
      calParameter(nz,&kstart,&kend,&nzSub,rankZ,&biasZ,D->N);

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=(D->minYSub-D->minYDomain)+biasY;
      offset[2]=(D->minZSub-D->minZDomain)+biasZ;
      
      saveFieldComp(D->Ex,name,"/Ex",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Ey,name,"/Ey",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Ez,name,"/Ez",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Bx,name,"/Bx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->By,name,"/By",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Bz,name,"/Bz",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Jx,name,"/Jx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Jy,name,"/Jy",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Jz,name,"/Jz",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->RhoPair,name,"/RhoPair",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      break;
    case (Split-1)*3+3:
      nx=D->nx+5;
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);
      ny=D->ny+5;
      calParameter(ny,&jstart,&jend,&nySub,rankY,&biasY,D->M);
      nz=D->nz+5;
      calParameter(nz,&kstart,&kend,&nzSub,rankZ,&biasZ,D->N);

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=(D->minYSub-D->minYDomain)+biasY;
      offset[2]=(D->minZSub-D->minZDomain)+biasZ;
      
      saveFieldComp(D->Ex,name,"/Ex",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Bx,name,"/Bx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->ExNow,name,"/ExNow",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->BxNow,name,"/BxNow",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Pr,name,"/Pr",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Pl,name,"/Pl",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Sr,name,"/Sr",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Sl,name,"/Sl",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Jx,name,"/Jx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Jy,name,"/Jy",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Jz,name,"/Jz",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->F,name,"/F",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->RhoPair,name,"/RhoPair",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      break;
    }	//End of switch (dimension)
}


/*
void saveDumpParticleResolHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend;
    int cnt,totalCnt,index,start,cntList[D->nSpecies],dataCnt=10;
    int minXSub,minYSub,minZSub,nxSub,nySub,nzSub;
    char name[100],dataName[100];
    double *data,px,py,pz,old,dx;
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
    hsize_t dimsf[2],count[2],offset[2],block[2],stride[2],a_dims,metaDim[1];
    herr_t ierr;

    nxSub=D->nxSub; nySub=D->nySub; nzSub=D->nzSub;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;
    minXSub=D->minXSub;
    minYSub=D->minYSub;
    minZSub=D->minZSub;

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    ierr=H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
    //create file
    sprintf(name,"dumpResolParticle%d.h5",iteration);
    file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
    ierr=H5Pclose(plist_id);

    recv = (int *)malloc(nTasks*sizeof(int ));

    for(s=0; s<D->nSpecies; s++)
    {
      cnt=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)          {
            p=particle[i][j][k].head[s]->pt;
            while(p)   {
              cnt++;
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
      cntList[s]=totalCnt;

      //file space
      dimsf[0]=totalCnt;
      dimsf[1]=dataCnt;
      filespace=H5Screate_simple(2,dimsf,NULL);
      sprintf(dataName,"%d",s);
      dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

      dx=1.0/(double)(D->resolX);
      if(totalCnt>0)
      {
        data = (double *)malloc(cnt*dataCnt*sizeof(double ));

        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)  {
              p=particle[i][j][k].head[s]->pt;
              while(p)    {
                if(D->resolHigh==ON)  {
                  if(D->resolX>1)  {
                    old=p->oldX-istart+minXSub;
                    data[index*dataCnt+0]=(p->x+i-p->oldX)*dx+old;
                    old=p->oldY-jstart+minYSub;
                    data[index*dataCnt+1]=(p->y+j-p->oldY)*dx+old;
                    old=p->oldZ-kstart+minZSub;
                    data[index*dataCnt+2]=(p->z+k-p->oldZ)*dx+old;
                    px=calHighInter(p->p1Old2,p->p1Old1,p->p1,D->resolX);
                    py=calHighInter(p->p2Old2,p->p2Old1,p->p2,D->resolX);
                    pz=calHighInter(p->p3Old2,p->p3Old1,p->p3,D->resolX);            
                  } else {
                    px=p->p1;
                    py=p->p2;
                    pz=p->p3;
                    old=p->oldX-istart+minXSub;
                    data[index*dataCnt+0]=(p->x+i-p->oldX)*dx+old;
                    old=p->oldY-jstart+minYSub;
                    data[index*dataCnt+1]=(p->y+j-p->oldY)*dx+old;
                    old=p->oldZ-kstart+minZSub;
                    data[index*dataCnt+2]=(p->z+k-p->oldZ)*dx+old;
                  }
                }  else {
                  if(D->resolX>1)  {
                    px=calLowInter(p->p1Old2,p->p1Old1,p->p1,D->resolX);
                    py=calLowInter(p->p2Old2,p->p2Old1,p->p2,D->resolX);
                    pz=calLowInter(p->p3Old2,p->p3Old1,p->p3,D->resolX);
                    data[index*dataCnt+0]=p->x+i-istart+minXSub;
                    data[index*dataCnt+1]=p->y+j-jstart+minYSub;
                    data[index*dataCnt+2]=p->z+k-kstart+minZSub;
                  } else {
                    px=p->p1Old1;
                    py=p->p2Old1;
                    pz=p->p3Old1;
                    data[index*dataCnt+0]=p->oldX-istart+minXSub;
                    data[index*dataCnt+1]=p->oldY-jstart+minYSub;
                    data[index*dataCnt+2]=p->oldZ-kstart+minZSub;
                  }
                }
                data[index*dataCnt+3]=px;
                data[index*dataCnt+4]=py;
                data[index*dataCnt+5]=pz;
                data[index*dataCnt+6]=p->index;
                data[index*dataCnt+7]=p->core;
                data[index*dataCnt+8]=p->weight;
                data[index*dataCnt+9]=p->charge;
                index++;
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
        H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,data);
        H5Pclose(plist_id);
        H5Sclose(memspace);
      
        free(data);
      }	else	; 	//End of totalCnt>0
      MPI_Barrier(MPI_COMM_WORLD);
//
//      //write meta
//      a_dims=1;
//      as_id=H5Screate_simple(1,&a_dims,NULL);
//      sprintf(dataName,"%dtotalCnt",s);
//      attr_id=H5Acreate2(dset_id,dataName,H5T_NATIVE_INT,as_id,H5P_DEFAULT,H5P_DEFAULT);
//      H5Awrite(attr_id,H5T_NATIVE_INT,&totalCnt);
//      H5Aclose(attr_id);
//      H5Sclose(as_id);
//
      H5Dclose(dset_id);

      H5Sclose(filespace);
    }	//End of nSpecies
    free(recv);
    H5Fclose(file_id);

    if(myrank==0)  {
      sprintf(dataName,"nSpecies");
      saveIntMeta(name,dataName,&D->nSpecies,1);
      for(s=0; s<D->nSpecies; s++)  {
        sprintf(dataName,"%dtotalCnt",s);
        saveIntMeta(name,dataName,&cntList[s],1);
      }
      saveIntMeta(name,"/nx",&D->nx,1);
      saveIntMeta(name,"/ny",&D->ny,1);
      saveIntMeta(name,"/nz",&D->nz,1);
      saveIntMeta(name,"/maxXDomain",&D->maxXDomain,1);
      saveIntMeta(name,"/minXDomain",&D->minXDomain,1);
      saveIntMeta(name,"/minYDomain",&D->minYDomain,1);
      saveIntMeta(name,"/minZDomain",&D->minZDomain,1);
    }  else	;
}
*/

double calHighInter(double y1,double y2,double y3,double resolX)
{
   double a,b,c,dx,x,result;

   a=0.5*(y1+y3)-y2;
   b=2.0*y2-1.5*y1-0.5*y3;
   c=y1;
   dx=1.0/resolX;
   x=1.5+dx*0.5;
   result=a*x*x+b*x+c;

   return result;
}   

double calLowInter(double y1,double y2,double y3,int resolX)
{
   double a,b,c,dx,x,result;

   a=0.5*(y1+y3)-y2;
   b=2.0*y2-1.5*y1-0.5*y3;
   c=y1;
   x=2.0-0.5*(resolX-1);
   result=a*x*x+b*x+c;
//lala
   return result;
}   

void saveEDumpHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int biasX,biasY,biasZ,rankX,rankY,rankZ;
    int nxSub,nySub,nzSub;
    double x,y,z,y1,y2,y3,a,b,c,dx,dy,dz;
    int offset[3];
    char name[100],name2[100];
    int *saveInt;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    hsize_t metaDim[1];
    herr_t status;
    void saveFieldComp();

    nxSub=D->nxSub;  nySub=D->nySub;  nzSub=D->nzSub;
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    metaDim[0]=1;

    sprintf(name,"resolE%d.h5",iteration);
    if(myrank==0)
    {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else ;

    rankX=myrank/(D->M*D->N);
    rankZ=(myrank%(D->M*D->N))/D->M;
    rankY=(myrank%(D->M*D->N))%D->M;

    if(myrank==0)  {
      saveIntMeta(name,"/minXDomain",&D->minXDomain,1);
      saveIntMeta(name,"/minYDomain",&D->minYDomain,1);
      saveIntMeta(name,"/minZDomain",&D->minZDomain,1);
      saveIntMeta(name,"/nx",&D->nx,1);
      saveIntMeta(name,"/ny",&D->ny,1);
      saveIntMeta(name,"/nz",&D->nz,1);
      saveDoubleMeta(name,"/dtRatio",&D->dtRatio,1);
    } else	;
    MPI_Barrier(MPI_COMM_WORLD);
    nx=D->nx;        ny=D->ny;        nz=D->nz;

    switch((D->fieldType-1)*3+D->dimension) {
    //1D
    case (Pukhov-1)*3+1:
      nx=D->nx+5; 
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=0;    offset[2]=0;
      
      saveFieldComp(D->Ex,name,"/Ex",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->Ey,name,"/Ey",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->Ez,name,"/Ez",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      break;
    //2D
    case (Pukhov-1)*3+2:
      nx=D->nx+5;
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);
      ny=D->ny+5;
      calParameter(ny,&jstart,&jend,&nySub,rankY,&biasY,D->M);

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=(D->minYSub-D->minYDomain)+biasY;
      offset[2]=0;
      
      saveFieldComp(D->Ex,name,"/Ex",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Ey,name,"/Ey",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Ez,name,"/Ez",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      break;
    }
}

void saveJDumpHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int biasX,biasY,biasZ,rankX,rankY,rankZ;
    int nxSub,nySub,nzSub;
    double x,y,z,y1,y2,y3,a,b,c,dx,dy,dz;
    int offset[3];
    char name[100],name2[100];
    int *saveInt;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    hsize_t metaDim[1];
    herr_t status;
    void saveFieldComp();
    void saveIntMeta();

    nxSub=D->nxSub;  nySub=D->nySub;  nzSub=D->nzSub;
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;
    metaDim[0]=1;

    sprintf(name,"resolJ%d.h5",iteration);
    if(myrank==0)
    {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else ;

    rankX=myrank/(D->M*D->N);
    rankZ=(myrank%(D->M*D->N))/D->M;
    rankY=(myrank%(D->M*D->N))%D->M;

    if(myrank==0)  {
      saveIntMeta(name,"/minXDomain",&D->minXDomain,1);
      saveIntMeta(name,"/minYDomain",&D->minYDomain,1);
      saveIntMeta(name,"/minZDomain",&D->minZDomain,1);
      saveIntMeta(name,"/nx",&D->nx,1);
      saveIntMeta(name,"/ny",&D->ny,1);
      saveIntMeta(name,"/nz",&D->nz,1);
      saveDoubleMeta(name,"/dtRatio",&D->dtRatio,1);
    }	else	;
    MPI_Barrier(MPI_COMM_WORLD);
    nx=D->nx;        ny=D->ny;        nz=D->nz;

    switch(D->dimension) {
    //1D
    case 1 :
      nx=D->nx+5; 
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=0;    offset[2]=0;
      
      saveFieldComp(D->Jx,name,"/Jx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->Jy,name,"/Jy",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->Jz,name,"/Jz",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      break;
    //2D
    case 2:
      nx=D->nx+5;  ny=D->ny+5;
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);
      calParameter(ny,&jstart,&jend,&nySub,rankY,&biasY,D->M);

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=(D->minYSub-D->minYDomain)+biasY;
      offset[2]=0;
      
      saveFieldComp(D->Jx,name,"/Jx",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Jy,name,"/Jy",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Jz,name,"/Jz",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      break;
    }
}


void saveBDumpHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int biasX,biasY,biasZ,rankX,rankY,rankZ;
    int nxSub,nySub,nzSub;
    double x,y,z,y1,y2,y3,a,b,c,dx,dy,dz;
    int offset[3];
    char name[100],name2[100];
    int *saveInt;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    hsize_t metaDim[1];
    herr_t status;
    void saveFieldComp();
    void saveIntMeta();

    nxSub=D->nxSub;    nySub=D->nySub;    nzSub=D->nzSub;
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    metaDim[0]=1;

    sprintf(name,"resolB%d.h5",iteration);
    if(myrank==0)
    {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else ;

    rankX=myrank/(D->M*D->N);
    rankZ=(myrank%(D->M*D->N))/D->M;
    rankY=(myrank%(D->M*D->N))%D->M;

    if(myrank==0)  {
      saveIntMeta(name,"/minXDomain",&D->minXDomain,1);
      saveIntMeta(name,"/minYDomain",&D->minYDomain,1);
      saveIntMeta(name,"/minZDomain",&D->minZDomain,1);
      saveIntMeta(name,"/nx",&D->nx,1);
      saveIntMeta(name,"/ny",&D->ny,1);
      saveIntMeta(name,"/nz",&D->nz,1);
      saveDoubleMeta(name,"/dtRatio",&D->dtRatio,1);
    }	else	;
    MPI_Barrier(MPI_COMM_WORLD);
    nx=D->nx;  ny=D->ny;  nz=D->nz;

    switch((D->fieldType-1)*3+D->dimension) {
    //1D
    case (Pukhov-1)*3+1:
      nx=D->nx+5;
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=0;      offset[2]=0;
      
      saveFieldComp(D->Bx,name,"/Bx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->By,name,"/By",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      saveFieldComp(D->Bz,name,"/Bz",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      break;
    //2D
    case (Pukhov-1)*3+2:
      nx=D->nx+5;     ny=D->ny+5;
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);
      calParameter(ny,&jstart,&jend,&nySub,rankY,&biasY,D->M);

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=(D->minYSub-D->minYDomain)+biasY;
      offset[2]=0;
      
      saveFieldComp(D->Bx,name,"/Bx",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->By,name,"/By",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Bz,name,"/Bz",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      break;
    }
}

void saveParticleComp_Double(double *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet)
{
  int i,j,k;
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id,tic_id;
  herr_t status;
  hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
  hsize_t dimsf[1],count[1],offset[1];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);

  dimsf[0]=totalCnt;
  filespace=H5Screate_simple(1,dimsf,NULL);

  count[0]=cnt;
  offset[0]=offSet;
//  if(cnt>0)
    memspace=H5Screate_simple(1,count,NULL);

  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,data);
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

//  if(cnt>0)
    H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void saveParticleComp_Int(int *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet)
{
  int i,j,k;
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id,tic_id;
  herr_t status;
  hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
  hsize_t dimsf[1],count[1],offset[1];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);

  dimsf[0]=totalCnt;
  filespace=H5Screate_simple(1,dimsf,NULL);

  count[0]=cnt;
  offset[0]=offSet;
//  if(cnt>0)
    memspace=H5Screate_simple(1,count,NULL);

  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT,memspace,subfilespace,plist_id,data);
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);
 
//  if(cnt>0)
    H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void calParameter(int nx,int *istart,int *iend,int *nxSub,int rankX,int *biasX,int L)
{
  if(L==1) {
    *istart=0;
    *iend+=3;
    *biasX=0;
    *nxSub=nx;
  } else  {
    if(rankX==0)  {
      *istart=0;
      *nxSub+=2;
      *biasX=0;
    }  else if(rankX==L-1)  {
      *iend+=3;
      *nxSub+=3;
      *biasX=2;
    } else
      *biasX=2;
  }
}

void deleteFieldData(Domain *D,double ***field,int *nxSub,int *nySub,int *nzSub,int rankX,int rankY)
{
   int i,j,k;

   if(rankX==0)  *nxSub-=2;
   else if(rankX==D->L-1)  *nxSub-=3;
   if(rankY==0)  *nySub-=2;
   else if(rankY==D->M-1)  *nySub-=3;

   for(i=0; i<*nxSub+5; i++)
   {
     for(j=0; j<*nySub+5; j++)
       free(field[i][j]);
     free(field[i]);
   }
   free(field);
}


void saveFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet)
{
    int ii,i,j,k,start;
    double *field;
    FILE *out;
    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[3],count[3],offset[3];

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//    H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//    MPI_Barrier(MPI_COMM_WORLD);

    file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);

    dimsf[0]=ny;
    dimsf[1]=nx;
    dimsf[2]=nz;
    filespace=H5Screate_simple(3,dimsf,NULL);

    count[0]=nySub;
    count[1]=nxSub;
    count[2]=nzSub;
    offset[0]=offSet[1];
    offset[1]=offSet[0];
    offset[2]=offSet[2];
    memspace=H5Screate_simple(3,count,NULL);

    field = (double *)malloc(nxSub*nySub*nzSub*sizeof(double ));

    dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    subfilespace=H5Dget_space(dset_id);
    H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);

    start=0;
      for(j=jstart; j<jend; j++)
        for(i=istart; i<iend; i++)
        {
          for(k=kstart; k<kend; k++)
            field[start+k-kstart]=data[i][j][k];
          start+=nzSub;     
        } 

    plist_id=H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
//    H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
	 
    H5Pclose(plist_id);
    H5Sclose(subfilespace);
    H5Dclose(dset_id);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Fclose(file_id);
    free(field);
	 
}

void MPI_saveIntArray(int *data,char *fileName,char *dataName,int offSet)
{
    int i;
    FILE *out;
    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[1],count[1],offset[1];

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//    H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//    MPI_Barrier(MPI_COMM_WORLD);

    file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);

    dimsf[0]=nTasks;
    filespace=H5Screate_simple(1,dimsf,NULL);

    count[0]=1;
    offset[0]=offSet;
    memspace=H5Screate_simple(1,count,NULL);

    dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    subfilespace=H5Dget_space(dset_id);
    H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);

    plist_id=H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_INT,memspace,subfilespace,plist_id,data);
    H5Pclose(plist_id);
    H5Sclose(subfilespace);
    H5Dclose(dset_id);

    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Fclose(file_id);
}

void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}
      
void saveDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}
      


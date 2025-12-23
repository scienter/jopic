#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include <unistd.h>

void saveFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void saveCoordHDF(Domain *D,char *fileName);
void Efield_xdmf(int dimension,char *fileName,int nx,int ny,int nz);
void Bfield_xdmf(int dimension,char *fileName,int nx,int ny,int nz);

void saveEFieldHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int nxSub,nySub,nzSub;
    int offset[3];
    char name[100],dataName[100],fileName[100];
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;

    nxSub=D->nxSub; nySub=D->nySub; nzSub=D->nzSub;
    istart=D->istart; iend=D->iend;
    jstart=D->jstart; jend=D->jend;
    kstart=D->kstart; kend=D->kend;
    nx=D->nx; ny=D->ny; nz=D->nz;

    sprintf(name,"fieldE%d.h5",iteration);
    if(myrank==0)     {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
      saveIntMeta(name,"/nx",&nx,1);
      saveIntMeta(name,"/ny",&ny,1);
      saveIntMeta(name,"/nz",&nz,1);
    } else      ;
    MPI_Barrier(MPI_COMM_WORLD);
	 	 
    offset[0]=D->minXSub-D->minXDomain;
    offset[1]=D->minYSub-D->minYDomain;
    offset[2]=D->minZSub-D->minZDomain;

    switch(D->fieldType) {
    case Split:
      saveFieldComp(D->Ex,name,"/Ex",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
//		sleep(1);
      saveFieldComp(D->Pr,name,"/Pr",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
//		sleep(1);
      saveFieldComp(D->Pl,name,"/Pl",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
//		sleep(1);
//      saveFieldComp(D->F,name,"/F",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
//      MPI_Barrier(MPI_COMM_WORLD);
      if(myrank==0)   {
        saveCoordHDF(D,name);
        sprintf(fileName,"fieldE%d",iteration);
        Efield_xdmf(D->dimension,fileName,nx,ny,nz);
        printf("%s\n",name); 
      }	 else  ;
      MPI_Barrier(MPI_COMM_WORLD);

      break;
    case Yee:
    case Pukhov:
      saveFieldComp(D->Ex,name,"/Ex",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Ey,name,"/Ey",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Ez,name,"/Ez",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      //saveFieldComp(D->F,name,"/F",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);

      if(myrank==0)   {
        saveCoordHDF(D,name);
        sprintf(fileName,"fieldE%d",iteration);
        Efield_xdmf(D->dimension,fileName,nx,ny,nz);
        printf("%s\n",name); 
      }	 else  ;
      break;
    }   //End of switch (fieldType)

}

void saveBFieldHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int nxSub,nySub,nzSub;
    int offset[3];
    char name[100],dataName[100],fileName[100];
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;

    nxSub=D->nxSub; nySub=D->nySub; nzSub=D->nzSub;
    istart=D->istart; iend=D->iend;
    jstart=D->jstart; jend=D->jend;
    kstart=D->kstart; kend=D->kend;
    nx=D->nx; ny=D->ny; nz=D->nz;

    sprintf(name,"fieldB%d.h5",iteration);
    if(myrank==0)     {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
      saveIntMeta(name,"/nx",&nx,1);
      saveIntMeta(name,"/ny",&ny,1);
      saveIntMeta(name,"/nz",&nz,1);
    } else      ;
    MPI_Barrier(MPI_COMM_WORLD);

    offset[0]=D->minXSub-D->minXDomain;
    offset[1]=D->minYSub-D->minYDomain;
    offset[2]=D->minZSub-D->minZDomain;

    switch(D->fieldType) {

    case Split:
      saveFieldComp(D->Bx,name,"/Bx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Sr,name,"/Sr",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Sl,name,"/Sl",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);

      if(myrank==0)   {
        saveCoordHDF(D,name);
        sprintf(fileName,"fieldB%d",iteration);
        Efield_xdmf(D->dimension,fileName,nx,ny,nz);
        printf("%s\n",name); 
      }	 else  ;
      MPI_Barrier(MPI_COMM_WORLD);
      break;

    case Yee:
    case Pukhov:
      saveFieldComp(D->Bx,name,"/Bx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->By,name,"/By",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);
      saveFieldComp(D->Bz,name,"/Bz",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      MPI_Barrier(MPI_COMM_WORLD);

      if(myrank==0)   {
        saveCoordHDF(D,name);
        sprintf(fileName,"fieldB%d",iteration);
        Efield_xdmf(D->dimension,fileName,nx,ny,nz);
        printf("%s\n",name); 
      }	 else  ;
      break;
    }   //End of switch (dimension)
}

void Efield_xdmf(int dimension,char *fileName,int nx,int ny,int nz)
{
    FILE *xmf = 0;
    char name[100];
    const char *Names[] = {"/Ex","/Ey","/Ez"};
    int i;
    //
     // Open the file and write the XML description of the mesh..
     //
    sprintf(name,"%s.xmf",fileName);
    xmf = fopen(name,"w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");

    switch (dimension)  {
    case 2 :
      fprintf(xmf, "     <Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"%d %d\"/>\n", ny,nx);
      fprintf(xmf, "     <Geometry GeometryType=\"VXVY\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s.h5:/X\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ny);
      fprintf(xmf, "        %s.h5:/Y\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");
      for(i=0; i<3; i++)
      {
        fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",Names[i]);
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx,ny);
        fprintf(xmf, "        %s.h5:/%s\n",fileName,Names[i]);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;
    case 3 :
      //3D
      fprintf(xmf, "     <Topology TopologyType=\"3DRectMesh\" NumberOfElements=\"%d %d %d\"/>\n",ny,nx,nz);
      fprintf(xmf, "     <Geometry GeometryType=\"VXVYVZ\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nz);
      fprintf(xmf, "        %s.h5:/Z\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s.h5:/X\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ny);
      fprintf(xmf, "        %s.h5:/Y\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");
      for(i=0; i<3; i++)
      {
        fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",Names[i]);
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nz,nx,ny);
        fprintf(xmf, "        %s.h5:/%s\n",fileName,Names[i]);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;
    }
}

void Bfield_xdmf(int dimension,char *fileName,int nx,int ny,int nz)
{
    FILE *xmf = 0;
    char name[100];
    const char *Names[] = {"/Bx","/By","/Bz"};
    int i;
    //
     // Open the file and write the XML description of the mesh..
     //
    sprintf(name,"%s.xmf",fileName);
    xmf = fopen(name,"w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");

    switch (dimension)  {
    case 2 :
      fprintf(xmf, "     <Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"%d %d\"/>\n", ny,nx);
      fprintf(xmf, "     <Geometry GeometryType=\"VXVY\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s.h5:/X\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ny);
      fprintf(xmf, "        %s.h5:/Y\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");
      for(i=0; i<3; i++)
      {
        fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",Names[i]);
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx,ny);
        fprintf(xmf, "        %s.h5:/%s\n",fileName,Names[i]);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;
    case 3 :
      //3D
      fprintf(xmf, "     <Topology TopologyType=\"3DRectMesh\" NumberOfElements=\"%d %d %d\"/>\n",ny,nx,nz);
      fprintf(xmf, "     <Geometry GeometryType=\"VXVYVZ\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nz);
      fprintf(xmf, "        %s.h5:/Z\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s.h5:/X\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ny);
      fprintf(xmf, "        %s.h5:/Y\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");
      for(i=0; i<3; i++)
      {
        fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",Names[i]);
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nz,nx,ny);
        fprintf(xmf, "        %s.h5:/%s\n",fileName,Names[i]);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;
    }
}


// From x-y ascii data from xgrafix, reconstruct the density profile
#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "math.h"
#include "mpi.h"
#include "string.h"

void restoreFloatArray(char *fileName,char *dataName,double *data,int totalCnt);
void saveIntMeta(char *fileName,char *dataName,int *data);
void saveDataHDF3D(char *fileName,double ***ne,int numX,int numY,int numZ,double *xc,double *yc,double *zc);
void saveXmf(char *fileName,int numX,int numY,int numZ);
void restoreHyperSlab2D(char *fileName,char *dataName,double *data,int nx,int ny,int pick);
void restoreField2D(char *fileName,char *dataName,double **data,int nx,int ny);
void restoreFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz);

void main(int argc, char *argv[])
{
   FILE *in,*out;
   char fileName[100],dataName[100],fileType[100],name1[100],name2[100],**name,outFile[100];
   int i,j,k,n,nx,ny,nz,pick,mode,numMode,dimension,stX,stY,stZ;
   int initial,final,timeStep,step,numS,minIndex,totalStep,index;
   double pickY,*sum,rU,ex,ey,ez,bx,by,bz,x,divisionLambda,dx,dy;
   double cosP,sinP,angle,reangle,coef;
   double ***Ex,***Ey,*dataX,*dataY,*dataZ;

    int myrank, nTasks;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


   if(argc < 4)
   {
      printf("hdf_field mode dimension initial final timestep stX stY\n");
      printf("mode0 : filetype(field,density) dataName \n");
      printf("mode:1 2D => minIndex filetype dataName\n");
      printf("mode:1 1D (sum field) => initial final timeStep dx(4e-8) rU disionLambda\n");
      printf("mode2(sum) : fileType numS spc1 spc2 ...\n");
      printf("mode3(Ey) : fileType numS spc1 spc2 ...\n");
      printf("mode4(cyl) : fileType dataName angle\n");
      printf("mode5(Split) : outType (Ey,Ez,By,Bz)\n");
      exit(0);
   }
   mode=atoi(argv[1]);
   dimension=atoi(argv[2]);
   initial=atoi(argv[3]);
   final=atoi(argv[4]);
   timeStep=atoi(argv[5]);
   stX=atoi(argv[6]);
   stY=atoi(argv[7]);
   stZ=stY;

   switch (mode) {
   
   case 5 :
     sprintf(dataName,"%s",argv[8]);
     if(strstr(dataName,"Ey") || strstr(dataName,"Bz")) {
        sprintf(fileName,"fieldE%d.h5",initial);
        sprintf(name1,"Pr");
        sprintf(name2,"Pl");
     } else if(strstr(dataName,"Ez") || strstr(dataName,"By")) {
        sprintf(fileName,"fieldB%d.h5",initial);
        sprintf(name1,"Sr");
        sprintf(name2,"Sl");
     } else {
        printf("outType does not matched.\n");
        exit(0);
     }
     if(fopen(fileName,"r")==NULL)  {
       printf("%s is not exited.\n",fileName);
       exit(0);
     } else ;

     if(strstr(dataName,"Ey") || strstr(dataName,"Ez")) coef=1.0;
     else if(strstr(dataName,"By") || strstr(dataName,"Bz")) coef=-1.0;
    
     saveIntMeta(fileName,"nx",&nx);
     saveIntMeta(fileName,"ny",&ny);
     saveIntMeta(fileName,"nz",&nz);

     if(dimension==2) 
     {       
       Ex=(double ***)malloc(nx*sizeof(double **));      
       Ey=(double ***)malloc(nx*sizeof(double **));      
       for(i=0; i<nx; i++)  {
         Ex[i]=(double **)malloc(ny*sizeof(double *));      
         Ey[i]=(double **)malloc(ny*sizeof(double *));      
         for(j=0; j<ny; j++)  {
           Ex[i][j]=(double *)malloc(nz*sizeof(double ));      
           Ey[i][j]=(double *)malloc(nz*sizeof(double ));      
         }
       }
       dataX=(double *)malloc(nx*sizeof(double));      
       dataY=(double *)malloc(ny*sizeof(double));      

       for(step=initial; step<=final; step+=timeStep)
       {
         if(strstr(dataName,"Ey") || strstr(dataName,"Bz")) {
           sprintf(fileName,"fieldE%d.h5",step);
         } else if(strstr(dataName,"Ez") || strstr(dataName,"By")) {
           sprintf(fileName,"fieldB%d.h5",step);
         }
         if(fopen(fileName,"r")==NULL)  {
           printf("%s is not exited.\n",fileName);
           exit(0);
         } else ;

         restoreFloatArray(fileName,"X",dataX,nx);
         restoreFloatArray(fileName,"Y",dataY,ny);
         restoreFieldComp(Ex,fileName,name1,nx,ny,nz);
         restoreFieldComp(Ey,fileName,name2,nx,ny,nz);

         for(i=0; i<nx; i++)
           for(j=0; j<ny; j++)
             for(k=0; k<nz; k++)
               Ex[i][j][k]+=Ey[i][j][k]*coef;

         sprintf(outFile,"field%s%d",dataName,step);
         out=fopen(outFile,"w");
         k=0;
         for(i=0; i<nx; i+=stX)
         {
           for(j=0; j<ny; j+=stY)
             fprintf(out,"%g %g %g\n",dataX[i],dataY[j],Ex[i][j][k]);
           fprintf(out,"\n");
         }
         fclose(out);

         printf("%s is made.\n",outFile);

       }

       for(i=0; i<nx; i++) {
         for(j=0; j<ny; j++) {
           free(Ex[i][j]);   
           free(Ey[i][j]);   
         }
         free(Ex[i]);
         free(Ey[i]);
       }
       free(Ex); free(Ey);
       free(dataX);   
       free(dataY);  
     }

     break;   
   
   
   
   case 4 :
     angle=atof(argv[10]);
     reangle=angle*3.14159265359/180.0;
     sprintf(fileType,"%s",argv[8]);     
     sprintf(fileName,"%s%d.h5",fileType,initial);
     if(fopen(fileName,"r")==NULL)  {
       printf("%s is not exited.\n",fileName);
       exit(0);
     } else ;
     saveIntMeta(fileName,"nx",&nx);
     saveIntMeta(fileName,"ny",&ny);
     saveIntMeta(fileName,"numMode",&numMode);

     Ex=(double ***)malloc(nx*sizeof(double **));      
     Ey=(double ***)malloc(nx*sizeof(double **));      
     for(i=0; i<nx; i++)  {
       Ex[i]=(double **)malloc(ny*sizeof(double *));      
       Ey[i]=(double **)malloc(ny*sizeof(double *));      
       for(j=0; j<ny; j++)  {
         Ex[i][j]=(double *)malloc(numMode*sizeof(double ));      
         Ey[i][j]=(double *)malloc(numMode*sizeof(double ));      
       }
     }
     dataX=(double *)malloc(nx*sizeof(double));      
     dataY=(double *)malloc(ny*sizeof(double));      

     for(step=initial; step<=final; step+=timeStep)
     {
       n=0;
       sprintf(fileName,"%s%d.h5",fileType,step);
       if(fopen(fileName,"r")==NULL)  {
         printf("%s is not exited.\n",fileName);
         exit(0);
       } else ;
       restoreFloatArray(fileName,"X",dataX,nx);
       restoreFloatArray(fileName,"Y",dataY,ny);
       sprintf(dataName,"%sR",argv[9]);
       restoreFieldComp(Ex,fileName,dataName,nx,ny,numMode);
       sprintf(dataName,"%sI",argv[9]);
       restoreFieldComp(Ey,fileName,dataName,nx,ny,numMode);
       printf("%s is done\n",fileName);

       sprintf(fileName,"%s%s%d_%d_%g",argv[9],fileType,step,0,angle);
       out=fopen(fileName,"w");
       for(i=0; i<nx; i+=stX)
       {
         for(j=0; j<ny; j+=stY)
           fprintf(out,"%g %g %g\n",dataX[i],dataY[j],Ex[i][j][0]);
         fprintf(out,"\n");
       }
       fclose(out);
       printf("%s is made.\n",fileName);

       for(k=1; k<numMode; k++) {
         cosP=cos(k*reangle);
         sinP=sin(k*reangle);
         sprintf(fileName,"%s%s%d_%d_%g",argv[9],fileType,step,k,angle);
         out=fopen(fileName,"w");
         for(i=0; i<nx; i+=stX)  {
           for(j=0; j<ny; j+=stY) 
             fprintf(out,"%g %g %g\n",dataX[i],dataY[j]
                        ,Ex[i][j][k]*cosP-Ey[i][j][k]*sinP);
           fprintf(out,"\n");
         }
         fclose(out);
         printf("%s is made.\n",fileName);
       }

       for(i=0; i<nx; i+=stX) 
         for(j=0; j<ny; j+=stY) 
           for(k=1; k<numMode; k++) {
             cosP=cos(k*angle);
             sinP=sin(k*angle);
             Ex[i][j][0]+=Ex[i][j][k]*cosP-Ey[i][j][k]*sinP;
           }
       sprintf(fileName,"%s%s%d_%s_%g",argv[9],fileType,step,"sum",angle);
       out=fopen(fileName,"w");
       for(i=0; i<nx; i+=stX)
       {
         for(j=0; j<ny; j+=stY)
           fprintf(out,"%g %g %g\n",dataX[i],dataY[j],Ex[i][j][0]);
         fprintf(out,"\n");
       }
       fclose(out);
       printf("%s is made.\n",fileName);
     }

     for(i=0; i<nx; i++) {
       for(j=0; j<ny; j++) {
         free(Ex[i][j]);   
         free(Ey[i][j]);   
       }
       free(Ex[i]);
       free(Ey[i]);
     }
     free(Ex); free(Ey);
     free(dataX);   
     free(dataY);  
     break;

   case 0 :
     sprintf(fileName,"%s%d.h5",argv[8],initial);
     if(fopen(fileName,"r")==NULL)  {
       printf("%s is not exited.\n",fileName);
       exit(0);
     } else ;
     sprintf(dataName,"%s",argv[9]);
     saveIntMeta(fileName,"nx",&nx);
     saveIntMeta(fileName,"ny",&ny);
     saveIntMeta(fileName,"nz",&nz);

     Ex=(double ***)malloc(nx*sizeof(double **));      
     for(i=0; i<nx; i++)  {
       Ex[i]=(double **)malloc(ny*sizeof(double *));      
       for(j=0; j<ny; j++)  
         Ex[i][j]=(double *)malloc(nz*sizeof(double ));      
     }

     dataX=(double *)malloc(nx*sizeof(double));      
     dataY=(double *)malloc(ny*sizeof(double));      
     dataZ=(double *)malloc(nz*sizeof(double));      

     for(step=initial; step<=final; step+=timeStep)
     {
       sprintf(fileName,"%s%d.h5",argv[8],step);
       if(fopen(fileName,"r")==NULL)  {
         printf("%s is not exited.\n",fileName);
         exit(0);
       } else ;

       restoreFloatArray(fileName,"X",dataX,nx);
       if(dimension>1) restoreFloatArray(fileName,"Y",dataY,ny);
//       if(dimension>2) restoreFloatArray(fileName,"Z",dataZ,nz);
       if(dimension>2) for(i=0; i<nz; i++) dataZ[i]=dataY[i];

       restoreFieldComp(Ex,fileName,dataName,nx,ny,nz);

       sprintf(fileName,"%s%d_%s_XY",argv[8],step,argv[9]);
       k=0;
       if(dimension==3)	k=nz/2;	
       else ;

       if(dimension>1) {
         //saving file in XY plane
         out=fopen(fileName,"w");
         for(i=0; i<nx; i+=stX)
         {
           for(j=0; j<ny; j+=stY)
             fprintf(out,"%g %g %g\n",dataX[i],dataY[j],Ex[i][j][k]);
           fprintf(out,"\n");
         }
         fclose(out);
         printf("%s is made.\n",fileName);
       } else ;

       if(dimension==3) { 
         j=ny/2;
         //saving file in XY plane
         sprintf(fileName,"%s%d_%s_XZ",argv[8],step,argv[9]);
         out=fopen(fileName,"w");
         for(i=0; i<nx; i+=stX)
         {
           for(k=0; k<nz; k+=stY)
             fprintf(out,"%g %g %g\n",dataX[i],dataZ[k],Ex[i][j][k]);
           fprintf(out,"\n");
         }
         fclose(out);
         printf("%s is made.\n",fileName);
       } else ;

//       //saving file in center
//       sprintf(fileName,"cen%s%d_%s",argv[8],step,argv[9]);
       
//       if(dimension==1) {       k=0; j=0; }
//       else if(dimension==2) {       k=0; j=ny/2; }
//       else if(dimension==3)	{ k=nz/2; j=ny/2; }
//       else ;
//       out=fopen(fileName,"w");
//       for(i=0; i<nx; i++)
//         fprintf(out,"%.10g %g\n",dataX[i],Ex[i][j][k]);
//       fclose(out);
//       printf("%s is made.\n",fileName);


     }

     for(i=0; i<nx; i++)
       for(j=0; j<ny; j++)
         free(Ex[i][j]);   
     for(i=0; i<nx; i++)
       free(Ex[i]);   
     free(Ex);
     free(dataX); free(dataY); free(dataZ);

     break;

   case 1 :
       minIndex=atoi(argv[8]);
       
       sprintf(fileName,"%s%d.h5",argv[9],initial);
       if(fopen(fileName,"r")==NULL)  {
         printf("%s is not exited.\n",fileName);
         exit(0);
       } else ;
       sprintf(dataName,"%s",argv[10]);
       saveIntMeta(fileName,"nx",&nx);
       saveIntMeta(fileName,"ny",&ny);
       saveIntMeta(fileName,"nz",&nz);

       Ex=(double ***)malloc(nx*sizeof(double **));      
       for(i=0; i<nx; i++)  {
         Ex[i]=(double **)malloc(ny*sizeof(double *));      
         for(j=0; j<ny; j++)  
           Ex[i][j]=(double *)malloc(nz*sizeof(double ));      
       }
       dataX=(double *)malloc(nx*sizeof(double));      
       dataY=(double *)malloc(ny*sizeof(double));      
	   totalStep=(final-initial)/timeStep+1;
       sum=(double *)malloc(totalStep*sizeof(double));      
	   for(i=0; i<totalStep; i++) sum[i]=0.0;
	   
       restoreFloatArray(fileName,"X",dataX,nx);
//       dx=dataX[stX]-dataX[0]; 
       if(dimension>1) restoreFloatArray(fileName,"Y",dataY,ny);
//       dy=dataY[stY]-dataY[0];

       for(step=initial; step<=final; step+=timeStep)
       {
         sprintf(fileName,"%s%d.h5",argv[9],step);
         if(fopen(fileName,"r")==NULL)  {
           printf("%s is not exited.\n",fileName);
           exit(0);
         } else ;
         sprintf(dataName,"%s",argv[10]);

         restoreFieldComp(Ex,fileName,dataName,nx,ny,nz);
		 index=(step-initial)/timeStep;
         for(i=nx-minIndex; i<nx; i+=stX)
           for(j=0; j<ny; j+=stY)
             sum[index]+=Ex[i][j][0]*Ex[i][j][0];
         printf("%s is done.\n",fileName);
       }

	   sprintf(name1,"laserEnergy");
	   out=fopen(name1,"w");
	   for(i=0; i<totalStep-1; i++)
         fprintf(out,"%d %g %g\n",i*timeStep+initial,sum[i],sum[i]-sum[i+1]);	    
	   fclose(out);
       printf("%s is made.\n",name1);

	    
       for(i=0; i<nx; i++)
         for(j=0; j<ny; j++)
           free(Ex[i][j]);   
       for(i=0; i<nx; i++)
         free(Ex[i]);   
       free(Ex);
       free(dataX); if(dimension>1) free(dataY);   
	   
     break;

   case 2 :
     sprintf(fileType,"%s",argv[8]);
     numS=atoi(argv[9]);
     name=(char **)malloc(numS*sizeof(char *));      
     for(n=0; n<numS; n++)
       name[n]=(char *)malloc(100*sizeof(char ));      
     
     for(n=0; n<numS; n++) {
       sprintf(name[n],"%s",argv[n+10]);      
     }
     
     sprintf(fileName,"%s%s%d.h5",name[0],fileType,initial);
     if(fopen(fileName,"r")==NULL)  {
       printf("%s is not exited.\n",fileName);
       exit(0);
     } else ;
     sprintf(dataName,"%s",name[0]);
     saveIntMeta(fileName,"nx",&nx);
     saveIntMeta(fileName,"ny",&ny);
     saveIntMeta(fileName,"nz",&nz);

     Ex=(double ***)malloc(nx*sizeof(double **));      
     Ey=(double ***)malloc(nx*sizeof(double **));      
     for(i=0; i<nx; i++)  {
       Ex[i]=(double **)malloc(ny*sizeof(double *));      
       Ey[i]=(double **)malloc(ny*sizeof(double *));      
       for(j=0; j<ny; j++)  {
         Ex[i][j]=(double *)malloc(nz*sizeof(double ));      
         Ey[i][j]=(double *)malloc(nz*sizeof(double ));      
       }
     }
     dataX=(double *)malloc(nx*sizeof(double));      
     dataY=(double *)malloc(ny*sizeof(double));      
     dataZ=(double *)malloc(nz*sizeof(double));      

     for(step=initial; step<=final; step+=timeStep)
     {
       n=0;
       sprintf(fileName,"%s%s%d.h5",name[n],fileType,step);
       if(fopen(fileName,"r")==NULL)  {
         printf("%s is not exited.\n",fileName);
         exit(0);
       } else ;
       sprintf(dataName,"%s",name[n]);
       restoreFloatArray(fileName,"X",dataX,nx);
       restoreFloatArray(fileName,"Y",dataY,ny);
       if(dimension==3)
         restoreFloatArray(fileName,"Z",dataZ,nz);
       restoreFieldComp(Ex,fileName,name[n],nx,ny,nz);
       printf("%s is done\n",name[n]);

       for(n=1; n<numS; n++) {
         sprintf(fileName,"%s%s%d.h5",name[n],fileType,step);
         if(fopen(fileName,"r")==NULL)  {
           printf("%s is not exited.\n",fileName);
           exit(0);
         } else ;
         restoreFieldComp(Ey,fileName,name[n],nx,ny,nz);
         for(i=0; i<nx; i++)
           for(j=0; j<ny; j++)
             for(k=0; k<nz; k++)
               Ex[i][j][k]+=Ey[i][j][k];
         printf("%s is done\n",name[n]);
       }
       if(dimension>1) {
         sprintf(fileName,"%s%d_XY",fileType,step);
         out=fopen(fileName,"w");
         k=nz/2;
         for(i=0; i<nx; i+=stX)  {
           for(j=0; j<ny; j+=stY)
             fprintf(out,"%g %g %g\n",dataX[i],dataY[j],Ex[i][j][k]);
           fprintf(out,"\n");
         }
         fclose(out);
         printf("%s is made.\n",fileName);
       } else ;
       if(dimension==3) {
         sprintf(fileName,"%s%d_XZ",fileType,step);
         out=fopen(fileName,"w");
         j=ny/2;
         for(i=0; i<nx; i+=stX)  {
           for(k=0; k<nz; k+=stZ)
             fprintf(out,"%g %g %g\n",dataX[i],dataZ[k],Ex[i][j][k]);
           fprintf(out,"\n");
         }
         fclose(out);
         printf("%s is made.\n",fileName);
       }
     }

     for(i=0; i<nx; i++) {
       for(j=0; j<ny; j++) {
         free(Ex[i][j]);   
         free(Ey[i][j]);   
       }
       free(Ex[i]);
       free(Ey[i]);
     }
     free(Ex); free(Ey);
     free(dataX); free(dataY); free(dataZ); 
     for(n=0; n<numS; n++)  free(name[n]);
     free(name); 
     break;

   case 3 :
     sprintf(fileType,"%s",argv[8]);
     numS=atoi(argv[9]);
     name=(char **)malloc(numS*sizeof(char *));      
     for(n=0; n<numS; n++)
       name[n]=(char *)malloc(100*sizeof(char ));      
     
     for(n=0; n<numS; n++) {
       sprintf(name[n],"%s",argv[n+10]);      
     }
     
     sprintf(fileName,"%s%d.h5",fileType,initial);
     if(fopen(fileName,"r")==NULL)  {
       printf("%s is not exited.\n",fileName);
       exit(0);
     } else ;
     sprintf(dataName,"%s",name[0]);
     saveIntMeta(fileName,"nx",&nx);
     saveIntMeta(fileName,"ny",&ny);
     saveIntMeta(fileName,"nz",&nz);

     if(dimension==2) 
     {       
       Ex=(double ***)malloc(nx*sizeof(double **));      
       Ey=(double ***)malloc(nx*sizeof(double **));      
       for(i=0; i<nx; i++)  {
         Ex[i]=(double **)malloc(ny*sizeof(double *));      
         Ey[i]=(double **)malloc(ny*sizeof(double *));      
         for(j=0; j<ny; j++)  {
           Ex[i][j]=(double *)malloc(nz*sizeof(double ));      
           Ey[i][j]=(double *)malloc(nz*sizeof(double ));      
         }
       }
       dataX=(double *)malloc(nx*sizeof(double));      
       dataY=(double *)malloc(ny*sizeof(double));      

       for(step=initial; step<=final; step+=timeStep)
       {
         n=0;
         sprintf(fileName,"%s%d.h5",fileType,step);
         if(fopen(fileName,"r")==NULL)  {
           printf("%s is not exited.\n",fileName);
           exit(0);
         } else ;
         sprintf(dataName,"%s",name[n]);
         restoreFloatArray(fileName,"X",dataX,nx);
         restoreFloatArray(fileName,"Y",dataY,ny);
         restoreFieldComp(Ex,fileName,name[n],nx,ny,nz);
         printf("%s is done\n",name[n]);

         for(n=1; n<numS; n++) {
           sprintf(fileName,"%s%d.h5",fileType,step);
           if(fopen(fileName,"r")==NULL)  {
             printf("%s is not exited.\n",fileName);
             exit(0);
           } else ;
           restoreFieldComp(Ey,fileName,name[n],nx,ny,nz);
           for(i=0; i<nx; i++)
             for(j=0; j<ny; j++)
               for(k=0; k<nz; k++)
                 Ex[i][j][k]+=Ey[i][j][k];
           printf("%s is done\n",name[n]);
         }

         sprintf(fileName,"%s%d_Ey",fileType,step);
         out=fopen(fileName,"w");
         k=0;
         for(i=0; i<nx; i+=stX)
         {
           for(j=0; j<ny; j+=stY)
             fprintf(out,"%g %g %g\n",dataX[i],dataY[j],Ex[i][j][k]);
           fprintf(out,"\n");
         }
         fclose(out);

         printf("%s is made.\n",fileName);

       }

       for(i=0; i<nx; i++) {
         for(j=0; j<ny; j++) {
           free(Ex[i][j]);   
           free(Ey[i][j]);   
         }
         free(Ex[i]);
         free(Ey[i]);
       }
       free(Ex); free(Ey);
       free(dataX);   
       free(dataY);  
       for(n=0; n<numS; n++)  free(name[n]);
       free(name); 

     }

     break;

   }

    MPI_Finalize();

}
//lala
void restoreFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz)
{
  int i,j,k,start;
  double *field;
  char name[100];
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[3],count[3],offset[3];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//  H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//  MPI_Barrier(MPI_COMM_WORLD);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[0]=ny;
  dimsf[1]=nx;
  dimsf[2]=nz;
  filespace=H5Screate_simple(3,dimsf,NULL);

  count[0]=ny;
  count[1]=nx;
  count[2]=nz;
  offset[0]=0;
  offset[1]=0;
  offset[2]=0;
  memspace=H5Screate_simple(3,count,NULL);

  field = (double *)malloc(nx*ny*nz*sizeof(double ));

  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);

  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
  start=0;
  for(j=0; j<ny; j++)
    for(i=0; i<nx; i++)
    {
      for(k=0; k<nz; k++)
        data[i][j][k]=field[start+k];
      start+=nz;
    }
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
  free(field);
}


/*
void restoreField2D(char *fileName,char *dataName,double **data,int nx,int ny)
{
  hid_t file_id,dset_id,dataspace,memspace,filespace;
  hsize_t dimsf[2],offset[2],count[2],offset_out[2],count_out[2];
  herr_t status;
  double *field;
  int start,i,j;

  file_id=H5Fopen(fileName,H5F_ACC_RDONLY,H5P_DEFAULT);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  dataspace=H5Dget_space(dset_id);

  //define memory dataspace
  dimsf[1]=nx;
  dimsf[0]=ny;
  filespace=H5Screate_simple(2,dimsf,NULL);
 
  //define hyperslab in dataset
  offset[0]=0;
  offset[1]=0;
  count[1]=nx;
  count[0]=ny;
  memspace=H5Screate_simple(2,count,NULL);

  field=(double *)malloc(nx*ny*sizeof(double));      


//  status=H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offset,NULL,count,NULL);

//  status=H5Dread(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

  //define memory hyperslab  
//  offset_out[0]=0;
//  offset_out[1]=0;
//  count_out[0]=nx;
//  count_out[1]=ny;
//  status=H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset_out,NULL,count_out,NULL);

  //read data from hyperslab in the file into the hyperslab in memory and display 
  status=H5Dread(dset_id,H5T_NATIVE_FLOAT,memspace,dataspace,H5P_DEFAULT,field);

  start=0;
  for(j=0; j<ny; j++)
  {
    for(i=0; i<nx; i++)
      data[i][j]=field[start+i];
    start+=nx;
  }
       
  free(field);
  H5Dclose(dset_id);
  H5Sclose(dataspace);
  H5Sclose(memspace);
  H5Fclose(file_id);
}
*/

void restoreFloatArray(char *fileName,char *dataName,double *data,int totalCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=totalCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDONLY,H5P_DEFAULT);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);

  filespace=H5Screate_simple(1,metaDim,NULL);
  status=H5Dread(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}


void saveIntMeta(char *fileName,char *dataName,int *data)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=1;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}


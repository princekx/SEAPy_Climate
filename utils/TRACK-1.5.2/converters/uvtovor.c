#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "drscdf.h"

#define  TOL  1.0e-4

#define PI    3.14159265358979323846
#define PI180 0.017453292519943295
#define RE    6.371e+6

#define MAXCHR 100

void main(int argc, char *argv[])

{

   int i,j, k;
   int ierr=0;
   int ndim;
   int nx1, ny1, nt1;
   int nx2, ny2, nt2;
   int nl1, nl2;
   int dtx1, dty1, dtx2, dty2;
   int dtt1, dtt2;
   int dl1, dl2;
   int i1, i2;


   char fn1[MAXCHR];
   char fn2[MAXCHR];
   char fno[MAXCHR];

   char name[16];
   char unit[40], type[8];
   char time[8], date[8], title[80];
   char desc[120];

   float toff;
   float x11,x12, y11,y12;
   float x21,x22, y21,y22;
   float t11,t12, t21,t22;
   float lev1, lev2;
   float dtime, dlng, dlat;
   float *lng=NULL, *lat=NULL;
   float *nlng=NULL, *nlat=NULL;
   float *uarr=NULL, *varr=NULL;
   float *vor=NULL;
   float *xmean=NULL;
   float *ymean=NULL;
   float *dvdln=NULL;
   float *dudlt=NULL;
   float *atime=NULL;
   float r2;

   double arg1, arg2;

   if(argc != 5){

       printf("***Usage***, uvtovor [input file U (no extension)] [input file V (no extension)] [output file VOR (no extension)] [Time Offset]\n\n");

       exit(1);

   }

   strcpy(fn1, argv[1]);
   strcpy(fn2, argv[2]);

   if((ierr = Aslun(20, fn1, 21, fn1, IDRS_READ)) != IDRS_SUCCESS){

      printf("***ERROR***, can't open DRS files for READ, error code %d\n", ierr);

      exit(1);

   }

   ierr = Cluvdb();
   ierr = Setname(" ", "ua", " ", " ", " ");
   ierr = Inqdict(20, IDRS_GETFIRSTVAR);

   ierr = Getname(desc, name, title, unit, date, time, type, &ndim);

   printf("The first data file should have the zonal wind component\n");
   printf("The nature of the data for the first file is:-\r\n"
          "Name: %s, \n Title: %s, \n Units: %s, \n Data Type: %s, \n Var. No. %d.\n\n",
           name, title, unit, type, ndim);


   ierr=Getedim(1, desc, name, title, unit, &dtx1, &nx1, &x11, &x12);
   ierr=Getedim(2, desc, name, title, unit, &dty1, &ny1, &y11, &y12);
   ierr=Getedim(3, desc, name, title, unit, &dtt1, &nt1, &t11, &t12);
   ierr=Getedim(4, desc, name, title, unit, &dl1, &nl1, &lev1, &lev1);


   lng = (float *)calloc(nx1, sizeof(float));
   if(!lng){
      printf("****ERROR****, assigning memory\n");
      exit(1);
   }

   nlng = (float *)calloc(nx1, sizeof(float));
   if(!nlng){
      printf("****ERROR****, assigning memory\n");
   }  


   if(dtx1 == IDRS_EQUALLY_SPACED){

      dlng = (x12 - x11) / (nx1 - 1);

      for(i=0; i<nx1; i++){
          *(lng+i) = x11 + i * dlng;

      }

   }

   else{
      ierr = Setname(" ", "longitude", " "," ", " ");
      Getcdim(1, desc, name, title, unit, &dtx1, nx1, lng, &nx1);

   }

   for(i=0; i<nx1-1; i++) *(nlng + i) = (lng[i] + lng[i+1]) * 0.5;
   *(nlng + nx1 - 1) = lng[nx1-1] + (360.0 - lng[nx1-1] + lng[0]) * 0.5;


   lat = (float *)calloc(ny1, sizeof(float));
   if(!lat){
      printf("****ERROR****, assigning memory\n");
      exit(1);
   } 

   nlat = (float *)calloc(ny1, sizeof(float));
   if(!nlat){
      printf("****ERROR****, assigning memory\n");
      exit(1);
   }


   if(dty1 == IDRS_EQUALLY_SPACED){

      dlat = (y12 - y11) / (ny1 - 1);

      for(i=0; i<ny1; i++){
          *(lat+i) = y11 + i * dlat;
      }

   }

   else{
      ierr = Setname(" ", "latitude", " "," ", " ");
      Getcdim(2, desc, name, title, unit, &dty1, ny1, lat, &ny1);

   }


   for(i=0; i<ny1-1; i++) *(nlat + i) = (lat[i] + lat[i+1]) * 0.5;


   if((ierr = Aslun(30, fn2, 31, fn2, IDRS_READ)) != IDRS_SUCCESS){

      printf("***ERROR***, can't open DRS files for READ, error code %d\n", ierr);
      exit(1);

   }

   ierr = Cluvdb();
   ierr = Setname(" ", "va", " ", " ", " ");
   ierr = Inqdict(30, IDRS_GETFIRSTVAR);

   ierr = Getname(desc, name, title, unit, date, time, type, &ndim);

   printf("The second data file should have the meridional wind component\n");
   printf("The nature of the data for the second file is:-\r\n"
          "Name: %s, \n Title: %s, \n Units: %s, \n Data Type: %s, \n Var. No. %d.\n\n",
           name, title, unit, type, ndim);

   ierr=Getedim(1, desc, name, title, unit, &dtx2, &nx2, &x21, &x22);
   ierr=Getedim(2, desc, name, title, unit, &dty2, &ny2, &y21, &y22);
   ierr=Getedim(3, desc, name, title, unit, &dtt2, &nt2, &t21, &t22);
   ierr=Getedim(4, desc, name, title, unit, &dl2, &nl2, &lev2, &lev2);


   if(((dtx1 != dtx2) || (nx1 != nx2)) || (fabs(x11 - x21) > TOL) || (fabs(x12 - x22) > TOL)){
       printf("Input files are incompatable.\n");
       free(lng); free(lat);
       free(nlng); free(nlat);
       exit(1);
   }
   if((dty1 != dty2) || (ny1 != ny2) || (fabs(y11 - y21) > TOL) || (fabs(y12 - y22) > TOL)){
       printf("Input files are incompatable.\n");
       free(lng); free(lat);
       free(nlng); free(nlat);
       exit(1);
   }
   if((dtt1 != dtt2) || (nt1 != nt2) || (fabs(t11 - t21) > TOL) || (fabs(t12 - t22) > TOL)){
       printf("Input files are incompatable.\n");
       free(lng); free(lat);
       free(nlng); free(nlat);
       exit(1);
   }

   if((nl1 != nl2) || (fabs(lev1 - lev1) > TOL) ){
       printf("Input files are incompatable.\n");
       free(lng); free(lat);
       free(nlng); free(nlat);
       exit(1);
   }



   uarr = (float *)calloc(nx1*ny1, sizeof(float));
   if(!uarr){
      printf("****ERROR****, assigning memory\n");
      exit(1);
   }

   varr = (float *)calloc(nx1*ny1, sizeof(float));
   if(!varr){
      printf("****ERROR****, assigning memory\n");
      exit(1);
   }

   vor = (float *)calloc(nx1*(ny1-1), sizeof(float));
   if(!vor){
      printf("****ERROR****, assigning memory\n");
      exit(1);
   }

   xmean = (float *)calloc(nx1*ny1, sizeof(float));
   if(!xmean){
      printf("****ERROR****, assigning memory\n");
      exit(1);
   }

   ymean = (float *)calloc(nx1*(ny1-1), sizeof(float));
   if(!ymean){
      printf("****ERROR****, assigning memory\n");
      exit(1);
   }

   dvdln = (float *)calloc(nx1*(ny1-1), sizeof(float));
   if(!dvdln){
      printf("****ERROR****, assigning memory\n");
      exit(1);
   }

   dudlt = (float *)calloc(nx1*(ny1-1), sizeof(float));
   if(!dudlt){
      printf("****ERROR****, assigning memory\n");
      exit(1);
   }

   atime = (float *)calloc(nt1, sizeof(float));
   if(!atime){
      printf("****ERROR****, assigning memory\n");
      exit(1);
   }

   dtime = (t12 - t11) / (nt1 - 1);

   sscanf(argv[4], "%f", &toff);

   for(i=0; i<nt1; i++) *(atime+i) = i * dtime + t11 + toff;

   strcpy(fno, argv[3]);

   if((ierr = Aslun(40, fno, 41, fno, IDRS_CREATE)) != IDRS_SUCCESS){

      printf("***ERROR***, can't open DRS files for CREATE, error code %d\n", ierr);
      exit(1);

   }

   ierr = Cluvdb();
   ierr = Setname(" ", "longitude", " ", " ", " ");
   ierr = Putvdim(40, nx1, nlng, &i1, &i2);

   ierr = Cluvdb();
   ierr = Setname(" ", "latitude", " ", " ", " ");
   ierr = Putvdim(40, ny1-1, nlat, &i1, &i2);

   ierr = Cluvdb();
   ierr = Setname(" ", "time", " ", " ", " ");
   ierr = Putvdim(40, nt1, atime, &i1, &i2);

   ierr = Cluvdb();
   ierr = Setname(" ", "level", " ", " ", " ");
   ierr = Putvdim(40, 1, &lev1, &i1, &i2);

   if((ierr = Aslun(50, "new_va", 51, "new_va", IDRS_CREATE)) != IDRS_SUCCESS){

      printf("***ERROR***, can't open DRS files for CREATE, error code %d\n", ierr);
      exit(1);

   }

   ierr = Cluvdb();
   ierr = Setname(" ", "longitude", " ", " ", " ");
   ierr = Putvdim(50, nx1, nlng, &i1, &i2);

   ierr = Cluvdb();
   ierr = Setname(" ", "latitude", " ", " ", " ");
   ierr = Putvdim(50, ny1, lat, &i1, &i2);

   ierr = Cluvdb();
   ierr = Setname(" ", "time", " ", " ", " ");
   ierr = Putvdim(50, nt1, atime, &i1, &i2);

   ierr = Cluvdb();
   ierr = Setname(" ", "level", " ", " ", " ");
   ierr = Putvdim(50, 1, &lev1, &i1, &i2);



   for(i=0; i<nt1; i++){


      r2 = i * dtime + t11;
      ierr=Cluvdb();
      ierr = Setname(" ","ua"," "," "," ");
      ierr = Setdim(1,"longitude"," ",nx1,x11,x12);
      ierr = Setdim(2,"latitude"," ",ny1,y11,y12);
      ierr = Setdim(3,"time"," ",1,r2,r2);
      ierr = Getdat(20, uarr, nx1*ny1*4);
      ierr=Cluvdb ();
      ierr = Setname(" ","va"," "," "," ");
      ierr = Setdim(1,"longitude"," ",nx1,x11,x12);
      ierr = Setdim(2,"latitude"," ",ny1,y11,y12);
      ierr = Setdim(3,"time"," ",1,r2,r2);
      ierr = Getdat(30, varr, nx1*ny1*4);
      for(j=0; j<ny1; j++){
         for(k=0; k<nx1-1; k++){
             xmean[k+j*nx1] = (uarr[k+j*nx1] + uarr[k+1+j*nx1]) * 0.5;

         }
      }

      for(j=0; j < ny1; j++) {

          xmean[nx1-1+j*nx1] = (uarr[nx1-1+j*nx1] + uarr[j*nx1]) * 0.5;
      }

      for(j=0; j<nx1; j++){
         for(k=0; k<ny1-1; k++){
             ymean[k*nx1+j] = (varr[k*nx1+j] + varr[(k+1)*nx1+j])*0.5;
         }

      }

      for(j=0; j<nx1-1; j++){
          dlng = (lng[j+1] - lng[j]) * PI180;
          for(k=0; k < ny1-1; k++){
             dvdln[k*nx1+j] = (ymean[k*nx1+j+1] - ymean[k*nx1+j]) / dlng;
          }

      }

      if(lng[0] < TOL) dlng = (360.0 - lng[nx1-1]) * PI180;
      else dlng = (lng[0] + 360.0 - lng[nx1-1]) * PI180;

      for(j=0; j < ny1-1; j++) {
          dvdln[nx1-1+j*nx1] = (ymean[j*nx1] - ymean[nx1-1+j*nx1]) / dlng;
      }


      for(j=0; j<nx1; j++){
         for(k=0; k < ny1-1; k++){
/*            dudlt[k*nx1 + j] = (xmean[k*nx1 + j] * cos(lat[k] * PI180) - xmean[(k+1)*nx1 + j] * cos(lat[k+1] * PI180)) / ((lat[k] - lat[k+1]) * PI180); */
            dudlt[k*nx1 + j] = (xmean[(k+1)*nx1 + j] * cos(lat[k+1] * PI180) - xmean[k*nx1 + j] * cos(lat[k] * PI180)) / ((lat[k+1] - lat[k]) * PI180);
          }


      }


      for(j=0; j<ny1-1; j++){

          for(k=0; k < nx1; k++){

              vor[j*nx1+k] = (dvdln[j*nx1+k] - dudlt[j*nx1+k]) / (RE * cos(nlat[j] * PI180));

          }

      }

      ierr=Cluvdb ();
      ierr = Setname(" ","vor","vorticity","s-1"," ");
      ierr = Setvdim(1," ","longitude"," "," ",nlng[0],nlng[nx1-1]);
      ierr = Setvdim(2," ","latitude"," "," ",nlat[0],nlat[ny1-2]);
      ierr = Setvdim(3," ","time"," "," ", atime[i], atime[i]);
      ierr = Setvdim(4," ","level"," "," ", lev1, lev1);
      ierr = Putdat(40, vor);

      if(dtx1 == IDRS_EQUALLY_SPACED){

         ierr=Cluvdb ();
         ierr = Setname(" ","va","northerly wind","m/s"," ");
         ierr = Setdim(1,"longitude"," ",nx1,lng[0],lng[nx1-1]);
         ierr = Setdim(2,"latitude"," ",ny1,lat[0],lat[ny1-1]);
         ierr = Setdim(3,"time"," ",1,atime[i], atime[i]);
         ierr = Setdim(4,"level","mb",1, lev1, lev1);
         ierr = Putdat(50, varr);

      }

      else{

         ierr=Cluvdb ();
         ierr = Setname(" ","va","northerly wind","m/s"," ");
         ierr = Setvdim(1," ","longitude"," "," ",lng[0],lng[nx1-1]);
         ierr = Setvdim(2," ","latitude"," "," ",lat[0],lat[ny1-1]);
         ierr = Setvdim(4," ","time"," "," ", atime[i], atime[i]);
         ierr = Setvdim(3," ","level"," "," ", lev1, lev1);
         ierr = Putdat(50, varr);

      }

    }

    ierr = Cllun(20);
    ierr = Cllun(30);
    ierr = Cllun(40);

    free(lng); free(nlng);
    free(lat); free(nlat);
    free(uarr); free(varr);
    free(vor);
    free(xmean);
    free(ymean);
    free(dvdln);
    free(dudlt);
    free(atime);

    return;

}

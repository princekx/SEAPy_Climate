#include <Stdio.h>
#include <Math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "grid.h"
#include "files_out.h"
#include "file_handle.h"
#include "pp.h"
#include "mem_er.h"
#include "m_values.h"
#include "utf.h"
#include "netcdf_info.h"
#include "geo_values.h"

#define  TOLCON       1.0e-20
#define  THRESH       0.0
#define  MDIFF        1.0e+6


/* function to extract fields from a file with a chosen sampling, applying
   user chosen transforms and filling missing value holes.                 */

extern GRID *gr;
extern float *abuf;
extern char *ihead[NUMLINES];
extern int nlines;
extern int form;
extern int eqsw;
extern int utfv;
extern int frnum;
extern int i_utf4x, i_utf4y;
extern int std_x, std_y;
extern PP *pph;

extern char *fext;

extern int iext;


int rf(FILE * ,FILE * , int , int );
void writef(FILE * , float * , int );
int powi(int , int );
void fill_holes(float * , float , float , int , int , int , int );
void interp_holes(float * , float , int , int , int , int );
void convert(FILE * , int , int , int );
void lat_interp(float * , float , float , float , float , int , int , int , int , int );
int missing(float, float, int );

int extract(FILE *fdat)

{

   int i, j, nf=0;
   int nn=0;
   int dim=0;
   int wf='a';
   int trans='n';
   int tr_type=0;
   int fr1=0, fri=0, frl=0;
   int ix=0, iy=0;
   int interp=0;
   int il1=0, il2=0;
   int isin=0;
   int ifr = 0;
   int iyy = 0;
   int inv_utf=0;
   int nframe=0;
   int imiss='n';
   int icmp=0;

   float lat1=0.0, lat2=0.0, diff=MDIFF, dd, dlat=0.0;
   float nmiss=0.0;

   float *aa=NULL, *a1=NULL, *a2=NULL;
   float *plvor=NULL;
   float plv=0.0;
   float dum=0.0;

   char tsub[MAXCHR];

   off_t chrnum=0, place1, place2;

   float tr_val=1.0;
   float mval=MVAL, mdval=MVAL;

   FILE *fsub=NULL;


   strncpy(tsub, EXTRACT, MAXCHR);
   if(iext) strcpy(strstr(tsub, EXTENSION), fext);

   nlines = 0;

   for(i=0; i < NUMLINES; i++){
       ihead[i] = (char *)calloc(MAXCHR, sizeof(char));
       mem_er((ihead[i] == NULL) ? 0 : 1, MAXCHR * sizeof(char));
   }


   if(form < 2 || form == 4) {ix = std_x; iy = std_y; dim = std_x * std_y;}
   else if (form == 2){
      if(utfv == 3) {ix = gr->ix; iy = (eqsw) ? gr->iy + 1 : gr->iy; dim=gr->ix * iy;}
      else if(utfv == 4) {ix = i_utf4x; iy = i_utf4y; dim = i_utf4x * i_utf4y;}

      if(!gr->h_inv) {inv_utf = 1; gr->h_inv= 1;}

      if(inv_utf){
         aa = (float *)calloc(dim, sizeof(float));
         mem_er((aa == NULL) ? 0 : 1, dim * sizeof(float));

      }

   }

   else if(form == 3) {ix = pph->lbnpt; iy = pph->lbrow; dim = pph->lbrow * pph->lbnpt;}

   printf("What start, rate and end frames are required, fr1, fri, frl\n\n");
   scanf("%d %d %d", &fr1, &fri, &frl);

   printf("Do you want all frames in data file or just those frames   \r\n"
          " selected or all frames at the specified sampling rate from\r\n"
          " the starting frame, 'a', 's', or 'i'.                     \n\n");
 

   scanf("\n");
   wf = getchar();

   if(wf == 'a' ){fr1 = 1; fri = 1; frl = frnum;}

   else if(wf == 's'){
      if(fri < 1) fri = 1;
      if(frl > frnum) frl = frnum;
   }

   else if(wf == 'i'){
      if(fri < 1) fri = 1;
      frl = frnum;
   }

   else {

      printf("***WARNING***, incorrect decriptor defaulting to 'a' (all).\n\n");

      wf = 'a';
      fr1 = 1; fri = 1; frl = frnum;

   }

   printf("***INFORMATION***, frame counters for computing time average or time average removal are:\r\n"
          "                   Start = %d, Interval = %d, Last = %d \n\n", fr1, fri, frl);


   printf("If there are missing data values, delineated by some extreme value,     \r\n"
          "do you want no hole filling      '0'                                    \r\n"
          "do you want the holes filled by, '1' simple extraploation, or           \r\n"
          "                                 '2' interpolation or                   \r\n"
          "                                 '3' are the holes between latitudes.   \r\n"
          "                                 '4' replace missing value.             \n\n");

   scanf("%d", &interp);
   if(interp < 1 || interp > 4) interp = 0;


   if(interp){

      printf("****WARNING****, the grid should not be modified, i.e. by adding/removing \r\n"
             "                 lines of longitude/latitude.                             \n\n");


      if(form == 4){

         if(((NETCDF_INFO *)fdat)->imiss){
             printf("****INFORMATION****, missing data value is %e\n\n", ((NETCDF_INFO *)fdat)->missing);
             mdval = ((NETCDF_INFO *)fdat)->missing;
             imiss = 'y';

         }

         else {
             printf("****WARNING****, no missing values, nothing to do.                        \r\n"
                    "                 do you still want to interpolate/extrapolate? 'y' or 'n'.\n\n");
             scanf("\n");
             if(getchar() == 'n'){
                interp = 0;
                imiss = 'n';
             }
             else {
                imiss = 'y';
                printf("What is the default missing data value?\n\n");
                scanf("%f", &mdval);
             }
         }

      }

      else {
      
         printf("Are there missing data values, delineated by some extreme value. 'y' or 'n'\n\n");
         scanf("\n");
         imiss=getchar();
         if(imiss == 'y'){
            printf("What is the default missing data value?\n\n");
            scanf("%f", &mdval);
         }

      }

      if(imiss == 'y'){

         printf("What is the missing data value to test against?\n\n");
         scanf("%f", &mval);

         printf("How do you want to compare with missing value?\r\n"
                "'0' -- equality.                              \r\n"
                "'1' -- less than.                             \r\n"
                "'2' -- greater than.                          \n\n");
         scanf("%d", &icmp);

      }

      if(interp == 3){
         printf("What are the two latitudes to interpolate between, l1, l2\n\n");
         scanf("%f %f", &lat1, &lat2);

/* find the nearest grid lines associated with chosen latitudes */

         diff = MDIFF;

         for(i=0; i < iy; i++){
             if((dd=fabs(*(gr->ygrid + i) - lat1)) < diff){diff = dd; il1 = i;}
         }

         lat1 = *(gr->ygrid + il1);

         diff = MDIFF;

         for(i=0; i < iy; i++){
             if((dd=fabs(*(gr->ygrid + i) - lat2)) < diff){diff = dd; il2 = i;}
         }

         lat2 = *(gr->ygrid + il2);


         if(il1 > il2){
            iyy = il1;
            il1 = il2; 
            il2 = iyy;
            diff = lat1;
            lat1 = lat2;
            lat2 = diff;

         }

         dlat = lat1 - lat2;

         printf("Do you want interpolation wrt to Sin of latitude, 'y', or 'n'\n\n");
         scanf("\n");
         if(getchar() == 'y'){
            lat1 = sin(lat1 * FP_PI);
            lat2 = sin(lat2 * FP_PI);
            dlat = lat1 - lat2;
            isin = 1;       
         }


      }

      else if(interp == 4) {
         printf("What new value is required?\n\n");
         scanf("%f", &nmiss);
      }



   }

   else mval = MVAL;


   printf("Do you want to apply a transformation, e.g. multiply by a constant, 'y' or 'n'\n\n");
   scanf("\n");
   if((trans = getchar()) == 'y'){

      trans = 1;

      printf("What form of transform do you want:              \r\n"
             "  '1' add a constant.          (f+a)             \r\n"
             "  '2' subtract a constant.     (f-a)             \r\n"
             "  '3' multiply by a constant.  (f*a)             \r\n"
             "  '4' divide by a constant.    (f/a)             \r\n"
             "  '5' power                    (f^a)             \r\n"
             "  '6' lower threshold.         (f<a)             \r\n"
             "  '7' upper threshold.         (f>a)             \r\n"
             "  '8' log_10                   (log_10(f))       \r\n"
             "  '9' log_e                    (log_e(f))        \r\n"
             " '10' square                   (f*f)             \r\n"
             " '11' reciprocal               (1/f)             \r\n"
             " '12' exponential              (exp(f))          \r\n"
             " '13' absolute value           (abs(f))          \r\n"
             " '14' remove planetary vorticity from absolute.  \n\n");

     scanf("%d", &tr_type);
     if(tr_type < 1 || tr_type > 14){

        printf("****WARNING****, incorrect specifier, no transform perormed.\n\n");

     }

     if(tr_type >= 1 && tr_type <= 7){
        printf("What constant value do you want to use for the transform.\n\n");
        scanf("%f", &tr_val);
        if(tr_type == 4 && fabs(tr_val) <= TOLCON){
           printf("****WARNING****, constant too small no division performed.\n\n");
           tr_val = 1.0;
        }

     }

     if(tr_type == 14){
       plvor = (float *)calloc(iy, sizeof(float));
       mem_er((plvor == NULL) ? 0 : 1, iy*sizeof(float));
       for(i=0; i < iy; i++) *(plvor + i) = OMEGA * sin(FP_PI * *(gr->ygrid + i));  
     } 

     if(!interp){

        if(form == 4){

           if(((NETCDF_INFO *)fdat)->imiss){
              printf("****INFORMATION****, missing data value is %e\n\n", ((NETCDF_INFO *)fdat)->missing);
              mval = ((NETCDF_INFO *)fdat)->missing;
              imiss = 'y';

           }

         }

         else {

            printf("Are there missing data values, delineated by some extreme value. 'y' or 'n'\n\n");

            scanf("\n");
            imiss=getchar();
         }

         if(imiss == 'y'){

            printf("What is the missing data value?\n\n");
            scanf("%f", &mval);

            printf("How do you want to compare with missing value?\r\n"
                   "'0' -- equality.                              \r\n"
                   "'1' -- less than.                             \r\n"
                   "'2' -- greater than.                          \n\n");
            scanf("%d", &icmp);

         }

      }

   }

   else trans = 0;


   printf("***INFORMATON***, extracting fields from file........\r\n\n"
          "                                                     \r\n\n"
          " Please Wait, this may take some time ..............   \n\n");




/* open file for writing the selected fields */

   if(form != 4) fsub = open_file(tsub, "w+");
   else fsub = (FILE *)nc_define((NETCDF_INFO *)fdat, tsub);

   nframe = 0;

/* position file pointer */

   while(nf <= frl){

        if(form != 4){

           if(!nf) {

              place1 = ftello(fdat); 
              nlines = 0;

              rf(fdat, NULL, dim, 1);

              if(fr1 == 1){
                 ++nf;
                 for(i=0; i<nlines; i++) fprintf(fsub, "%s", ihead[i]);
              }

 
              place2 = ftello(fdat);
              chrnum = place2 - place1;

              if(fr1 > 1){

                 fseeko(fdat, (fr1-2)*chrnum, ORIGIN);

                 nlines = 0;

                 rf(fdat, NULL, dim, 1);

                 for(i=0; i<nlines; i++) fprintf(fsub, "%s", ihead[i]);


                 nf = fr1;           

              }

              nf += fri;

 
           }

           else {

              fseeko(fdat, (fri-1)*chrnum, ORIGIN);


              if(rf(fdat, NULL, dim, 1)) break;

              for(i=0; i<nlines; i++) fprintf(fsub, "%s", ihead[i]);

              nf += fri;

           }

        }

        else {

           if(!nf) nf = fr1;

           ((NETCDF_INFO *)fdat)->iframe = nf - 1;
           
           if(rf(fdat, NULL, dim, 1)) break; 

           ((NETCDF_INFO *)fsub)->iframe = nf - 1;
           ((NETCDF_INFO *)fsub)->nframe = nframe;
      
           nf += fri;

        }

        ++nframe;

        if(interp){

/* if utf and field is N-S invert to S-N */

           if(inv_utf){

              for(i=0; i < iy; i++){

                  a1 = aa + i * ix;
                  a2 = abuf + (iy - i - 1) * ix;

                  memcpy(a1, a2, ix * sizeof(float));

               }

               memcpy(abuf, aa, dim * sizeof(float));

           }

           if(interp == 1) {fill_holes(abuf, mval, mdval, ix, iy, ifr, icmp); ifr = 1;}
           else if(interp == 2){interp_holes(abuf, mval, ix, iy, ifr, icmp); ifr = 1;}
           else if(interp == 3) lat_interp(abuf, dlat, lat1, lat2, mval, il1, il2, ix, isin, icmp);

           else if(interp == 4){
               for(i=0; i < dim; i++){
                   if(missing(abuf[i], mval, icmp)) abuf[i] = nmiss;
               }
           }

/* need to check if we need to invert back. */

           if(inv_utf){

              for(i=0; i < iy; i++){

                  a1 = aa + i * ix;
                  a2 = abuf + (iy - i - 1) * ix;

                  memcpy(a1, a2, ix * sizeof(float));

               }

               memcpy(abuf, aa, dim * sizeof(float));

           }


        }


        if(trans){


           switch(tr_type){
                case 1:
                  for(i=0; i < dim; i++){
                     if(!missing(*(abuf + i), mval, icmp)) *(abuf + i) += tr_val;
                  }
                  break;
                case 2:
                  for(i=0; i < dim; i++) {
                     if(!missing(*(abuf + i), mval, icmp)) *(abuf + i) -= tr_val;
                  }
                  break;
                case 3:
                  for(i=0; i < dim; i++) {
                     if(!missing(*(abuf + i), mval, icmp)) *(abuf + i) *= tr_val;
                  }
                  break;
                case 4:
                  for(i=0; i < dim; i++) {
                     if(!missing(*(abuf + i), mval, icmp)) *(abuf + i) /= tr_val;
                  }
                  break;
                case 5:
                  for(i=0; i < dim; i++) {
                     if(*(abuf + i) < 0.0 && fabs(modff(tr_val, &dum)) > 0.0){
                        printf("****WARNING****, negative argument not allowed for fractional power.\n\n");
                        exit(1);
                     }
                     if(!missing(*(abuf + i), mval, icmp)) *(abuf + i) = pow(*(abuf + i), tr_val);
                  }
                  break;
                case 6:
                  for(i=0; i < dim; i++) {
                      if(*(abuf + i) <= tr_val && !missing(*(abuf + i), mval, icmp)) *(abuf + i) = THRESH;
                  }
                  break;
                case 7:
                  for(i=0; i < dim; i++) {
                      if(*(abuf + i) >= tr_val && !missing(*(abuf + i), mval, icmp)) *(abuf + i) = THRESH;
                  }
                  break;
                case 8:
                  for(i=0; i < dim; i++) {
                     if(*(abuf + i) <= 0.0){
                        printf("****WARNING****, negative argument not allowed for log10 function.\n\n");
                        exit(1);
                     }
                     if(!missing(*(abuf + i), mval, icmp)) *(abuf + i) = log10(*(abuf + i));
                  }
                  break;
                case 9:
                  for(i=0; i < dim; i++) {
                     if(*(abuf + i) <= 0.0){
                        printf("****WARNING****, negative argument not allowed for log function.\n\n");
                        exit(1);
                     }
                     if(!missing(*(abuf + i), mval, icmp)) *(abuf + i) = log(*(abuf + i));
                  }
                  break;
                case 10:
                  for(i=0; i < dim; i++) {
                      if(!missing(*(abuf + i), mval, icmp)) *(abuf + i) = *(abuf + i) * *(abuf + i);
                  }
                  break;
                case 11:
                  for(i=0; i < dim; i++) {
                     if(fabs(*(abuf + i)) <= TOLCON){
                        printf("****WARNING****, argument close to zero for reciprocal.\n\n");
                     }
                     if(!missing(*(abuf + i), mval, icmp)) *(abuf + i) = 1.0 / *(abuf + i);
                  }
                  break;
                case 12:
                  for(i=0; i < dim; i++) {
                      if(!missing(*(abuf + i), mval, icmp)) *(abuf + i) = exp(*(abuf + i));
                  }
                  break;
                case 13:
                  for(i=0; i < dim; i++) {
                      if(!missing(*(abuf + i), mval, icmp)) *(abuf + i) = fabs(*(abuf + i));
                  }
                  break;
                case 14:
                  for(i=0; i < iy; i++){
                      plv = *(plvor + i);
                      for(j=0; j < ix; j++)
                          if(!missing(*(abuf + i), mval, icmp)) *(abuf + i * ix + j) -= plv;
                  } 

                  break;
           }


        }


        writef(fsub, abuf, dim);

        ++nn;

        if(nf > frl) break;

   }

   printf("\n");

   if(interp == 1) fill_holes(NULL, 0.0, 0.0, 0, 0, -1, 0); 
   else if(interp == 2) interp_holes(NULL, 0.0, 0, 0, -1, 0);
 
   printf("****INFORMATION****, number of fields extracted is %d\n\n", nn);

   printf("Do you want to perform a conversion to binary data, 'y' or 'n'\n\n");
   scanf("\n");
   if(getchar() == 'y') {

      if(form != 4) fseeko (fsub, (off_t)0, FSTART);
      convert(fsub, 1, 1, nn);

   }

   if(form != 4) close_file(fsub, tsub);
   else netcdf_close((NETCDF_INFO *)fsub);

   if(plvor) free(plvor);
   free(abuf);
   for(i=0; i<NUMLINES; i++) free(ihead[i]);

   if(inv_utf)free(aa);

   return nn;

}

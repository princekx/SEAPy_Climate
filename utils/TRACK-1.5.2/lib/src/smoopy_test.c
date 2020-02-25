#include <Stdio.h>
#include "reg_dat.h"
#include "bisp.h"

/* calling program to test the C-fortran version of smoopy */

int smoopy_c(float * , float * , float * , int * , int * , float * , float , int );
void smoopy_setup(int );

float *xgrid, *ygrid;
float xmn, ymn, xmx, ymx;
float *ap;

int x1u, x2u, y1u, y2u;
int ix, iy;
int delb=0, bs='n';

struct sp_dat cokn;

void main()

{

   int i, j, k, mx, my, mxy;
   int iopt, ier;
   int nx, ny;
   int idif=1;
   int n=20;

   float *tx, *ty, *c;
   float s[5] = {10., 0.5, 0.22, 0.1, 0.};
   float fp, b;
   float d1[2];
   float f1, f2;

/*   extern float bisp_(); */

   extern struct rspline rs;

   FILE *read;

   if(!(read = fopen("datain", "r"))){

      printf("***error*** opening file for read\n");
      exit(1);

   }

/* read in x grid data */

   fscanf(read, "%d", &mx);

   xgrid = (float * )calloc(mx, sizeof(float));
   if(xgrid == NULL) printf("***error*** assigning storage for x grid data\n");

   for(i=0; i < mx; i++) fscanf(read, "%f", xgrid+i);

/* read in y grid data */

   fscanf(read, "%d", &my);

   ygrid = (float * )calloc(my, sizeof(float));
   if(ygrid == NULL) printf("***error*** assigning storage for y grid data\n");

   for(i=0; i < my; i++) fscanf(read, "%f", ygrid+i);

   ix = mx;
   iy = my;
   x1u = 1;
   y1u = 1;
   x2u = ix;
   y2u = iy;
   xmn = *(xgrid + x1u - 1);
   ymn = *(ygrid + y1u - 1);
   xmx = *(xgrid + x2u - 1);
   ymx = *(ygrid + y2u - 1);

/* read in field data */

   mxy = mx * my;

   ap = (float * )calloc(mxy, sizeof(float));
   if(ap == NULL) printf("***error*** assigning storage for field\n");

   for(i=0; i < my; i++){

       for(j=0; j < mx; j++) fscanf(read, "%f", ap + i*mx + j);

   }

   if(fclose(read) != 0){

     printf("***error*** closing file after read\n");
     exit(1);

   }

   smoopy_setup(0);

/* assign storage for knots and coefficients */

   cokn.tx = (float * )calloc(rs.nmax, sizeof(float));
   cokn.ty = (float * )calloc(rs.nmax, sizeof(float));
   cokn.c = (float * )calloc(rs.mxy, sizeof(float));
   if(cokn.tx == NULL || cokn.ty == NULL || cokn.c == NULL)
     printf("***error*** assigning storage for knot data\n");


   iopt = 0;

   for(j=4; j < 5; j++){

       ier = smoopy_c(&fp, s[j], iopt);

       printf("%f  %d\n", s[j], ier);
       printf("%e \n", fp);
       for(i=0; i < nx; i++) printf("%f ", *(tx + i));
       printf("\n");
       for(i=0; i < ny; i++) printf("%f ", *(ty + i));
       printf("\n\n");
       for(i=0; i < mx; i++){
          for(k=0; k < my; k++){

             bisp_(&b,d1,tx, &nx, ty, &ny, c, &rs.ncof, xgrid+i, ygrid+k,&idif);
             printf("%f %f %f %f %f\n", *(xgrid+i), *(ygrid+k), b, d1[0], d1[1]); 
/*             printf("%f ", b);*/
          }
          printf("\n");
       }
/*       for(i=0; i < n; i++){
          f1 = *xgrid + ((float)i/(n-1))*(*(xgrid + mx -1) - (*xgrid));
          for(k=0; k < n; k++){
              f2 = *ygrid + ((float)k/(n-1))*(*(ygrid + my -1) - (*ygrid));
              bisp_(&b,d1,tx, &nx, ty, &ny, c, &rs.ncof, &f1, &f2,&idif);
              printf("%f %f %f\n",f1, f2, b);
          }
       }  
       printf("\n\n");  */

       iopt = 1;

   }


   return;


}

void MAIN_(){}

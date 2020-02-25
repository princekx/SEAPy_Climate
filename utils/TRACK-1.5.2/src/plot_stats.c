#include <Stdio.h>
#include <stdlib.h>
#include "statistic.h"
#include "mem_er.h"
#include "proj.h"
#include "grid.h"
#include "file_handle.h"

#define  AVEC   2        /* max. number of vector attributes */
#define  UNDFF  0.0
#define  MAXCHR 100

/* function to plot statistical results on contour maps using uniras */

void stat_type();
void extract_sdt(struct tot_stat * , float * , float * , int * ,int , int , PROJ , GRDP * , GRDP * );
struct tot_stat *read_stats(FILE * );

#ifdef  NOUNDERSCORE


void statspl(float *, float *, float *, float *, int *, int *, int *, int *,
             int *, int *, int *, int *, int *, float *, float *, int *,
             float *, float *, float *, float *, float *, float *, int *,
             int *, float *, int *, int *, float *);

#else

void statspl_(float *, float *, float *, float *, int *, int *, int *, int *,
             int *, int *, int *, int *, int *, float *, float *, int *,
             float *, float *, float *, float *, float *, float *, int *,
             int *, float *, int *, int *, float *);

#endif

extern GRID *gr, *gr1, *gr2;
extern CNTRY *cm, *cm1, *cm2;
extern int aniso;

extern int tl, gof;
extern float period;

extern float xmn, ymn, xmx, ymx;


void plot_stats(struct tot_stat *trsav, int ptp, int npx, int npy)

{

   int st=0, stt=0, sty=-1, nmp=0;
   int str=0, stor=sizeof(float);
   int i, mv=AVEC;
   int ipj;
   int *ity=NULL, *fty=NULL;
   int *kern=NULL;
   int c1=0, c2=0, ca=0, cb=0;
   int ad_stat='n';
   int th_ty=0;
   int strep='n';
   int ireg=0;

   float *sp=NULL, *smth=NULL;
   float *xd, *yd;
   float ct;
   float *zp=NULL;
   float *xp=NULL, *yp=NULL;
   float *grid=NULL, *tdgrd=NULL;
   float xa1, xa2, ya1, ya2;
   float fd_th, td_th;
   float fd_val=0.0, td_val=0.0;
   float f1_val=0.0, f2_val=0.0;


   float tset=0.0;

   struct tot_stat *st1=NULL, *st2=NULL;

   struct pt_stat *pst=NULL;

   char sfile[MAXCHR];

   FILE *sff=NULL;

   PROJ proj={NULL, NULL};
   
   GRDP g1={0, 0, 0.0, 0.0, 0.0, 0.0};
   GRDP g2={0, 0, 0.0, 0.0, 0.0, 0.0};
   GRDP g22={0, 0, 0.0, 0.0, 0.0, 0.0};

/* ############################### */

   if(tl == 'y'){

     printf("****INFORMATION****, grid has been translated.\n\n");

     printf("Do you want to translate the statistics by the same amount, 'y' or 'n'?\n\n");
     scanf("\n");

     if(getchar() == 'y'){

        if(gof > 0) tset = (*(gr->xgrid + gof) - *(gr->xgrid));
        else if(gof < 0) tset = -(*(gr->xgrid + gr->ix - 1) - *(gr->xgrid + gr->ix + gof - 1));


        for(i=0; i < trsav->ptnum; i++){

            pst = trsav->ptst + i;

            pst->xs += tset;

            if(pst->xs > *(gr->xgrid + gr->ix - 1)) pst->xs -= period;
            else if(pst->xs < 0.) pst->xs += period;

        }

     }

  }

/* ++++++++++++++++++++++++++++ */

   printf("Do you want to use densities from another source for mean \r\n"
          "attribute surpression, 'y' or 'n'. This can be useful for \r\n"
          "difference plots.                                         \n\n");
   scanf("\n");
   if((ad_stat = getchar()) == 'y'){

      printf("You can use 1 or 2 alternative statistics files.   \r\n"
             "If two distinct files are used the max function is \r\n"
             "used to determine which value to use.              \n\n");

      printf("What is the first file to load?\n\n");
      scanf("%s", sfile);

      sff = open_file(sfile, "r");
      st1 = read_stats(sff);
      close_file(sff, "r"); 

      printf("What is the second file to load?\n\n");
      scanf("%s", sfile);

      sff = open_file(sfile, "r");
      st2 = read_stats(sff);
      close_file(sff, "r");  

      if(st1->ptnum != trsav->ptnum || st2->ptnum != trsav->ptnum){
         printf("****ERROR****, number of estimation points must be the same for\r\n"
                "               each set of statistics. Exiting.                \n\n");
         exit(1);
      }

      printf("Do you want the feature and/or track density replaced in the\r\n"
             "statistics by the max values from the alternative source, or\r\n"
             "placed in the spare statistics slots                        \r\n"
             "Useful for difference plots. 'y', 'n' or 's'.               \n\n");
      scanf("\n");
      strep = getchar(); 

      if(strep == 's'){
        trsav->sm[15] = trsav->sm[4];
        trsav->kern[15] = trsav->kern[4];
        trsav->sm[16] = trsav->sm[7];
        trsav->kern[16] = trsav->kern[7];
      }

   }

   printf("Input feature density threshold.\n\n");
   scanf("%f", &fd_th);

   printf("Input track density threshold.\n\n");
   scanf("%f", &td_th);

   printf("Currently the mean attribute statistics are suppressed for  \r\n"
          "seperate feature density and track density thresholds.      \r\n"
          "Do you want to make the suppression for only feature density\r\n"
          "or only track density, 'y' or 'n'.                          \n\n");

   scanf("\n");
   if(getchar() == 'y'){

      printf("Feature or Track density? '1' or '2'\n");
      
      scanf("%d", &th_ty);

      if(th_ty == 1)
         printf("****INFORMATION****, using feature density thresholding only.\n\n");
      else
         printf("****INFORMATION****, using track density thresholding only.\n\n");

   }

   else {

      printf("****INFORMATION****, using standard density thresholding, feature and track.\n\n");

   }

   xp = (float *)calloc(trsav->ptnum, sizeof(float));
   mem_er((xp == NULL) ? 0 : 1, trsav->ptnum * sizeof(float));

   yp = (float *)calloc(trsav->ptnum, sizeof(float));
   mem_er((yp == NULL) ? 0 : 1, trsav->ptnum * sizeof(float));

   zp = (float *)calloc(AVEC * trsav->ptnum, sizeof(float));
   mem_er((zp == NULL) ? 0 : 1, AVEC * trsav->ptnum * sizeof(float));
   
   g1.prty = gr->prty;
   g1.prgr = gr->prgr;

   if(g1.prgr) {

      g1.alat = gr1->alat; g1.alng = gr1->alng;
      g1.sp1= gr1->sp1; g1.sp2 = gr1->sp2;
   }

   xa1 = trsav->xa1;
   xa2 = trsav->xa2;
   ya1 = trsav->ya1;
   ya2 = trsav->ya2;

   printf("Statistics are defined on the region X:- %f -- %f            \r\n"
          "                                     Y:- %f -- %f            \r\n"
          "The current region is                X:- %f -- %f            \r\n"
          "                                     Y:- %f -- %f            \n\n", 
                                         xa1, xa2, ya1, ya2, xmn, xmx, ymn, ymx);

   printf("Do you want to use the defined region '0' \r\n"
          "                   the current region '1' \r\n"
          "                 or choose the region '2' \n\n");
   scanf("%d", &ireg);

   if(ireg == 1) {xa1 = xmn; xa2 = xmx; ya1 = ymn; ya2 = ymx;}
   else if(ireg == 2){
      printf("What is the new region for X, xa1, xa2?\n\n");
      scanf("%f %f", &xa1, &xa2);
      printf("What is the new region for Y, ya1, ya2?\n\n");
      scanf("%f %f", &ya1, &ya2);
   }


   printf("do you wish to change the projection, 'y' or 'n'    \r\n\n"
          "this assumes that the statistics have been computed   \r\n"
          "on the sphere so that only forward conversions from   \r\n"
          "cylindrical Plate Carre are possible.                 \n\n");

   proj_report(gr->prty, gr->prgr);

   scanf("\n");
   if((ipj = getchar()) == 'y') {

       proj = proj_group(-1, -1);
       
       g2.prty = gr->prty;
       g2.prgr = gr->prgr;
       if(g2.prgr){
          g2.alat = gr1->alat; g2.alng = gr1->alng;
	  g2.sp1= gr1->sp1; g2.sp2 = gr1->sp2;
	  
	  if(gr2){
	     g22.prty = gr->prty;
             g22.prgr = gr->prgr;
	     g22.alat = gr2->alat; g22.alng = gr2->alng;
	  }
       }

   }
   
   


   if(ipj == 'y' && !(g2.prgr) && !(g1.prgr)){

     if(proj.prj1 != NULL){

       (*proj.prj1)(&xa1, 'x', 0);
       (*proj.prj1)(&ya1, 'y', 0);
       (*proj.prj1)(&xa2, 'x', 0);
       (*proj.prj1)(&ya2, 'y', 0);

     }

     if(g2.prty > 0){

       (*proj.prj2)(&xa1, 'x', 1);
       (*proj.prj2)(&ya1, 'y', 1);
       (*proj.prj2)(&xa2, 'x', 1);
       (*proj.prj2)(&ya2, 'y', 1);

     }

   }


   for(i=0; i < trsav->ptnum; i++){


       pst = trsav->ptst + i;

       if(ad_stat == 'y'){

          f1_val = (st1->ptst + i)->stat3;
          f2_val = (st2->ptst + i)->stat3;
          fd_val = (f1_val < f2_val) ? f1_val : f2_val;

          f1_val = (st1->ptst + i)->stat6;
          f2_val = (st2->ptst + i)->stat6;
          td_val = (f1_val < f2_val) ? f1_val : f2_val;

          if(strep == 'y') {pst->stat3 = fd_val; pst->stat6 = td_val;} 
          else if(strep == 's') {pst->spare1 = fd_val; pst->spare2 = td_val;}

       }

       else {fd_val = pst->stat3; td_val = pst->stat6;}

       if(th_ty == 1){

          if(fd_val < fd_th){
             (pst->stat1).mean = (pst->stat1).var = UNDFF;
             (pst->stat2).mean = (pst->stat2).var = UNDFF;
             (pst->stat7).xcomp = (pst->stat7).ycomp = UNDFF;
             pst->stat9 = UNDFF;
             pst->stat10 = UNDFF;
             (pst->stat11).xcomp = (pst->stat11).ycomp = UNDFF;
             pst->stat12 = UNDFF;
             pst->stat13 = UNDFF;
          }  

          if(fd_val < fd_th) pst->stat8 = UNDFF;

       }

       else if(th_ty == 2){

          if(td_val < td_th){
             (pst->stat1).mean = (pst->stat1).var = UNDFF;
             (pst->stat2).mean = (pst->stat2).var = UNDFF;
             (pst->stat7).xcomp = (pst->stat7).ycomp = UNDFF;
             pst->stat9 = UNDFF;
             pst->stat10 = UNDFF;
             (pst->stat11).xcomp = (pst->stat11).ycomp = UNDFF;
             pst->stat12 = UNDFF;
             pst->stat13 = UNDFF;
          }  

          if(td_val < td_th) pst->stat8 = UNDFF;

       }

       else {

          if(fd_val < fd_th){
             (pst->stat1).mean = (pst->stat1).var = UNDFF;
             (pst->stat2).mean = (pst->stat2).var = UNDFF;
             (pst->stat7).xcomp = (pst->stat7).ycomp = UNDFF;
             pst->stat9 = UNDFF;
             pst->stat10 = UNDFF;
             (pst->stat11).xcomp = (pst->stat11).ycomp = UNDFF;
             pst->stat12 = UNDFF;
             pst->stat13 = UNDFF;
          }  

          if(td_val < td_th) pst->stat8 = UNDFF; 

       }


       xd = xp + c1;
       yd = yp + c1;

       *xd = pst->xs;
       *yd = pst->ys;

       if(ipj == 'y'){

          pst->ptyp = 1;

          switch(g2.prgr){
               case 0:
                  proj_point(proj, xd, yd, &g1, &g2);                            
                  pst->ptyp = 0;
                  ++c1;
                  break;
               case 1: case 2:
                  if(proj_point(proj, xd, yd, &g1, &g2)) 
                     {++c1; pst->ptyp = 0;}
                  else if(gr2){
                     proj_point(proj, xd, yd, &g1, &g22);
                     ++c2;
                     xp[trsav->ptnum - c2] = *xd;
                     yp[trsav->ptnum - c2] = *yd;

                   }
                   break;
          }

       }

       else ++c1;

   }

   if(st1){
     free(st1->ptst);
     free(st2->ptst);
     free(st1);
     free(st2);
   }

   printf("***INFORMATION***, the current plotting grid is nx=%d by ny=%d grid points\n\n", npx, npy);
   printf("Do you wish to change the plotting grid dimensions, 'y' or 'n'\n");
   scanf("\n");
   if(getchar() == 'y'){

      printf("Input new plotting grid dimensions, nx and ny\n\n");
      scanf("%d %d", &npx, &npy);

   }

   grid = (float *)calloc(npx*npy, sizeof(float));
   mem_er((grid == NULL) ? 0 : 1, npx*npy * sizeof(float));

   tdgrd = (float *)calloc(npx*npy, sizeof(float));
   mem_er((grid == NULL) ? 0 : 1, npx*npy * sizeof(float));

   ct = (*(gr->xgrid + gr->ix - 1) - *(gr->xgrid)) / 4.0;

   if(ptp){

      printf("which set do you want plotted, those available are:-\r\n");
      stat_type();
      printf("input zero or negative integer to finish input\n");

      scanf("%d", &sty);

/* put data into array for transfer to plotting routine */

      while(sty > 0){

           if(sty > STNM)

              printf("***ERROR***, incorrect stat identifier, continuing\n");

           else {

              ++st;
              ++stt;

              if(!nmp) {

                 nmp = trsav->ptnum;

                 fty = (int * )malloc_initl(sizeof(int));
                 mem_er((fty == NULL) ? 0 : 1, sizeof(int));

                 if(sty == 9 || sty == 13){

                   nmp += trsav->ptnum;

                   sp = (float * )calloc(nmp, stor);  
                   mem_er((sp == NULL) ? 0 : 1, nmp * stor);

                   *fty = 2;

                 }

                 else{

                   sp = (float * )calloc(nmp, stor);  
                   mem_er((sp == NULL) ? 0 : 1, nmp * stor);

                   *fty = 1;

                 }

                 ity = (int * )malloc_initl(sizeof(int));
                 mem_er((ity == NULL) ? 0 : 1, sizeof(int));

                 *ity = sty;

                 smth = (float * )malloc_initl(sizeof(float));
                 mem_er((smth == NULL) ? 0 : 1, sizeof(float));

                 kern = (int * )malloc_initl(sizeof(int));
                 mem_er((kern == NULL) ? 0 : 1, sizeof(int));

                 if(sty == 9 || sty == 13) ++stt;
                 extract_sdt(trsav, sp, smth, kern, stt, sty, proj, &g1, &g2);

              }

              else {

                 nmp += trsav->ptnum;

                 fty = (int * )realloc_n(fty, st * sizeof(int));
                 mem_er((fty == NULL) ? 0 : 1, st * sizeof(int));


                 if(sty == 9 || sty == 13){

                   nmp += trsav->ptnum;
                   str = nmp * stor;

                   sp = (float * )realloc_n(sp, str);  
                   mem_er((sp == NULL) ? 0 : 1, str);

                   *(fty + st - 1) = 2;

                 }

                 else{

                   str = nmp * stor;

                   sp = (float * )realloc_n(sp, str);  
                   mem_er((sp == NULL) ? 0 : 1, str);

                   *(fty + st - 1) = 1;

                 }

                 ity = (int * )realloc_n(ity, st * sizeof(int));
                 mem_er((ity == NULL) ? 0 : 1, st * sizeof(int));

                *(ity + st - 1) = sty;

                 smth = (float * )realloc_n(smth, st * sizeof(float));
                 mem_er((smth == NULL) ? 0 : 1, st * sizeof(float));

                 kern = (int * )realloc_n(kern, st * sizeof(int));
                 mem_er((kern == NULL) ? 0 : 1, st * sizeof(int));

                 if(sty == 9 || sty == 13) ++stt;
                 extract_sdt(trsav, sp, smth+st-1, kern+st-1, stt, sty, proj, &g1, &g2);

              }

          }

          printf("input next identifier\n");

          scanf("%d", &sty);


      }


   }


   else {

retry:

       printf("which stat type do you require, those available are:-\n");
       stat_type();

       scanf("%d", &sty);

       if(sty && sty <= STNM){

          st = 1;
          stt = 1;

          nmp = trsav->ptnum;

          fty = (int * )malloc_initl(sizeof(int));
          mem_er((fty == NULL) ? 0 : 1, sizeof(int));

          if(sty == 9 || sty == 13){

            sp = (float * )calloc(2*nmp, stor);  
            mem_er((sp == NULL) ? 0 : 1, 2*nmp * stor);

            *fty = 2;

          }

          else{

             sp = (float * )calloc(nmp, stor);  
             mem_er((sp == NULL) ? 0 : 1, nmp * stor);

             *fty = 1;

          } 


          ity = (int * )malloc_initl(sizeof(int));
          mem_er((ity == NULL) ? 0 : 1, sizeof(int));

          *ity = sty;

          smth = (float * )malloc_initl(sizeof(float));
          mem_er((smth == NULL) ? 0 : 1, sizeof(float));

          kern = (int * )malloc_initl(sizeof(int));
          mem_er((kern == NULL) ? 0 : 1, sizeof(int));

          if(sty == 9 || sty == 13) ++stt;
          extract_sdt(trsav, sp, smth, kern, stt, sty, proj, &g1, &g2);

        }

        else {

           printf("***ERROR***, no such stat identifier, try again\n");

           goto retry;

        }

      
   }

   cb = c1;

   if(st){

     if(!(g2.prgr))

#ifdef   NOUNDERSCORE

        statspl(xp, yp, sp, zp, &ca, &cb, &mv, &trsav->ptnum, &st, &stt, ity, fty, 
                cm->cmi,cm->cmxg, cm->cmyg, &cm->dcm, &xa1, &xa2, &ya1, &ya2, &ct,  
                grid, &npx, &npy, smth, kern, &gr->prgr, tdgrd);

#else

        statspl_(xp, yp, sp, zp, &ca, &cb, &mv, &trsav->ptnum, &st, &stt, ity, fty, 
                 cm->cmi,cm->cmxg, cm->cmyg, &cm->dcm, &xa1, &xa2, &ya1, &ya2, &ct,  
                 grid, &npx, &npy, smth, kern, &gr->prgr, tdgrd);

#endif

      else {

         xa1 = *(gr1->xgrid);
         xa2 = *(gr1->xgrid + gr1->ix - 1);
         ya1 = *(gr1->ygrid);
         ya2 = *(gr1->ygrid + gr1->iy - 1);

         printf("\nCurrent region is X:- %f -- %f            \r\n"
                "                  Y:- %f -- %f            \r\n"
                "Do you want a different region, 'y' or 'n'\n\n", xa1, xa2, ya1, ya2);

         scanf("\n");
         if(getchar() == 'y'){
            printf("What is the new region for X, xa1, xa2?\n\n");
            scanf("%f %f", &xa1, &xa2);
            printf("What is the new region for Y, ya1, ya2?\n\n");
            scanf("%f %f", &ya1, &ya2);

         }


#ifdef   NOUNDERSCORE

         statspl(xp, yp, sp, zp, &ca, &cb, &mv, &trsav->ptnum, &st, &stt, ity, fty, 
                 cm1->cmi,cm1->cmxg, cm1->cmyg, &cm1->dcm, &xa1, &xa2, &ya1, &ya2,   
                 &ct,grid, &npx, &npy, smth, kern, &gr->prgr, tdgrd);

#else

         statspl_(xp, yp, sp, zp, &ca, &cb, &mv, &trsav->ptnum, &st, &stt, ity, fty, 
                  cm1->cmi,cm1->cmxg, cm1->cmyg, &cm1->dcm, &xa1, &xa2, &ya1, &ya2,   
                  &ct,grid, &npx, &npy, smth, kern, &gr->prgr, tdgrd);

#endif

         if(c2) {

             ca = c1;
             cb = trsav->ptnum - c1;
             xa1 = *(gr2->xgrid);
             xa2 = *(gr2->xgrid + gr2->ix - 1);
             ya1 = *(gr2->ygrid);
             ya2 = *(gr2->ygrid + gr2->iy - 1);

             printf("\nCurrent region is X:- %f -- %f            \r\n"
                    "                  Y:- %f -- %f            \r\n"
                    "Do you want a different region, 'y' or 'n'\n\n", xa1, xa2, ya1, ya2);

             scanf("\n");
             if(getchar() == 'y'){
                printf("What is the new region for X, xa1, xa2?\n\n");
                scanf("%f %f", &xa1, &xa2);
                printf("What is the new region for Y, ya1, ya2?\n\n");
                scanf("%f %f", &ya1, &ya2);

             }


#ifdef   NOUNDERSCORE

             statspl(xp+c1, yp+c1, sp, zp, &ca, &cb, &mv, &trsav->ptnum, &st, &stt, ity, 
                     fty, cm2->cmi,cm2->cmxg, cm2->cmyg, &cm2->dcm, &xa1, &xa2,   
                     &ya1,&ya2, &ct, grid,&npx, &npy, smth, kern, &gr->prgr,
                     tdgrd);

#else

             statspl_(xp+c1, yp+c1, sp, zp, &ca, &cb, &mv, &trsav->ptnum, &st, &stt, ity, 
                      fty, cm2->cmi,cm2->cmxg, cm2->cmyg, &cm2->dcm, &xa1, &xa2,   
                      &ya1,&ya2, &ct, grid,&npx, &npy, smth, kern, &gr->prgr, 
                      tdgrd);

#endif



         }

      }

   }

   else {

       printf("*****ERROR*****, no statistics chosen for plotting.\n\n");

   }

   free(sp);
   free(ity);
   free(fty);
   free(smth);
   free(kern);
   free(xp);
   free(yp);
   free(zp);


   free(grid);
   free(tdgrd);

   if(g1.prty != g2.prty || g1.prgr != g2.prgr) proj_group(g1.prgr, g1.prty);

   return;

}

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                         SUPPORT ROUTINES

   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */


void stat_type()

{

    printf(" '1',  mean strength                         \r\n"
           " '2',  standard deviation of strength        \r\n"
           " '3',  mean speed                            \r\n"
           " '4',  standard deviation of speed           \r\n"
           " '5',  phenomena density                     \r\n"
           " '6',  genesis                               \r\n"
           " '7',  lysis                                 \r\n"
           " '8',  track density                         \r\n"
           " '9',  mean velocity                         \r\n"
           " '10', mean lifetime                         \r\n"
           " '11', growth/decay rate                     \r\n"
           " '12', feature anistropy                     \r\n"
           " '13', feature orientation                   \r\n"
           " '14', tendency                              \r\n"
           " '15', mean area (not currently avalilable). \r\n"
           " '16', additional or derived data            \r\n"
           " '17', additional or derived data            \n\n");

    return;

}


void extract_sdt(struct tot_stat *trsav, float *sp, float *sm, int *kty, int asp, int st, PROJ proj, GRDP *g1, GRDP *g2)

{

   int nstp;
   int i;
   int c1=0, c2=0;

   float *ss=NULL, *s7=NULL;
   float xcmp, ycmp;

   struct pt_stat *pst;

   nstp = (st == 9 || st == 13) ? (asp-2)*trsav->ptnum : (asp-1)*trsav->ptnum;

   ss = sp + nstp;

   switch(st){

      case 1:
         for(i=0; i < trsav->ptnum; i++) {
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = (pst->stat1).mean;
           }
           else {*(ss+c2) = (pst->stat1).mean; ++c2;}
         }
         *sm = trsav->sm[0];
         *kty = trsav->kern[0];
         break;
      case 2:
         for(i=0; i < trsav->ptnum; i++){
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = (pst->stat1).var;
           }
           else {*(ss+c2) = (pst->stat1).var; ++c2;}
         }
         *sm = trsav->sm[1];
         *kty = trsav->kern[1];
         break;
      case 3:
         for(i=0; i < trsav->ptnum; i++){
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = (pst->stat2).mean;
           }
           else{*(ss+c2) = (pst->stat2).mean; ++c2;}
         }
         *sm = trsav->sm[2];
         *kty = trsav->kern[2];
         break;
      case 4:
         for(i=0; i < trsav->ptnum; i++) {
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = (pst->stat2).var;
           }
           else{*(ss+c2) = (pst->stat2).var; ++c2;}
         }
         *sm = trsav->sm[3];
         *kty = trsav->kern[3];
         break;
      case 5:
         for(i=0; i < trsav->ptnum; i++) {
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = pst->stat3;
           }
           else {*(ss+c2) = pst->stat3; ++c2;}
         } 
         *sm = trsav->sm[4];
         *kty = trsav->kern[4];
         break;
      case 6:
         for(i=0; i < trsav->ptnum; i++) {
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = pst->stat4;
           }
           else {*(ss+c2) = pst->stat4; ++c2;}
         } 
         *sm = trsav->sm[5];
         *kty = trsav->kern[5];
         break;
      case 7:
         for(i=0; i < trsav->ptnum; i++) {
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = pst->stat5;
           }
           else {*(ss+c2) = pst->stat5; ++c2;}
         }
         *sm = trsav->sm[6];
         *kty = trsav->kern[6];
         break;
      case 8:
         for(i=0; i < trsav->ptnum; i++){
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = pst->stat6;
           }
           else {*(ss+c2) = pst->stat6; ++c2;}
         }
         *sm = trsav->sm[7];
         *kty = trsav->kern[7];
         break;
      case 9:
         s7 = ss + trsav->ptnum;
         for(i=0; i < trsav->ptnum; i++) {

             pst = trsav->ptst + i;

             xcmp = (pst->stat7).xcomp; 
             ycmp = (pst->stat7).ycomp;         

/* problem here if both hemispheres are plotted */

             if(g2->prgr) vect_proj(pst->xs, pst->ys, &xcmp, &ycmp, proj, g1, g2);

             if(pst->ptyp){
                ++c1;
                *(ss + trsav->ptnum - c1) = xcmp;
                *(s7 + trsav->ptnum - c1) = ycmp;
             }
             else {
                *(ss + c2) = xcmp;
                *(s7 + c2) = ycmp;
                ++c2;
             }
         }
         *sm = trsav->sm[8];
         *kty =  trsav->kern[8];
         break; 
      case 10:
         for(i=0; i < trsav->ptnum; i++) {
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = pst->stat8;
           }
           else {*(ss+c2) = pst->stat8; ++c2;}
         }
         *sm = trsav->sm[9];
         *kty = trsav->kern[9];
         break;
      case 11:
         for(i=0; i < trsav->ptnum; i++) {
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = pst->stat9;
           }
           else {*(ss+c2) = pst->stat9; ++c2;}
         }
         *sm = trsav->sm[10];
         *kty = trsav->kern[10];
         break;
      case 12:
         for(i=0; i < trsav->ptnum; i++) {
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = pst->stat10;
           }
           else {*(ss+c2) = pst->stat10; ++c2;}
         }
         *sm = trsav->sm[11];
         *kty = trsav->kern[11];
         break;
      case 13:
         s7 = ss + trsav->ptnum;
         for(i=0; i < trsav->ptnum; i++) {

             pst = trsav->ptst + i;

             xcmp = (pst->stat11).xcomp; 
             ycmp = (pst->stat11).ycomp;     

             if(g2->prgr) vect_proj(pst->xs, pst->ys, &xcmp, &ycmp, proj, g1, g2);

             if(pst->ptyp){
                ++c1;
                *(ss + trsav->ptnum - c1) = xcmp;
                *(s7 + trsav->ptnum - c1) = ycmp;
             }
             else {
                *(ss + c2) = xcmp;
                *(s7 + c2) = ycmp;
                ++c2;
             }
         }
         *sm = trsav->sm[12];
         *kty =  trsav->kern[12];
         break;
      case 14:
         for(i=0; i < trsav->ptnum; i++) {
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = pst->stat12;
           }
           else {*(ss+c2) = pst->stat12; ++c2;}
         }
         *sm = trsav->sm[13];
         *kty = trsav->kern[13];
         break;
      case 15:
         for(i=0; i < trsav->ptnum; i++) {
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = pst->stat13;
           }
           else {*(ss+c2) = pst->stat13; ++c2;}
         }
         *sm = trsav->sm[14];
         *kty = trsav->kern[14];
         break;
      case 16:
         for(i=0; i < trsav->ptnum; i++) {
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = pst->spare1;
           }
           else {*(ss+c2) = pst->spare1; ++c2;}
         }
         *sm = trsav->sm[15];
         *kty = trsav->kern[15];
         break;
      case 17:
         for(i=0; i < trsav->ptnum; i++) {
           pst = trsav->ptst + i;
           if(pst->ptyp){
              ++c1;
              *(ss + trsav->ptnum - c1) = pst->spare2;
           }
           else {*(ss+c2) = pst->spare2; ++c2;}
         }
         *sm = trsav->sm[16];
         *kty = trsav->kern[16];
         break;
      default:
         printf("***ERROR***, no such stat identifier for file %s\n", __FILE__);
         exit(1);
   }


   return;


}

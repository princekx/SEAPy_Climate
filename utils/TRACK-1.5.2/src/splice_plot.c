#include <Stdio.h>
#include <stdlib.h>
#include "proj.h"
#include "grid.h"
#include "mem_er.h"
#include "splice.h"

/* function to produce plotting arrays and plot a combined set of track data */

#ifdef   NOUNDERSCORE

void featpl(float *, float *, int *, int *, float *, float *, float *, int *, 
             int *, int *, float *, float *, int *, int *, float *, float *,
             float *, float *, int *, float *, float *, float *, float *,
             int *, int *, int *, int *, float *, float *, float *);

#else

void featpl_(float *, float *, int *, int *, float *, float *, float *, int *, 
             int *, int *, float *, float *, int *, int *, float *, float *,
             float *, float *, int *, float *, float *, float *, float *,
             int *, int *, int *, int *, float *, float *, float *);

#endif

extern int delb;

extern int aniso;

extern int tl, gof;

extern float w1, w2, dmax, phimax;
extern float xmn, ymn, xmx, ymx;

extern float period;

extern GRID *gr;
extern CNTRY *cm;

extern GRID *gr1, *gr2;
extern CNTRY *cm1, *cm2;

extern int nfld, nf;
extern int *nfwpos;

void splice_plot(struct tot_tr *all_tr, int tot_fet, int tr_count, int izm, int ntr, int trtyp)

{

     int i, j, count=0;
     int *tid=NULL;
     int tr_orig=tr_count;

     int k;
     int ipj;
     int ta;
     int c1=0, c2=1;
     int arr;
     int ix1=0, ix2=0, iy1=0, iy2=0;
     int pl_vec=0;
     int tr_mode='a';
     int totf, trid;
     int f1, f2;
     int rfr='n';
     int fnc1, fnc2;
     int iadd=0, iff=0;
     int igrwth=0;

     float pmax=0., pmin=0.;

     float tset=0.0;

     float *xf=NULL, *yf=NULL, *zf=NULL;
     float *xor=NULL, *yor=NULL, *zor=NULL;

     float xa, ya, za, xt, yt;

     struct fet_pt_tr *atr;
     struct tot_tr *altr;

     PROJ proj={NULL, NULL};
     
     GRDP g1={0, 0, 0.0, 0.0, 0.0, 0.0};
     GRDP g2={0, 0, 0.0, 0.0, 0.0, 0.0};
     GRDP g22={0, 0, 0.0, 0.0, 0.0, 0.0};

     if(tl == 'y'){

        printf("****INFORMATION****, grid has been translated.\n\n");

        printf("Do you want to translate the tracks by the same amount, 'y' or 'n'?\n\n");
        scanf("\n");

        if(getchar() == 'y'){

           if(gof > 0) tset = (*(gr->xgrid + gof) - *(gr->xgrid));
           else if(gof < 0) tset = -(*(gr->xgrid + gr->ix - 1) - *(gr->xgrid + gr->ix + gof - 1));


           for(i=0; i < tr_count; i++){

               altr = all_tr + i;

               for(j=0; j < altr->num; j ++){

                  atr = (altr->trpt) + j;

                  atr->xf += tset;

                  if(atr->xf > *(gr->xgrid + gr->ix - 1)) atr->xf -= period;
                  else if(atr->xf < 0.) atr->xf += period;

                   

               }

           }

        }

     }

     if(nf && trtyp != 'v'){

        printf("There are %d additional fields available, specify which one to use \r\n"
               "for intensity or input '0' for default.                            \n\n", nf);

        scanf("%d", &iadd);
        if(iadd < 0 || iadd > nf) {
           printf("****ERROR****, additional field %d identifier not known, using default.\n\n", iadd);
           iadd = 0;
        }
        if(iadd){
           iff = 0;
           for(i=0; i < iadd; i++){
              if(*(nfwpos + i)) iff += 3;
              else iff += 1;
           }
           --iff;
        }

     }

     else if(trtyp == 'v'){
        printf("Do you want to use tendency or growth rate as the intensity value, \r\n"
               "input '1' for tendency or '2' for growth rate.                     \n\n");
        scanf("%d", &igrwth);
     }

     printf("Do you want all tracks plotted, or individual tracks, 'a' or 'i'\n");
     scanf("\n");

     tr_mode = getchar();

     switch(tr_mode){
          case 'a':
            totf = tot_fet;
            break;
          case 'i':
            trid = tr_orig+1;
            while(trid > tr_orig){
            
              printf("Which track ID, refer to track output file ff_trs.***\n");
              scanf("%d", &trid);

              if(trid > tr_orig)
                 printf("****WARNING****, track identifier to large,    \r\n"
                        "                 number of tracks is %d,        \r\n"
                        "                 choose a different identifier \n\n", tr_orig);
            }

            totf = (all_tr + trid - 1)->num;
            tr_count = 1;
            break;
          default:
            printf("****ERROR****, incorrect specifier for plotting tracks,\r\n"
                   "               in file %s, defaulting to 'all'\n", __FILE__);
            tr_mode = 'a';
            totf = tot_fet;
            break;
            
     }

     xf = (float *)calloc(totf, sizeof(float));
     mem_er((xf == NULL) ? 0 : 1, totf * sizeof(float));

     yf = (float *)calloc(totf, sizeof(float));
     mem_er((yf == NULL) ? 0 : 1, totf * sizeof(float));

     zf = (float *)calloc(totf, sizeof(float));
     mem_er((zf == NULL) ? 0 : 1, totf * sizeof(float));

     tid = (int *)calloc(totf, sizeof(int));
     mem_er((tid == NULL) ? 0 : 1, totf * sizeof(int));

     if(aniso == 'y' && !gr->prgr){

        printf("do you want the orientation vectors plotted also, '1' for yes or '0' for no\n\n");

        scanf("%d", &pl_vec);

        if(pl_vec){

           xor = (float *)calloc(totf, sizeof(float));
           mem_er((xor == NULL) ? 0 : 1, totf * sizeof(float));

           yor = (float *)calloc(totf, sizeof(float));
           mem_er((yor == NULL) ? 0 : 1, totf * sizeof(float));

           zor = (float *)calloc(2 * totf, sizeof(float));
           mem_er((zor == NULL) ? 0 : 1, 2 * totf * sizeof(float));

        }

     }

     g1.prty = gr->prty;
     g1.prgr = gr->prgr;

     if(g1.prgr) {

        g1.alat = gr1->alat; g1.alng = gr1->alng;
	g1.sp1= gr1->sp1; g1.sp2 = gr1->sp2;

        ix1 = gr1->ox1u;
        ix2 = gr1->ox2u;
        iy1 = gr1->oy1u;
        iy2 = gr1->oy2u;

     }

     printf("do you wish to change the projection, 'y' or 'n'\n");

     proj_report(g1.prty, g1.prgr);

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
     
     

     for(i=0; i < tr_count; i++){

         altr = all_tr + i;

         if(tr_mode == 'i') altr += trid - 1;

         if(altr->num > 0) {pmin = pmax = altr->trpt->zf; break;}

     }

newtrack:

     printf("Do you wish to restrict the frame range to be plotted, 'y' or 'n'\n\n");

     scanf("\n");
     if((rfr=getchar()) == 'y'){

        fnc1 = 0;
        fnc2 = 0;

        for(i=0; i < tr_count; i++){

            altr = all_tr + i;

            if(altr->trpt != NULL){

              if(!fnc1) fnc1 = fnc2 = (altr->trpt + altr->num - 1)->fr_id;

              atr = altr->trpt;

              if((atr + (altr->num - 1))->fr_id > fnc2) fnc2 = (atr + (altr->num - 1))->fr_id;
              if(atr->fr_id < fnc1) fnc1 = atr->fr_id; 

            }


        }

        printf("The available frames are %d --> %d\n\n", fnc1, fnc2);

        printf("Input the frame range, f1 < f2\n\n");

        scanf("%d %d", &f1, &f2);


        if((f1 > f2) || (f1 < 1) || (f2 > fnc2) ){
           printf("****WARNING****, incorrect frame range input, reverting to full range\n\n");
        rfr = 'n';

        }



     }

     count = 0;
     c1=0, c2=1;
     
     for(i=0; i < tr_count; i++){

           altr = all_tr + i;

           if(tr_mode == 'i') altr += trid - 1;

           for(j=0; j < altr->num; j++){

               atr = (altr->trpt) + j;

               if(rfr == 'y'){

                  if(atr->fr_id < f1 || atr->fr_id > f2) continue;

               }

               ta = i+1;
               xa = atr->xf;
               ya = atr->yf;

               if(igrwth == 1) za = atr->tend;
               else if(igrwth == 2) za = atr->gwthr;
               else za = (iadd) ? atr->add_fld[iff] : atr->zf;

               if(ipj == 'y'){
 
                 switch(g1.prgr){
                    case 0:
                       switch(g2.prgr){
                          case 0:
                            proj_point(proj, &xa, &ya, &g1, &g2);
                            tid[count] = ta;
                            xf[count] = xa;
                            yf[count] = ya;
                            zf[count] = za;
                            if(pl_vec){
                              xor[count] = atr->or_vec[0];
                              yor[count] = atr->or_vec[1];
                            }
                            break;
                          case 1: case 2:
                            if(proj_point(proj, &xa, &ya, &g1, &g2)){

                              tid[c1] = ta;
                              xf[c1] = xa;
                              yf[c1] = ya;
                              zf[c1] = za;
                              if(pl_vec){
                                xor[c1] = atr->or_vec[0];
                                yor[c1] = atr->or_vec[1];
                                vect_proj(atr->xf, atr->yf, xor+c1, yor+c1, proj, &g1, &g2);
                              }

                              ++c1;
                            }
                            else if(gr2){
                              proj_point(proj, &xa, &ya, &g1, &g22);
                              tid[tot_fet-c2] = ta;
                              xf[tot_fet-c2] = xa;
                              yf[tot_fet-c2] = ya;
                              zf[tot_fet-c2] = za;
                              if(pl_vec){
                                xor[tot_fet-c2] = atr->or_vec[0];
                                yor[tot_fet-c2] = atr->or_vec[1];
                                vect_proj(atr->xf, atr->yf, xor+tot_fet-c2, yor+tot_fet-c2, proj, &g1, &g22);
                              }
                              ++c2;
                            }
                            break;

                       }
                       break;

                    case 1: case 2:

                       switch(g2.prgr){
                          case 0:

                            proj_point(proj, &xa, &ya, &g1, &g2);
                            tid[count] = ta;
                            xf[count] = xa;
                            yf[count] = ya;
                            zf[count] = za;
                            break;
                          case 1: case 2:
                            xt = xa;
                            yt = ya;
                            if(proj_point(proj, &xt, &yt, &g1, &g2)){
                              tid[c1] = ta;
                              xf[c1] = xt;
                              yf[c1] = yt;
                              zf[c1] = za;
                              ++c1;
                            }
                            else if(gr2){
                              proj_point(proj, &xa, &ya, &g1, &g22);
                              tid[tot_fet-c2] = ta;
                              xf[tot_fet-c2] = xa;
                              yf[tot_fet-c2] = ya;
                              zf[tot_fet-c2] = za;
                              ++c2;
                            }
                            break;

                       }
                       break;

                 }


               }

               else{

                 tid[count] = ta;
                  xf[count] = xa;
                  yf[count] = ya;
                  zf[count] = za;
                  if(pl_vec){
                    xor[count] = atr->or_vec[0];
                    yor[count] = atr->or_vec[1];
                  }

               }


               if(za > pmax) pmax = za;
               if(za < pmin) pmin = za;

               ++count;

           }

     }

     tot_fet = count;

     if(ipj == 'n') c1 = count;

     printf("the min-max range for the strengths is %f to %f\n", pmin, pmax);

    switch(g2.prgr){
       case 0:
          arr = 0;
          k = -1;

          if(g1.prgr){
             xmn = *(gr->xgrid + ix1 - 1);
             xmx = *(gr->xgrid + ix2 - 1);
             ymn = *(gr->ygrid + iy1 - 1);
             ymx = *(gr->ygrid + iy2 - 1);
          }

#ifdef  NOUNDERSCORE

          featpl(gr->xgrid, gr->ygrid, &gr->ix, &gr->iy, xf, yf, zf, tid, 
                  &tot_fet, cm->cmi, cm->cmxg, cm->cmyg, &cm->dcm , &k, &w1,
                  &w2, &dmax, &phimax, &izm, &xmn, &ymn, &xmx, &ymx, &delb,
                  &arr, &ntr, &pl_vec, xor, yor, zor);

#else

          featpl_(gr->xgrid, gr->ygrid, &gr->ix, &gr->iy, xf, yf, zf, tid, 
                  &tot_fet, cm->cmi, cm->cmxg, cm->cmyg, &cm->dcm , &k, &w1,
                  &w2, &dmax, &phimax, &izm, &xmn, &ymn, &xmx, &ymx, &delb,
                  &arr, &ntr, &pl_vec, xor, yor, zor);

#endif
          break;
       case 1: case 2:
          arr = 0;
          k = -1;
          xa = *(gr1->xgrid);
          xt = *(gr1->xgrid + gr1->ix - 1);
          ya = *(gr1->ygrid);
          yt = *(gr1->ygrid + gr1->iy - 1); 

          printf("\nCurrent region is X:- %f -- %f            \r\n"
                 "                  Y:- %f -- %f            \r\n"
                 "Do you want a different region, 'y' or 'n'\n\n", xa, xt, ya, yt);

          scanf("\n");
          if(getchar() == 'y'){
             printf("What is the new region for X, xa, xt?\n\n");
             scanf("%f %f", &xa, &xt);
             printf("What is the new region for Y, ya, yt?\n\n");
             scanf("%f %f", &ya, &yt);

          }

#ifdef  NOUNDERSCORE
             
          featpl(gr1->xgrid, gr1->ygrid, &gr1->ix, &gr1->iy, xf, yf, zf,  
                  tid, &c1, cm1->cmi, cm1->cmxg, cm1->cmyg, &cm1->dcm , &k, 
                  &w1,&w2, &dmax, &phimax, &izm, &xa, &ya, &xt, &yt, &delb, 
                  &arr, &ntr, &pl_vec, xor, yor, zor);

#else

          featpl_(gr1->xgrid, gr1->ygrid, &gr1->ix, &gr1->iy, xf, yf, zf,  
                  tid, &c1, cm1->cmi, cm1->cmxg, cm1->cmyg, &cm1->dcm , &k, 
                  &w1,&w2, &dmax, &phimax, &izm, &xa, &ya, &xt, &yt, &delb, 
                  &arr, &ntr, &pl_vec, xor, yor, zor);

#endif

          if(gr2){
             arr = 1;
             k = -1;
             xa = *(gr2->xgrid);
             xt = *(gr2->xgrid + gr2->ix - 1);
             ya = *(gr2->ygrid);
             yt = *(gr2->ygrid + gr2->iy - 1);  

             printf("\nCurrent region is X:- %f -- %f            \r\n"
                    "                  Y:- %f -- %f            \r\n"
                    "Do you want a different region, 'y' or 'n'\n\n", xa, xt, ya, yt);

             scanf("\n");
             if(getchar() == 'y'){
                printf("What is the new region for X, xa, xt?\n\n");
                scanf("%f %f", &xa, &xt);
                printf("What is the new region for Y, ya, yt?\n\n");
                scanf("%f %f", &ya, &yt);

             }

#ifdef  NOUNDERSCORE
            
             featpl(gr2->xgrid, gr2->ygrid, &gr2->ix, &gr2->iy, xf+c1, yf+c1,
                     zf+c1,tid+c1, &c2, cm2->cmi, cm2->cmxg, cm2->cmyg, 
                     &cm2->dcm, &k, &w1,&w2, &dmax, &phimax, &izm, &xa, &ya,
                     &xt, &yt, &delb, &arr, &ntr, &pl_vec, xor+c1, yor+c1, 
                     zor+2*c1);

#else

             featpl_(gr2->xgrid, gr2->ygrid, &gr2->ix, &gr2->iy, xf+c1, yf+c1,
                     zf+c1,tid+c1, &c2, cm2->cmi, cm2->cmxg, cm2->cmyg, 
                     &cm2->dcm, &k, &w1,&w2, &dmax, &phimax, &izm, &xa, &ya,
                     &xt, &yt, &delb, &arr, &ntr, &pl_vec, xor+c1, yor+c1, 
                     zor+2*c1);

#endif

          }

     }



     if(tr_mode == 'i'){

        printf("Do you want to plot another track, 'y' or 'n'\n");
        scanf("\n");
        if(getchar() == 'y'){

           trid = tr_orig + 1;
           while(trid > tr_orig){

             printf("Which track ID, refer to track output file ff_trs.***\n");
             scanf("%d", &trid);

             if(trid > tr_orig)
                printf("****WARNING****, track identifier to large,    \r\n"
                       "                 number of tracks is %d,        \r\n"
                       "                 choose a different identifier \n\n", tr_orig);

           }


           totf = (all_tr + trid - 1)->num;


           xf = (float *)realloc_n(xf, totf*sizeof(float));
           mem_er((xf == NULL) ? 0 : 1, totf*sizeof(float));

           yf = (float *)realloc_n(yf, totf*sizeof(float));
           mem_er((yf == NULL) ? 0 : 1, totf* sizeof(float));

           zf = (float *)realloc_n(zf, totf*sizeof(float));
           mem_er((zf == NULL) ? 0 : 1, totf* sizeof(float));

           tid = (int *)realloc_n(tid, totf*sizeof(int));
           mem_er((tid == NULL) ? 0 : 1, totf* sizeof(int));

           if(pl_vec){

              xor = (float *)realloc_n(xor, totf*sizeof(float));
              mem_er((xor == NULL) ? 0 : 1, totf* sizeof(float));

              yor = (float *)realloc_n(yor, totf*sizeof(float));
              mem_er((yor == NULL) ? 0 : 1, totf* sizeof(float));

              zor = (float *)realloc_n(zor, 2*totf*sizeof(float));
              mem_er((zor == NULL) ? 0 : 1, 2*totf* sizeof(float));

           }


           goto newtrack;

        }


     }



     if(g1.prty != g2.prty) {

        printf("***IMPORTANT***, you must now transform back to original grid\n\n");

        proj_group(g1.prgr, g1.prty);

     }


     free(xf);
     free(yf);
     free(zf);
     free(tid);

     if(pl_vec){

        free(xor);
        free(yor);
        free(zor);

     }


     return;

}

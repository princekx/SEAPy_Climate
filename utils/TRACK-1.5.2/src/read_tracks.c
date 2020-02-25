#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include "complex.h"
#include "splice.h"
#include "mem_er.h"

/* function to read in track data obtained by combining different sets
   of track data                                                        */

#define  FIELDWD  2000

int wc(char * );
void time2fr_id(struct tot_tr * , int );

float sum_wt=0.0, sum_per=0.0;
int iper_num=0;
int gppr=-1, ippr=-1;

extern int aniso;
extern int nfld, nf;
extern int *nfwpos;

struct tot_tr *read_tracks(FILE *tsf, int *tr_count, int *gpr, int *ipr, int fld, float *alat, float *alng)

{

    int i, j, k;
    int awt=0;
    int wron=0;
    int itime=0;
    int npos=0;

    char fldr[FIELDWD];
    char *fldt=NULL;

    struct tot_tr *all_tr=NULL, *altr=NULL;
    struct fet_pt_tr *atr=NULL;


    iper_num = 0;

    fgets(fldr, FIELDWD, tsf);

    if(wc(fldr) == 1){
       sscanf(fldr, "%d", &aniso);
       if(aniso == 1) aniso = 'y';
       else aniso = 'n';

       fgets(fldr, FIELDWD, tsf);

       if(strstr(fldr, "WT_INFO")) {
          sscanf(fldr, "%*s %d %f", &awt, &sum_wt);
	  fgets(fldr, FIELDWD, tsf);
          sscanf(fldr, "%d %d", gpr, ipr);

       }
       else if(strstr(fldr, "PER_INFO")) {
          sscanf(fldr, "%*s %d %f", &iper_num, &sum_per);
	  fgets(fldr, FIELDWD, tsf);
          sscanf(fldr, "%d %d", gpr, ipr);
       }
       else if(wc(fldr) == 1) {
          sscanf(fldr, "%d", &awt);
	  fgets(fldr, FIELDWD, tsf);
          sscanf(fldr, "%d %d", gpr, ipr);
       }
       else sscanf(fldr, "%d %d", gpr, ipr);

    }
 
    else {sscanf(fldr, "%d %d", gpr, ipr); aniso = 'n';}
    
    gppr = *gpr;
    ippr = *ipr;

    if(*gpr) {
      fgets(fldr, FIELDWD, tsf);
      sscanf(fldr, "%f %f", alat, alng);
    }
    
    fgets(fldr, FIELDWD, tsf);
    sscanf(fldr, "%*s %d %*s %d %d", tr_count, &nf, &nfld);

    if(!nf){
    
       printf("***INFORMATION***, there are no additional fields associated with this track data.\n\n");
	
    }
    else {
       if(nfwpos) free(nfwpos);
       nfwpos = (int *)calloc(nf, sizeof(int));
       mem_er((nfwpos == NULL) ? 0 : 1, nf * sizeof(int));
       fldt = strchr(fldr, '&') + 1;
       for(i=0; i < nf; i++) {sscanf(fldt, "%1d", nfwpos + i); ++fldt; if(*(nfwpos + i)) ++npos;}
       printf("****INFORMATION****, there are %d additional fields associated with this track data.\r\n"
              "                     of which %d have positional information.                       \n\n", 
              nf, npos);
    }

    printf("****INFORMATION****, number of tracks is %d\n", *tr_count);

    all_tr = (struct tot_tr * )calloc(*tr_count, sizeof(struct tot_tr));
    mem_er((all_tr == NULL) ? 0 : 1, *tr_count * sizeof(struct tot_tr));

    for(i=0; i < *tr_count; i++){

        altr = all_tr + i;

        altr->awt = awt;

        fgets(fldr, FIELDWD, tsf);

        sscanf(fldr, "%*s %d %*s %10ld", &altr->trid, &altr->time);

        if(altr->time) itime = 1;

        fgets(fldr, FIELDWD, tsf);
        sscanf(fldr, "%*s %d", &altr->num);

        altr->trpt = (struct fet_pt_tr * )calloc(altr->num, sizeof(struct fet_pt_tr));
        mem_er((altr->trpt == NULL) ? 0 : 1, altr->num * sizeof(struct fet_pt_tr));

        if(fld == 's'){

           if(aniso == 'y' && awt){

              for(j=0; j < altr->num; j++){

                  atr = (altr->trpt) + j;
		  
	          fgets(fldr, FIELDWD, tsf);
                  atr->nfm = 0;
		  
		  if(nfld){
		     atr->add_fld = (float *)calloc(nfld, sizeof(float));
		     mem_er((atr->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));

                     if(altr->time)

                        sscanf(fldr, "%10ld %f %f %e", &atr->time, &atr->xf, &atr->yf, &atr->zf);

                     else

                        sscanf(fldr, "%d %f %f %e", &atr->fr_id, &atr->xf, &atr->yf, &atr->zf);
		     
		     fldt = strchr(fldr, '&') + 1;
		     for(k=0; k < nfld; k++) {
		         sscanf(fldt, "%e", atr->add_fld + k);
		         fldt = strchr(fldt, '&') + 1;		  
		     }
		     
		     sscanf(fldt, "%f %f %f %e %f %d", &(atr->sh_an), &(atr->or_vec[0]), &(atr->or_vec[1]), &(atr->area), &atr->wght, &atr->nfm);

		  }
		  
		  else {
		  
                     if(altr->time)

                        sscanf(fldr, "%10ld %f %f %e %f %f %f %e %f %d", &atr->time, &atr->xf, &atr->yf, &atr->zf, &(atr->sh_an), &(atr->or_vec[0]), &(atr->or_vec[1]), &(atr->area), &atr->wght, &atr->nfm);

                     else

                        sscanf(fldr, "%d %f %f %e %f %f %f %e %f %d", &atr->fr_id, &atr->xf, &atr->yf, &atr->zf, &(atr->sh_an), &(atr->or_vec[0]), &(atr->or_vec[1]), &(atr->area), &atr->wght, &atr->nfm);
					  
		  }

              }

           }

           else if(aniso == 'y' && !awt){

              for(j=0; j < altr->num; j++){

                  atr = (altr->trpt) + j;
		  
                  fgets(fldr, FIELDWD, tsf);
                  atr->nfm = 0;
		  
		  if(nfld){
		     atr->add_fld = (float *)calloc(nfld, sizeof(float));
		     mem_er((atr->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));

                     if(altr->time) 

                        sscanf(fldr, "%10ld %f %f %e", &atr->time, &atr->xf, &atr->yf, &atr->zf);

                     else

                        sscanf(fldr, "%d %f %f %e", &atr->fr_id, &atr->xf, &atr->yf, &atr->zf);

		     fldt = strchr(fldr, '&') + 1;
		     for(k=0; k < nfld; k++) {
		         sscanf(fldt, "%e", atr->add_fld + k);
		         fldt = strchr(fldt, '&') + 1;		  
		     }
		  
                     sscanf(fldt, "%f %f %f %e %d", &(atr->sh_an), &(atr->or_vec[0]), &(atr->or_vec[1]), &(atr->area), &atr->nfm);	  

		  }
		  
		  else {
		  
		     if(altr->time) 

                        sscanf(fldr, "%10ld %f %f %e %f %f %f %e %d", &atr->time, &atr->xf, &atr->yf, &atr->zf, &(atr->sh_an), &(atr->or_vec[0]), &(atr->or_vec[1]), &(atr->area), &atr->nfm);

                     else

                        sscanf(fldr, "%d %f %f %e %f %f %f %e %d", &atr->fr_id, &atr->xf, &atr->yf, &atr->zf, &(atr->sh_an), &(atr->or_vec[0]), &(atr->or_vec[1]), &(atr->area), &atr->nfm);
			
		  }

              }

           }

           else if(aniso == 'n' && awt){

              for(j=0; j < altr->num; j++){

                  atr = (altr->trpt) + j;
		  
		  fgets(fldr, FIELDWD, tsf);
                  atr->nfm = 0;
		  
		  if(nfld){
		     atr->add_fld = (float *)calloc(nfld, sizeof(float));
		     mem_er((atr->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));

                     if(altr->time)

                        sscanf(fldr, "%10ld %f %f %e", &atr->time, &atr->xf, &atr->yf, &atr->zf);

                     else

                        sscanf(fldr, "%d %f %f %e", &atr->fr_id, &atr->xf, &atr->yf, &atr->zf);

		     fldt = strchr(fldr, '&') + 1;
		     for(k=0; k < nfld; k++) {
		         sscanf(fldt, "%e", atr->add_fld + k);
		         fldt = strchr(fldt, '&') + 1;		  
		     }
		  
		     sscanf(fldt, "%f %d", &atr->wght, &atr->nfm);
		  
		  }
		  
		  else {
		  
		     if(altr->time)

                        sscanf(fldr, "%10ld %f %f %e %f %d", &atr->time, &atr->xf, &atr->yf, &atr->zf, &atr->wght, &atr->nfm);

                     else

                        sscanf(fldr, "%d %f %f %e %f %d", &atr->fr_id, &atr->xf, &atr->yf, &atr->zf, &atr->wght, &atr->nfm);
		  
		  }

              }


           }

           else {

              for(j=0; j < altr->num; j++){

                  atr = (altr->trpt) + j;
		  
		  fgets(fldr, FIELDWD, tsf);
                  atr->nfm = 0;
		  
		  if(nfld){
		     atr->add_fld = (float *)calloc(nfld, sizeof(float));
		     mem_er((atr->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));

                     if(altr->time)

                        sscanf(fldr, "%10ld %f %f %e", &atr->time, &atr->xf, &atr->yf, &atr->zf);

                     else

                        sscanf(fldr, "%d %f %f %e", &atr->fr_id, &atr->xf, &atr->yf, &atr->zf);
		     
		     fldt = strchr(fldr, '&') + 1;
		     for(k=0; k < nfld; k++) {
		         sscanf(fldt, "%e", atr->add_fld + k);
		         fldt = strchr(fldt, '&') + 1;		  
		     }
		  
                     sscanf(fldt, "%d", &atr->nfm);
		  
		  }
		  
		  else {
		  
		     if(altr->time)

                        sscanf(fldr, "%10ld %f %f %e %d", &atr->time, &atr->xf, &atr->yf, &atr->zf, &atr->nfm);

                     else

                        sscanf(fldr, "%d %f %f %e %d", &atr->fr_id, &atr->xf, &atr->yf, &atr->zf, &atr->nfm);
		  
		  }	  

              }

           }

        }

        else if(fld == 'v'){


           if(!wron) {
              printf("****WARNING****, older track, phase speed files may not be read correctely,\r\n"
                     "                 there should be 1 integer and 7 float fields for          \r\n"
                     "                 correct read.                                             \n\n");

              wron = 1;

           }

           if(!awt){

              for(j=0; j < altr->num; j++){

                  atr = (altr->trpt) + j;

                  fgets(fldr, FIELDWD, tsf);

                  if(altr->time)

                     sscanf(fldr, "%10ld %f %f %e %e %e %e %e", &atr->time, &atr->xf, &atr->yf, &atr->zf, &atr->gwthr, &atr->tend, &atr->vec[0], &atr->vec[1]);


                  else

                     sscanf(fldr, "%d %f %f %e %e %e %e %e", &atr->fr_id, &atr->xf, &atr->yf, &atr->zf, &atr->gwthr, &atr->tend, &atr->vec[0], &atr->vec[1]);

              }

           }

           else {

              for(j=0; j < altr->num; j++){

                  atr = (altr->trpt) + j;

                  fgets(fldr, FIELDWD, tsf);

                  if(altr->time)

                     sscanf(fldr, "%10ld %f %f %e %e %e %e %e %f", &atr->time, &atr->xf, &atr->yf, &atr->zf, &atr->gwthr, &atr->tend, &atr->vec[0], &atr->vec[1], &atr->wght);

                  else

                     sscanf(fldr, "%d %f %f %e %e %e %e %e %f", &atr->fr_id, &atr->xf, &atr->yf, &atr->zf, &atr->gwthr, &atr->tend, &atr->vec[0], &atr->vec[1], &atr->wght);

              }

           }

        }

    }

    if(itime){

/* convert time information to frame ID's */

       printf("Do you want to convert real time to frame Id's, 'y' or 'n'\n\n");
       scanf("\n");
       if(getchar() == 'y'){

          printf("****INFORMATION****, converting time info to frame ID's\n\n");

          time2fr_id(all_tr, *tr_count);

       }


    }

    return all_tr;

}


/* int wc(char *st)

{

    int nw=0, state=0;
    int i=0, c;

    while((c=st[i++]) != '\0'){

       if(c == ' ' || c == '\n' || c == '\t')

          state = 0;

       else if (!state){

          state = 1;
          ++nw;
       }


    }

    return nw;

}      */

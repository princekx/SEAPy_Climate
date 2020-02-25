#include <Stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include "grid.h"
#include "splice.h"
#include "file_handle.h"

/* function to write to file the track data obtained by combining
   different sets of track data                                        */

extern GRID *gr;

extern float sum_wt, sum_per;

extern int aniso;
extern int iper_num;
extern int nfld, nf;
extern int *nfwpos;
extern int gppr, ippr;

void meantrd(FILE *tsf, struct tot_tr *all_tr, int tr_count, int fld)

{

    int i, j, k;
    int trc=0;
    int awt=0;

    off_t pl=0;

    struct tot_tr *altr;
    struct fet_pt_tr *atr;

    if(aniso == 'y') fprintf(tsf, "%d\n", 1);
    else fprintf(tsf, "%d\n", 0);

    if(iper_num){
       fprintf(tsf, "PER_INFO %1d %12.5f\n", 1, sum_per);
    }

    else {

       for(i=0; i < tr_count; i++) {

           if((all_tr + i)->awt) {
             awt = (all_tr + i)->awt;
             fprintf(tsf, "WT_INFO %1d %12.5f\n", awt, sum_wt);
             break;

           }
       }
    }

    if(gppr < 0){
       fprintf(tsf, "%d %d\n", gr->prgr, gr->prty);
       if(gr->prgr) fprintf(tsf, "%f %f\n", gr->alat, gr->alng);
    }
    else {
       fprintf(tsf, "%d %d\n", gppr, ippr);
       if(gppr) {
           if(gppr != gr->prgr || ippr != gr->prty){
	      printf("****ERROR****, projections for tracks and grid mismatch.\n\n");
	      exit(1);
	   }
           fprintf(tsf, "%f %f\n", gr->alat, gr->alng);
	    
       } 
       
    }

    pl = ftello(tsf);

    fprintf(tsf, "TRACK_NUM  %8d ADD_FLD  %3d %3d &", tr_count, nf, nfld);
    if(nfwpos){
       for(i=0; i < nf; i++)fprintf(tsf, "%1d", *(nfwpos + i));
    }
    fprintf(tsf, "\n");

    for(i=0; i < tr_count; i++){

        altr = all_tr + i;

        awt = altr->awt;

        if(altr->num > 0){

           ++trc;

           if(altr->time) fprintf(tsf, "TRACK_ID  %d START_TIME %10ld\n", altr->trid, altr->time);
           else fprintf(tsf, "TRACK_ID  %d\n", altr->trid);

           fprintf(tsf, "POINT_NUM  %d\n", altr->num);

           if(fld == 's'){

              if(aniso == 'y' && awt){

                 for(j=0; j < altr->num; j++){

                     atr = (altr->trpt) + j;

                     if(altr->time)

                        fprintf(tsf, "%10ld %f %f %e ", atr->time, atr->xf, atr->yf, atr->zf);

                     else

                        fprintf(tsf, "%d %f %f %e ", atr->fr_id, atr->xf, atr->yf, atr->zf);

                     if(nfld){
		        fprintf(tsf, "& ");		     
		        for(k=0; k < nfld; k++) fprintf(tsf, "%e & ", *(atr->add_fld + k));
		     }
		     
		     fprintf(tsf, "%f %f %f %e %f ", atr->sh_an, atr->or_vec[0], atr->or_vec[1], atr->area, atr->wght);

                     if(atr->nfm)
                        fprintf(tsf, "%d\n", atr->nfm);
                     else
                        fprintf(tsf, "\n");

                 }

              }

              else if(aniso == 'y' && !awt){

                 for(j=0; j < altr->num; j++){

                     atr = (altr->trpt) + j;

                     if(altr->time)

                        fprintf(tsf, "%10ld %f %f %e ", atr->time, atr->xf, atr->yf, atr->zf);

                     else 

                        fprintf(tsf, "%d %f %f %e ", atr->fr_id, atr->xf, atr->yf, atr->zf);

                     if(nfld){
		        fprintf(tsf, "& ");		     
		        for(k=0; k < nfld; k++) fprintf(tsf, "%e & ", *(atr->add_fld + k));
		     }
		     			
		     fprintf(tsf, "%f %f %f %e ", atr->sh_an, atr->or_vec[0], atr->or_vec[1], atr->area);

                     if(atr->nfm)
                        fprintf(tsf, "%d\n", atr->nfm);
                     else
                        fprintf(tsf, "\n");

                 }

              }

              else if(aniso == 'n' && awt){

                 for(j=0; j < altr->num; j++){

                     atr = (altr->trpt) + j;

                     if(altr->time)

                        fprintf(tsf, "%10ld %f %f %e ", atr->time, atr->xf, atr->yf, atr->zf);

                     else 

                        fprintf(tsf, "%d %f %f %e ", atr->fr_id, atr->xf, atr->yf, atr->zf);

                     if(nfld){
		        fprintf(tsf, "& ");		     
		        for(k=0; k < nfld; k++) fprintf(tsf, "%e & ", *(atr->add_fld + k));
		     }
		     			
		     fprintf(tsf, "%f ", atr->wght);

                     if(atr->nfm)
                        fprintf(tsf, "%d\n", atr->nfm);
                     else
                        fprintf(tsf, "\n");

                 }

              }

              else {

                 for(j=0; j < altr->num; j++){

                     atr = (altr->trpt) + j;

                     if(altr->time)

                        fprintf(tsf, "%10ld %f %f %e ", atr->time, atr->xf, atr->yf, atr->zf);

			
                     else 

                        fprintf(tsf, "%d %f %f %e ", atr->fr_id, atr->xf, atr->yf, atr->zf);

                     if(nfld){
		        fprintf(tsf, "& ");		     
		        for(k=0; k < nfld; k++) fprintf(tsf, "%e & ", *(atr->add_fld + k));
		     }

                     if(atr->nfm)
                        fprintf(tsf, "%d\n",atr->nfm);
                     else
                        fprintf(tsf, "\n");

                 }
  
              }

           }

           else if(fld == 'v'){

              if(!awt){

                 for(j=0; j < altr->num; j++){

                     atr = (altr->trpt) + j;

                     if(altr->time)

                        fprintf(tsf, "%10ld %f %f %e %e %e %e %e\n", atr->time, atr->xf, atr->yf, atr->zf, atr->gwthr, atr->tend, atr->vec[0], atr->vec[1]);

                     else

                        fprintf(tsf, "%d %f %f %e %e %e %e %e\n", atr->fr_id, atr->xf, atr->yf, atr->zf, atr->gwthr, atr->tend, atr->vec[0], atr->vec[1]);

                 }

              }

              else {

                 for(j=0; j < altr->num; j++){

                     atr = (altr->trpt) + j;

                     if(altr->time)

                        fprintf(tsf, "%10ld %f %f %e %e %e %e %e %f\n", atr->time, atr->xf, atr->yf, atr->zf, atr->gwthr, atr->tend, atr->vec[0], atr->vec[1], atr->wght);

                     else

                        fprintf(tsf, "%d %f %f %e %e %e %e %e %f\n", atr->fr_id, atr->xf, atr->yf, atr->zf, atr->gwthr, atr->tend, atr->vec[0], atr->vec[1], atr->wght);

                 }

              }


           }

        }

    } 

    if(trc != tr_count){

       fseeko(tsf, pl, FSTART);
       fprintf(tsf, "TRACK_NUM  %8d", trc);

    }


    return;

}

#include <Stdio.h>
#include <Math.h>
#include "splice.h"
#include "m_values.h"

/* function to convert track lat-long's to Cartesians */

void convert_track(struct tot_tr *all_tr, int trn, int ifd, int ifdp)

{

   int i, j;
   struct tot_tr *altr=NULL;
   struct fet_pt_tr *at=NULL;

   double sn1, cn1, sn2, cn2;
   double thet, phi;

   for(i=0; i < trn; i++){

      altr = all_tr + i;

      for(j=0; j < altr->num; j++){

          at = altr->trpt + j;
	  
	  if(ifd){
	    phi = *(at->add_fld + ifdp) * FP_PI;
            thet = FP_PI2 - *(at->add_fld + ifdp + 1) * FP_PI;
	  }
	  else{
             phi = at->xf * FP_PI; 
             thet = FP_PI2 - (at->yf) * FP_PI;
	  }
	  
          if(thet < 0.) thet = 0.;

          sincos(phi, &sn1, &cn1);
          sincos(thet, &sn2, &cn2);

          at->pp[0] = sn2 * cn1;
          at->pp[1] = sn2 * sn1;
          at->pp[2] = cn2;

      }


   }

   return;

}

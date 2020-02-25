#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include <string.h>
#include "complex.h"
#include "m_values.h"
#include "mem_er.h"

/* Discrete fast cosine transform based on Jain, Fundamentals of Digital
   Image Processing                                                      */
   
#ifdef  FFTLIB

#ifdef  NOUNDERSCORE

void gpfa(double * , double * , double * , int * , int * , int * , int * , int *, int * );

#else

void gpfa_(double * , double * , double * , int * , int * , int * , int * , int *, int * );

#endif

#endif

int powi(int , int );
void dfft(complex * , int , int , int );

complex *yct=NULL;
complex *wct=NULL;

int *nj=NULL;

double *dtrig=NULL;

extern double *dfreal, *dfimag;

void dfct(complex *xt, int nn, int nexp, int inv, int itrig, int ifftyp)
{
  
   int i=0;
   int nxp=0;
   
   int jump=0, inc=1, lot=1;
   
   double si=0.0, ci=0.0;
   double sn2=0.0, sn22=0.0, psn2=0.0;
   
   complex w;
   
   if(inv != 1){

      if(inv != -1){

         printf("****ERROR****, transform direction must be either foreward (1)\r\n"
                "               or backward (-1).                              \n\n");
         exit(1);

      }

   }
   
   if(itrig > 0 && !ifftyp){
   
/* test settings */
      nxp = (int)(log((float)nn)/log(2.0));
      if(nxp != nexp){
         printf("****ERROR****, incompatable power of 2 for dct.\n\n");
	 exit(1);
      }
   
   }
   else if (itrig < 0){free(wct); free(yct); return;}
   
   psn2 = FPI / ((float)(2 * nn));
   sn22 = sqrt((float)(2 * nn));
   sn2 = 2.0 / sqrt((float)(2 * nn));

   if(!yct){
     yct = (complex *)calloc(nn, sizeof(complex));
     mem_er((yct == NULL) ? 0 : 1, nn * sizeof(complex));
   }
   
   if(!wct){
     wct = (complex *)calloc(nn, sizeof(complex));
     mem_er((wct == NULL) ? 0 : 1, nn * sizeof(complex));   
   }
   
   memset(yct, 0, nn*sizeof(complex));
   
   if(itrig > 0){
      if(inv == 1){
         for(i=0; i < nn; i++){
             sincos((float) i * psn2, &si, &ci);
	     comp((ci * sn2), (-si * sn2), (wct + i));
         }
      }
      else {
         for(i=0; i < nn; i++){
             sincos((float) i * psn2, &si, &ci);
	     comp((ci * sn22), (si * sn22), (wct + i));
         }      
      }
      wct->real /= SQ_R2;
      wct->imag /= SQ_R2;
   }
     
   if(inv == 1){  
   
      if(!ifftyp){
         for(i=0; i < nn / 2; i++){
             ccopy((xt + 2 * i), (yct + i));
             ccopy((xt + 2 * i + 1), (yct + nn - i - 1));
         }
     
         dfft(yct, nexp, inv, itrig);

         for(i=0; i < nn; i++) {
            cmul(*(wct + i), *(yct + i), &w);
            ccopy(&w, xt + i);
         }

      } 
      
      else{
      
         for(i=0; i < nn / 2; i++){
	     *(dfreal + i) = (xt + 2 * i)->real;
	     *(dfimag + i) = (xt + 2 * i)->imag;
             *(dfreal + nn - i - 1) = (xt + 2 * i + 1)->real;
	     *(dfimag + nn - i - 1) = (xt + 2 * i + 1)->imag;
         }
      
#ifdef  FFTLIB

#ifdef  NOUNDERSCORE
         gpfa(dfreal, dfimag, dtrig, &inc, &jump, &nn, &lot, &inv, nj);
#else
         gpfa_(dfreal, dfimag, dtrig, &inc, &jump, &nn, &lot, &inv, nj);
#endif

#endif

          for(i=0; i < nn; i++){
	      (yct + i)->real = *(dfreal + i);
	      (yct + i)->imag = -(*(dfimag + i));
              cmul(*(wct + i), *(yct + i), &w);
              ccopy(&w, xt + i);
	  }
    
      }

   }
   else{
   
      if(!ifftyp){
         for(i=0; i < nn; i++){
             cmul(*(xt + i), *(wct + i), &w);
	     comp(w.real, w.imag, (yct + i));
         } 
      
         dfft(yct, nexp, inv, itrig); 
	 
         for(i=0; i < nn / 2; i++){
             ccopy((yct + i), (xt + 2 * i));
	     ccopy((yct + nn - i - 1), (xt + 2 * i + 1));
         }
	 
      }
      
      else {
         for(i=0; i < nn; i++){
             cmul(*(xt + i), *(wct + i), &w);
	     *(dfreal + i) = w.real;
	     *(dfimag + i) = w.imag;
         }

	 
#ifdef  FFTLIB

#ifdef  NOUNDERSCORE
         gpfa(dfreal, dfimag, dtrig, &inc, &jump, &nn, &lot, &inv, nj);
#else
         gpfa_(dfreal, dfimag, dtrig, &inc, &jump, &nn, &lot, &inv, nj);
#endif

#endif

         xt->real = *dfreal / (float) nn;
         xt->imag = *dfimag / (float) nn;
         for(i=1; i <  nn / 2; i++){
	     (xt + 2 * i - 1)->real = *(dfreal + i) / (float) nn;
	     (xt + 2 * i - 1)->imag = *(dfimag + i) / (float) nn; 
	     (xt + 2 * i)->real = *(dfreal + nn - i) / (float) nn;
	     (xt + 2 * i)->imag = *(dfimag + nn - i) / (float) nn;
         }
	 (xt + nn - 1)->real = *(dfreal + (nn / 2)) / (float) nn;
	 (xt + nn - 1)->imag = *(dfimag + (nn / 2)) / (float) nn;
  
      }
      
   }
     

   return;
    
}

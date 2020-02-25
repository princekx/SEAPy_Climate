#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "splice.h"
#include "mem_er.h" 

/* program to do identification for multiple fields */

typedef struct _reg{
    float rlng1;
    float rlng2;
    float rlat1;
    float rlat2;
    int ig;
} REG;                      /* region structure */

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
int region(struct fet_pt_tr * , REG * , int, int );

extern int nfld, nff;
extern int *nfwpos;

int noheader=0;

int main(void)
{
    int i=0, j=0, k=0;
    int trnum=0;
    int nft=0;
    int nct=0, nmct=0;
    int nmat=0;
    int nth=0;
    int gpr=0, ipr=0;
    int ireg=0;
    int ninreg=0;

    int *ifld=NULL;
    int *ifncnt=NULL;
    int *imnmx=NULL;

    float alat=0.0, alng=0.0;

    float *fthr=NULL;

    float fadd=0.0;

    FILE *fin=NULL, *fout=NULL;
    FILE *fcmp=NULL;

    char filin[MAXCHR];
    char filout[MAXCHR];
    char buff[MAXCHR];

    struct tot_tr *tracks=NULL, *atr=NULL;
    struct fet_pt_tr *pt=NULL;
    
    REG reg={0.0, 0.0, 0.0, 0.0, 0};
 
/* read region data is present */   
    
    fcmp = fopen("region.dat", "r");
    if(fcmp){
       printf("****WARNING****, file exists with region data for sub-setting tracks,   \r\n"
              "                 this information may be used to restrict the analysis  \r\n"
              "                 to the specified region.                               \n\n");
       fgets(buff, MAXCHR, fcmp);
       sscanf(buff, "%f %f", &(reg.rlng1), &(reg.rlng2));
       fgets(buff, MAXCHR, fcmp);
       sscanf(buff, "%f %f", &(reg.rlat1), &(reg.rlat2));
       fgets(buff, MAXCHR, fcmp);
       
       if(reg.rlng2 < reg.rlng1) {
         printf("****WARNING****, stradeling the Greenwich meridian\n\n");
         reg.ig = 1;
       }
       ireg = 1;

       fclose(fcmp); fcmp = NULL;
       
    }

    printf("What track file is required?\n\n");
    scanf("%s", filin);

    strcpy(filout, filin);
    strcat(filout, ".mident");
    
    fin = fopen(filin, "r");
    if(!fin){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", filin);
       exit(1);
    }

    tracks = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

    fclose(fin);

    printf("What is the number of contiguous points to satisfy the thresholds?\n\n");
    scanf("%d", &nth);
    
    printf("How many fields are required to be tested?\n\n");
    scanf("%d", &nft);

    ifld = (int *)calloc(nft, sizeof(int));
    mem_er((ifld == NULL) ? 0 : 1, nft * sizeof(int));

    ifncnt = (int *)calloc(nft, sizeof(int));
    mem_er((ifncnt == NULL) ? 0 : 1, nft*sizeof(int));

    imnmx = (int *) calloc(nft, sizeof(int));
    mem_er((imnmx == NULL) ? 0 : 1, nft*sizeof(int));

    fthr = (float *) calloc(nft, sizeof(float));
    mem_er((fthr == NULL) ? 0 : 1, nft*sizeof(float));

    for(i=0; i < nft; i++){
       printf("What is the %d field required? Default field is '0'\n", i+1);
       scanf("%d", ifld + i);
       printf("\n"); 

       if(*(ifld + i) > nff){
          printf("****ERROR****, chosen added field Id to large.\n\n");
          exit(1);
       } 

       printf("Is a min or max test required, '0' for min or '1' for max.\n\n");
       scanf("%d", imnmx + i);
       if(*(imnmx + i) < 0 || *(imnmx + i) > 1){
          printf("****ERROR****, incorrect option chosen.\n\n");
          exit(1);
       }

       printf("What threshold is required for this field?\n");
       scanf("%f", fthr + i);

       if(*(ifld + i)){
          *(ifncnt + i) = 0;
          for(j=0; j < *(ifld + i); j++){
             if(*(nfwpos + j)) *(ifncnt + i) += 3;
             else *(ifncnt + i) += 1;
          }
          --(*(ifncnt + i));
       }

    }
    

    for(i=0; i < trnum; i++){
       atr = tracks + i;

       nct = 0;
       nmct = 0;

       for(j=0; j< atr->num; j++){
          pt = atr->trpt + j;
	  
	  nmat = 0;
	  
	  ninreg=0;
	  
          for(k=0; k < nft; k++){
              if(! *(ifld + k)){
	         if(ireg){
		    if(region(pt, &reg, 0, 0)) ++ninreg;
		 }
	      
                 if(! *(imnmx + k)){
                    if(pt->zf <= *(fthr + k)) ++nmat;
                 } 
                 else {
                    if(pt->zf >= *(fthr + k)) ++nmat;
                 }
              }
              else {
	      
	         if(ireg){
	           if(*(nfwpos + *(ifld + k) - 1)){
		      if(region(pt, &reg, *(ifld + k), *(ifncnt + k))) ++ninreg;
		   }
		   else {
		     if(region(pt, &reg, 0, 0)) ++ninreg;
		   }
		 
		 }
		 
                 fadd = pt->add_fld[*(ifncnt + k)]; 	

                 if(fadd > ADD_CHECK) continue;
                 if(! *(imnmx + k)){
                    if(fadd <= *(fthr + k)) ++nmat;
                 }
                 else {
 
                    if(fadd >= *(fthr + k)) {		    
		       ++nmat;
		    }
                 }

              }      
          }
	  
	  if(ireg){
	     if(nmat == nft && ninreg == nft) ++nct;
	     else nct = 0;	     	     
	  }
	  else{
             if(nmat == nft) ++nct;
             else nct = 0;	  
	  }
          
	  if(nct > nmct) nmct = nct;

       }      

       if(nmct < nth){ 
         for(j=0; j < atr->num; j++){
            free((atr->trpt + j)->add_fld);
         }
         free(atr->trpt);
         atr->num = 0;
       }

    }

    fout = fopen(filout, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", filout);
       exit(1);
    }

    meantrd(fout, tracks, trnum, 's', gpr, ipr, alat, alng);

    fclose(fout);    

    free(ifld);
    free(imnmx);
    free(fthr);
    
    
    for(i=0; i < trnum; i++){   
        atr = tracks + i;
        for(j=0; j < atr->num; j++){
            free((atr->trpt + j)->add_fld);
        }
        free(atr->trpt);
        atr->num = 0;
    }
    free(tracks);

    return 0;
}


int region(struct fet_pt_tr *fp, REG *reg, int loc, int lloc)
{

    int inreg=0;
    
    float xpos=0.0, ypos=0.0;
    
    if(!loc){
       xpos = fp->xf;
       ypos = fp->yf;
    }
    else{
       xpos = fp->add_fld[lloc-2];
       ypos = fp->add_fld[lloc-1];
    }
    
    
    if(xpos > ADD_CHECK) return inreg;
    
    if(! reg->ig){
      if((reg->rlng2 - xpos) * (xpos - reg->rlng1) >= 0.0 &&
         (reg->rlat2 - ypos) * (ypos - reg->rlat1) >= 0.0    ){
         inreg = 1;
      }
    }
    else {
      if((reg->rlat2 - ypos) * (ypos - reg->rlat1) >= 0.0 &&
         ((reg->rlng2 - xpos) * (xpos - 0.0) >= 0.0 || 
         (360.0 - xpos) * (xpos - reg->rlng1) >= 0.0)){
         inreg = 1;
      }

    }

    return inreg;
}

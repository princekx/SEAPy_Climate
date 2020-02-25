#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXCHR  500
#define FIELDWD 1000
#define LTMAX   65

/* program to compute the Siegel statistics */

long int new_time(long int , int );

int main(void)
{

    int i=0;
    
    int ierr=0;
    int nfc=0;
    int ipp=0;
    int step=0;
    int itm=0;
    int ii=0;
    int ifst=0;
    int nem=1;
    int iseed=0;
    int inn=0;
    int iout=0;
    int nfmax=0;
    
    int numpt[LTMAX];
    
    long int dnrb[LTMAX], dnra[LTMAX];
    long int snrb[LTMAX], snra[LTMAX];
    long int bnrb[LTMAX], bnra[LTMAX]; 
    
    long int dnb[LTMAX], dna[LTMAX];
    long int snb[LTMAX], sna[LTMAX];
    long int bnb[LTMAX], bna[LTMAX];    
    long int dabv[LTMAX], sabv[LTMAX], babv[LTMAX];
    
    long int dneff[LTMAX], dtsep[LTMAX];
    long int sneff[LTMAX], stsep[LTMAX];
    long int bneff[LTMAX], btsep[LTMAX];
    
    long int tim1=0, tim2=0, ttmp=0;
    long int fstart=0;
    
    FILE *fmean=NULL;
    FILE *ferr=NULL;
    
    char meanf[MAXCHR];
    char errf[MAXCHR];
    char der[MAXCHR];
    char stub[MAXCHR];
    char line[FIELDWD];
    
    float dmn[LTMAX], smn[LTMAX], bias[LTMAX];
    float dmn_std[LTMAX], smn_std[LTMAX], bias_std[LTMAX];
    float tst[LTMAX];
    
    double dee[LTMAX], see[LTMAX], bee[LTMAX];
    double stdee[LTMAX], stsee[LTMAX], stbee[LTMAX];
    double mmm=0.0, smm=0.0;
    
    float derr=0.0, serr=0.0, berr=0.0;
    
    for(i=0; i < LTMAX; i++){
        dmn[i] = smn[i] = bias[i] = dmn_std[i] = smn_std[i] = bias_std[i]=0.0;
        dnrb[i] = dnra[i] = snrb[i] = snra[i] = bnrb[i] = bnra[i] = 0;
        dnb[i] = dna[i] = snb[i] = sna[i] = bnb[i] = bna[i] = 0;
	dneff[i] = sneff[i] = bneff[i] = 0;
	dee[i] = see[i] = bee[i] = 0.0;
	dabv[i] = sabv[i] = babv[i] = 0;
	dtsep[i] = stsep[i] = btsep[i] = 0;
    }
    
    
    printf("What file contains the means and std statistics?\n\n");
    scanf("%s", meanf);
    
    printf("What type of error file is required:- \r\n"
           "Input '1' for ensembles               \r\n"
	   "Input '2' for ensemble means          \r\n"
	   "Input '3' for control                 \r\n"
	   "Input '4' for deterministic           \r\n"
	   "Input '5' for spread                  \n\n");

    scanf("%d", &ierr);
    
    switch(ierr) {
       case 1:
         sprintf(stub, "%s", "eps_err");
	 break;
       case 2:
         sprintf(stub, "%s", "mean_err");
	 break;
       case 3:
         sprintf(stub, "%s", "cntrl_err");
	 break;
       case 4:
         sprintf(stub, "%s", "det_err");
	 break;
       case 5:
         sprintf(stub, "%s", "sprd_err");
	 break;
       default:
         printf("*****ERROR*****, error type unkown.\n\n");
	 exit(1);
    }
    
    if(ierr < 1 || ierr > 5){
       printf("*****ERROR****, forecast error type not known.\n\n");
       exit(1);    
    }
    
    if(ierr == 1){
       printf("How many ensemble members are there?\n\n");
       scanf("%d", &nem);
    }

    fmean = fopen(meanf, "r");
    if(!fmean){
       printf("****ERROR****, can't open file %s\n", meanf);
       exit(1);
    }

/* read header */

    fgets(line, FIELDWD, fmean); 
    
/* read mean stats */

    if(ierr == 5){
    
       for(i=0; i < LTMAX; i++){
          fgets(line, FIELDWD, fmean);
          sscanf(line, "%f %d %f %f %f %f", tst+i, numpt+i, dmn+i, smn+i, dmn_std+i, smn_std+i); 
       }    

    }
    
    else{

       for(i=0; i < LTMAX; i++){
          fgets(line, FIELDWD, fmean);
          sscanf(line, "%f %d %f %f %f %f %f %f", tst+i, numpt+i, dmn+i, smn+i, bias+i, dmn_std+i, smn_std+i, bias_std+i); 
       }
    
    }

    fclose(fmean);
    
    printf("What is the directory with the error file?\n\n");
    scanf("%s", der);

    printf("What is the maximum number of error files to read?\n\n");
    scanf("%d", &nfmax);
    
    printf("What is the forecast seperation time in hours.\n\n");
    scanf("%d", &step);
    
    
    printf("Outout means and standard errors, '0', or seperation times, '1'\n\n");
    scanf("%d", &iout);
    
    if(iout < 0 || iout > 1){
       printf("****ERROR****, identifier %d is not correct,\n\n", iout);
       exit(1);
    }
    
    
    while(nfc <= nfmax){
       sprintf(errf, "%s/%s%d.dat", der, stub, nfc);

       ferr = fopen(errf, "r");
       if(!ferr){
          printf("****ERROR****, can't open file %s\n", errf);
          ++nfc;
          continue;
       }
       
       
       for(i=0; i < nem; i++){
       
          ifst = 0;
       
          while(fgets(line, FIELDWD, ferr) != NULL){
                if(strstr(line, "TTT")){
	           sscanf(line, "%*s %ld %d", &tim2, &ipp);
		   fgets(line, FIELDWD, ferr);
                   if(strstr(line, "TTT")) continue;
		  
		   if(ipp == i){

	     	      if(!ifst) {
		         ifst = 1;
			 itm = 1;
			 sscanf(line, "%d %f %f %f", &ii, &derr, &serr, &berr);
   
		         if(derr > dmn[ii]) {dabv[ii] = 1; ++dna[ii]; ++dnra[ii];}
		         else {dabv[ii] = 0; ++dnb[ii]; ++dnrb[ii];}
		         if(serr > smn[ii]) {sabv[ii] = 1; ++sna[ii]; ++snra[ii];}
		         else {sabv[ii] = 0; ++snb[ii]; ++snrb[ii];}
			 if(ierr < 5){
		           if(berr > bias[ii]) {babv[ii] = 1; ++bna[ii]; ++bnra[ii];}
		           else {babv[ii] = 0; ++bnb[ii]; ++bnrb[ii];}
			 }

                         tim1 = tim2;

                         continue;
		      }
		      else {
		         itm = 0;
		         ttmp = new_time(tim1, step);
		         if(ttmp == tim2) itm = 1;

                         tim1 = tim2;
		      } 

		   }
	        }
	      
	        if(ipp == i){

	           if(itm){
	              sscanf(line, "%d %f %f %f", &ii, &derr, &serr, &berr);

		      if(derr > dmn[ii]) {
		         ++dna[ii]; 
		         if(!dabv[ii]) {dabv[ii] = 1; ++dnra[ii];}
		      }
		      else {
		         ++dnb[ii]; 
		         if(dabv[ii]) {dabv[ii] = 0; ++dnrb[ii];}
		      }
	              if(serr > smn[ii]) {
		         ++sna[ii]; 
		         if(!sabv[ii]){sabv[ii] = 1; ++snra[ii];}
		      }
		      else {
		         ++snb[ii]; 
		         if(sabv[ii]){sabv[ii] = 0; ++snrb[ii];}
		      }
		      if(ierr < 5){
		         if(berr > bias[ii]) {
		            ++bna[ii]; 
		            if(!babv[ii]){babv[ii] = 1; ++bnra[ii];}
		         }
		         else {
		            ++bnb[ii]; 
		            if(babv[ii]) {babv[ii] = 0; ++bnrb[ii];}	
		         } 
		      }

	           }
		
		   else{
                      sscanf(line, "%d %f %f %f", &ii, &derr, &serr, &berr);

		      if(derr > dmn[ii]) {++dna[ii]; dabv[ii] = 1; ++dnra[ii];}
		      else {++dnb[ii]; dabv[ii] = 0; ++dnrb[ii];}
		      if(serr > smn[ii]) {++sna[ii]; sabv[ii] = 1; ++snra[ii];}
		      else {++snb[ii]; sabv[ii] = 0; ++snrb[ii];}
		      if(berr > bias[ii]) {++bna[ii]; babv[ii] = 1; ++bnra[ii];}
		      else {++bnb[ii]; babv[ii] = 0; ++bnrb[ii];}
		   
		   }
		   
	        }
             
          }
	  
	  fseek(ferr, 0L, SEEK_SET);
      
        }
       
       
        fclose(ferr);
       
        ++nfc;
    
    }

    for(i=0; i < LTMAX; i++){

        if(dna[i] + dnb[i]) {
	   mmm = (double)(2 * dna[i] * dnb[i]);
	   smm = (double)(dna[i] + dnb[i]);
	   dee[i] = (mmm / smm) + 1.0;
	   stdee[i] = sqrt(mmm * (mmm - smm) / smm * smm *(smm -1)); 
	}
	if(sna[i] + snb[i]) {
	   mmm = (double)(2 * sna[i] * snb[i]);
	   smm = (double)(sna[i] + snb[i]);
	   see[i] = (mmm / smm) + 1.0;
	   stsee[i] = sqrt(mmm * (mmm - smm) / smm * smm *(smm -1)); 
	}
	if(bna[i] + bnb[i]) {
	   mmm = (double)(2 * bna[i] * bnb[i]);
	   smm = (double)(bna[i] + bnb[i]);	
	   bee[i] = (mmm / smm) + 1.0;
	   stbee[i] = sqrt(mmm * (mmm - smm) / smm * smm *(smm -1)); 
	}
	
	if(dee[i] > 0.0) dneff[i] = (long int)(((double)(numpt[i] * (dnra[i] + dnrb[i])) / dee[i]) + 0.5);
	if(see[i] > 0.0) sneff[i] = (long int)(((double)(numpt[i] * (snra[i] + snrb[i])) / see[i]) + 0.5);
	if(bee[i] > 0.0) bneff[i] = (long int)(((double)(numpt[i] * (bnra[i] + bnrb[i])) / bee[i]) + 0.5);
	
	if(dneff[i]) dtsep[i] = (long int)(((double)(step * numpt[i]) / (double)dneff[i]) + 0.5);
	if(sneff[i]) stsep[i] = (long int)(((double)(step * numpt[i]) / (double)sneff[i]) + 0.5);
	if(bneff[i]) btsep[i] = (long int)(((double)(step * numpt[i]) / (double)bneff[i]) + 0.5);
	
	if(dtsep[i] < step) {dtsep[i] = step; dneff[i] = numpt[i];}
	if(stsep[i] < step) {stsep[i] = step; sneff[i] = numpt[i];}
	if(btsep[i] < step) {btsep[i] = step; bneff[i] = numpt[i];}
	
/*	printf("--%d %f %ld %ld %ld %ld %ld %ld\n", numpt[i], dee[i], dneff[i], dtsep[i], dnra[i], dnrb[i], dna[i], dnb[i]); */

	
        if(!iout){
	   if(dneff[i]) dmn_std[i] /= sqrt((double)dneff[i]);
	   if(sneff[i]) smn_std[i] /= sqrt((double)sneff[i]);
	   if(bneff[i]) bias_std[i] /= sqrt((double)bneff[i]);
	   
	   printf("%f %f %f %f %f %f %f\n", tst[i], dmn[i], smn[i], bias[i], dmn_std[i], smn_std[i], bias_std[i]);
	} 
	else
	   printf("%5.2f %ld %ld %ld\n", tst[i], dtsep[i], stsep[i], btsep[i]);

        
	
    }
    

    return 0;
}

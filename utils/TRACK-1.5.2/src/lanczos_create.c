#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Math.h>
#include "mem_er.h"
#include "m_values.h"
#include "file_handle.h"
#include "files_out.h"

#define TOLFREQ  1.0e-6
#define BIGPER   1.0e+6

extern char *fext;

extern int iext;

double **lanczos_create(int *nfilt, int *order, int *tsamp)

{
    int i,j, k;

    int nlanc=-1;
    int nor;
    int nfr;


    float *band=NULL, *freq=NULL, *tmp;
    float lowb;
    float samp;
    float fnor;
    float fi, bi;
    float df;
    float rsp;

    double **lweights=NULL, *ww;
    double *sum=NULL;
    double dk, sig, high, low;

    FILE *resp=NULL;
    FILE *fw=NULL;

    char lancr[MAXCHR];
    char lancw[MAXCHR];

    strncpy(lancr, LANCZOS_RESP, MAXCHR);
    strncpy(lancw, LANCZOS_W, MAXCHR);
    if(iext) {
       strcpy(strstr(lancr, EXTENSION), fext);
       strcpy(strstr(lancw, EXTENSION), fext);
    }

    printf("What is the number of data samples per time unit (e.g. day) for period?\n\n");
    scanf("%d", tsamp);

    if(*tsamp < 1){

       printf("****ERROR****, samples per unit time must be positive.\n\n");
       exit(1);

    }

    samp = (float) *tsamp;
    lowb = 2.0 / samp;

    printf("****INFORMATION****, create weights for a Lanczos time domain filter.\n\n");


    while(nlanc < 0){

        printf("How many filters/bands do you require?\n\n");

        scanf("%d", &nlanc);

        *nfilt = nlanc;

        if(nlanc < 0)
           printf("****ERROR***, number of filters must be > 0, try again.\n\n");

    }

    band = (float *)calloc(nlanc+1, sizeof(float));
    mem_er((band == NULL) ? 0 : 1, (nlanc+1) * sizeof(float));

    lweights = (double **)calloc(nlanc, sizeof(double *));
    mem_er((lweights == NULL) ? 0 : 1, nlanc * sizeof(double *));

    sum = (double *)calloc(nlanc, sizeof(double));
    mem_er((sum == NULL) ? 0 : 1, nlanc * sizeof(double));


    freq = (float *)calloc(nlanc+1, sizeof(float));
    mem_er((freq == NULL) ? 0 : 1, (nlanc+1) * sizeof(float));

    printf("What order do you want for the filters?                     \r\n"
           "Choose the order such that no>=1.3/(f2 - f1) for all the    \r\n" 
           "filters, for unit response (gain), where f1 and f2 are the  \r\n"
           "cuttoff frequencies. Also require order < data length.      \n\n");


    scanf("%d", &nor); 

    *order = nor;


/* assign memory for filter coeficients */

    for(i=0; i < nlanc; i++){

        *(lweights + i) = (double *)calloc(nor+1, sizeof(double));
        mem_er((*(lweights + i) == NULL) ? 0 : 1, (nor+1) * sizeof(double));

    }

/* define filter bands by period */

    printf("What are the filter bands in units of period\n\n");

    for(i=0; i<= nlanc; i++){

       printf("Cuttoff %d = ", i+1);
       scanf("%f", band+i);
       printf("\n");
       if(*(band + i) < lowb) *(band + i) = lowb;
       *(freq + i) = 1.0 / (samp * *(band+i));
       printf("Freq.=%f\n\n", *(freq + i));

    }

    for(i=1; i<=nlanc; i++) {

        if(1.3 /( *(freq + i - 1) - *(freq+i)) > (float)nor){

           printf("****WARNING****, the order chosen will not result in unit\r\n"
                  "                 response at the band center for the band\r\n"
                  "                 %f -- %f\r\n"
                  "                 require at least nor=%d\n\n", *(band + i - 1), *(band+i), (int)(1.3 /( *(freq + i - 1) - *(freq+i)))+1);
                 

        }

    }

/* initialize */

    for(i=0; i< nlanc; i++){

      tmp = freq + nlanc - i;

      *(*(lweights + i)) = 2.0 * (*(tmp - 1) - *tmp);
      *(sum + i) = *(*(lweights + i));

    }

    fnor = (float) nor;

/* calculate weights */

    for(i=1; i <= nor; i++){

        fi = (float) i;

        for(j=0; j < nlanc; j++){

           ww = (*(lweights + j) + i);

           tmp = freq + nlanc - j;

           dk = FPI * fi / fnor;
           sig = sin(dk) / dk;
           high = 2.0 * FPI * *(tmp - 1) * fi;
           low = 2.0 * FPI * *tmp * fi;
           *ww = sig * (sin(high) - sin(low)) /(FPI * fi);
           
           *(sum + j) += 2.0 * *ww;

        }

    }
 
    printf("***INFORMATION***, sum of weights for each band.\n\n");

    for(i=0; i< nlanc; i++) printf("% 10.5e ", *(sum + i));
    printf("\n\n");
    
    printf("***INFORMATION***, weights for each band, ordered low pass to high pass.\n\n");

    for(i=0; i<= nor; i++) {

        printf("%3d ", i);
        for(j=0; j < nlanc; j++)printf("% 10.5e ", *(*(lweights + j) + i));
        printf("\n");

    }

/* output weights to file */

    fw = open_file(lancw, "w");


    fprintf(fw, "%d %d %d\n", nlanc, nor, *tsamp);

    for(i=0; i<= nor; i++) {

        fprintf(fw, "%3d ", i);
        for(j=0; j < nlanc; j++)fprintf(fw, "% 10.5e ", *(*(lweights + j) + i));
        fprintf(fw, "\n");

    }

    close_file(fw, lancw);
    

/* output the filter responses to file */

    printf("****INFORMATION****, the sampled response functions will be output, \r\n"
           "                     ordered low pass to high pass, \r\n"
           "                     to the file:- \r\n" 
           "                     %s \r\n"
           "                     in the form:\r\n"
           "                     freq. period r r**2 r**3 r**4 r**10\r\n", lancr);

    printf("Input the number of sampled frequencies required for response.\n\n");
    scanf("%d", &nfr);

    df = 0.5 / (float)(nfr - 1);

    resp = open_file(lancr, "w");

    for(i=0; i < nlanc; i++){

        fprintf(resp, "Band: %f -- %f\n", *(freq + nlanc - i), *(freq + nlanc - i - 1));

        for(j=0; j < nfr; j++){

            rsp = **(lweights + i);
            fi = 2.0 * FPI * (float)j * df;
            if((float)j * df < TOLFREQ) bi = BIGPER;
            else bi = 1.0 / (samp * (float)j * df);

            for(k=1; k <= nor; k++)
                rsp += 2.0 * *(*(lweights + i) + k) * cos(fi * (float) k);


            fprintf(resp, "%e %e %e %e %e %e %e\n", (float)j * df, bi, rsp, rsp*rsp, pow(rsp, 3.0), pow(rsp, 4.0), pow(rsp, 10.0));
 


        }

    }

    close_file(resp, lancr);

    free(sum);
    free(band);
    free(freq);

    return lweights;

}

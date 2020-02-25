#include <stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "splice.h"

#define MAXCHR 500

int main(int argc, char **argv)
{

   int i, j, k;
   int trnum=0, ntheta=0, nr=0, nwfld=0;
   int ntrnum=0, nntheta=0, nnr=0, nnwfld=0;
   int idir=0, nidir=0, ifdir=0, nifdir=0;
   int ndsmth=0, nndsmth;
   int ifcnt=0, itpadd=0, iprojd=0;
   int nifcnt=0, nitpadd=0, niprojd=0;
   int igsmp=0, nigsmp=0;

   int idim=0;

/*   int nfil=0; */

   long int ptnum=0, nptnum=0;

   char freg1[MAXCHR], freg2[MAXCHR], frego[MAXCHR];
   char line[MAXCHR];

   float *ff=NULL, *flng=NULL, *flat=NULL;
   float *ffn=NULL, *flngn=NULL, *flatn=NULL;

   FILE *fin1=NULL, *fin2=NULL, *fout=NULL;

   if(argc != 4){
     printf("Using interactive input.\n");
     printf("Usage:- %s [1st stats file] [ 2nd stats file] [out stats file]\n", argv[0]);
     printf("2nd File - 1st File\n\n");

     printf("What is the first reg file to read?\n\n");
     scanf("%s", freg1);

     printf("What is the second reg file to read?\n\n");
     scanf("%s", freg2);

     printf("What is the output reg file to write?\n\n");
     scanf("%s", frego);

   }
   else{
     sscanf(argv[1], "%s", freg1);
     sscanf(argv[2], "%s", freg2);
     sscanf(argv[3], "%s", frego);
   }


/* open files for read and check headers are consistent */

   fin1 = fopen(freg1, "r");
   if(!fin1){
      printf("****ERROR****, can't open file %s for read.\n\n", freg1);
      exit(1);
   }

   fgets(line, 100, fin1);
   sscanf(line, "%d %ld %d %d %d %d %d %d %d %d %d %d", &trnum, &ptnum, &ntheta, &nr, &nwfld, &idir, &ndsmth, &ifcnt, &itpadd, &iprojd, &ifdir, &igsmp);

   fin2 = fopen(freg2, "r");
   if(!fin2){
      printf("****ERROR****, can't open file %s for read.\n\n", freg2);
      exit(1);
   }

   fgets(line, 100, fin2);
   sscanf(line, "%d %ld %d %d %d %d %d %d %d %d %d %d", &ntrnum, &nptnum, &nntheta, &nnr, &nnwfld, &nidir, &nndsmth, &nifcnt, &nitpadd, &niprojd, &nifdir, &nigsmp);

/* check headers */

   if((trnum != ntrnum) || (ptnum != nptnum) || (ntheta != nntheta) || (nr != nnr) || (nwfld != nnwfld) || 
      (idir != nidir) || (ndsmth != nndsmth) || (itpadd != nitpadd) || (iprojd != niprojd) || (ifdir != nifdir) ||
      (igsmp != nigsmp)) {
       printf("****ERROR****, reg file headers do not match.\n\n");
       exit(1);
   }

   if(ifcnt != nifcnt){
      printf("****WARNING****, using default storm location.\n\n");
      ifcnt = 0;
   }

   fout = fopen(frego, "w");
   if(!fout){
      printf("****ERROR****, can't open file %s for write.\n\n", frego);
      exit(1);
   }

   fprintf(fout, "%6d %10ld %6d %6d %2d %2d %3d %2d %2d %2d %2d %2d\n", trnum, ptnum, ntheta, nr, nwfld, idir, ndsmth, ifcnt, itpadd, iprojd, ifdir, igsmp);

   idim = nr * ntheta;
   ff = (float *)calloc(idim, sizeof(float));
   mem_er((ff == NULL) ? 0 : 1, idim*sizeof(float));
   ffn = (float *)calloc(idim, sizeof(float));
   mem_er((ffn == NULL) ? 0 : 1, idim*sizeof(float));


   flng = (float *)calloc(ntheta, sizeof(float));
   mem_er((flng == NULL) ? 0 : 1, ntheta*sizeof(float));
   flngn = (float *)calloc(ntheta, sizeof(float));
   mem_er((flngn == NULL) ? 0 : 1, ntheta*sizeof(float));

   flat = (float *)calloc(nr, sizeof(float));
   mem_er((flat == NULL) ? 0 : 1, nr*sizeof(float));
   flatn = (float *)calloc(nr, sizeof(float));
   mem_er((flatn == NULL) ? 0 : 1, nr*sizeof(float));

   fread(flng, ntheta * sizeof(float), 1, fin1);
   fscanf(fin1, "%*c");
   fread(flat, nr * sizeof(float), 1, fin1);
   fscanf(fin1, "%*c");
   fread(flngn, ntheta * sizeof(float), 1, fin2);
   fscanf(fin2, "%*c");
   fread(flatn, nr * sizeof(float), 1, fin2);
   fscanf(fin2, "%*c");

   fwrite(flng, ntheta * sizeof(float), 1, fout);
   fprintf(fout, "\n");
   fwrite(flat, nr * sizeof(float), 1, fout);
   fprintf(fout, "\n");

   for(i=0; i < ptnum; i++){
       for(j=0; j < nwfld; j++){
          fread(ff, idim * sizeof(float), 1, fin1);
          fscanf(fin1, "%*c");
          fread(ffn, idim * sizeof(float), 1, fin2);
          fscanf(fin2, "%*c");

          for(k=0; k < idim; k++){
              if(*(ff + k) > ADD_CHECK || *(ffn + k) > ADD_CHECK) *(ffn + k) = ADD_UNDEF;
              else *(ffn + k) -= *(ff + k);
          } 

          fwrite(ffn, idim * sizeof(float), 1, fout);
          fprintf(fout, "\n");
       }

   }

   fclose(fin1);
   fclose(fin2);
   fclose(fout);


   return 0;
}

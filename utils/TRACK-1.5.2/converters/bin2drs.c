#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "drscdf.h"

#define MAXCHR 100

void main(int argc, char *argv[])

{

    int i, j;
    int ix, iy, ntime;
    int dim;
    int ia;
    int ierr=0;
    int ig1, ig2;


    FILE *fin=NULL;
    FILE *drs1=NULL, *drs2=NULL;

    char str[MAXCHR];
    char outfil[MAXCHR];
    char dicf[MAXCHR];


    float *xg=NULL, *yg=NULL;
    float *field=NULL;
    float time;



    if(argc != 3){

       printf("***Usage***, bin2drs [input file] [output file (no extension)]\n\n");

       exit(1);

    }


    strcpy(outfil, argv[2]);
    strcpy(dicf, argv[2]);
    strcat(outfil, ".dat");
    strcat(dicf, ".dic");


    if((drs1=fopen(outfil, "r")) || (drs2=fopen(dicf, "r")) ){

       printf("****ERROR****, DRS files already exist for output filename:- \r\n"
              "               %s \r\n" 
              "               rename or remove files:-                      \r\n"
              "               %s \r\n"
              "               %s \r\n"
              "               and try again.\n\n", argv[2], outfil, dicf);

       if(drs1) fclose(drs1);
       if(drs2) fclose(drs2);
       exit(1);

    }

    fin = fopen(argv[1], "r");

    if(!fin){

       printf("****ERROR***, unable to open file:- \r\n"
              "              %s\n\n", argv[1]);
       exit(1);

    }

    fgets(str, 100, fin);
    sscanf(str, "%d %d %d", &ix, &iy, &ntime);
    dim = ix * iy;

    xg = (float *)calloc(ix, sizeof(float));
    if(!xg){

       printf("****ERROR****, unable to assign memory.\n\n");
       exit(1);

    }

    yg = (float *)calloc(iy, sizeof(float));
    if(!yg){

       printf("****ERROR****, unable to assign memory.\n\n");
       exit(1);

    }

    field = (float *)calloc(dim, sizeof(float));
    if(!field){

       printf("****ERROR****, unable to assign memory.\n\n");
       exit(1);

    }

    ia = 0;
    for(i=0; i < ix; i++) {
        fscanf(fin, "%f", xg+i);
        ++ia;
        if(ia == 10){fgets(str, 100, fin); ia = 0;}
    }
    if(ia) fgets(str, 100, fin);
    ia = 0;
    for(i=0; i < iy; i++) {
        fscanf(fin, "%f", yg+i);
        ++ia;
        if(ia == 10){fgets(str, 100, fin); ia = 0;}
    }
    if(ia) fgets(str, 100, fin);

/* open DRS files */

    if((ierr = Aslun(20, argv[2], 30, argv[2], IDRS_CREATE)) != IDRS_SUCCESS){

       printf("***ERROR***, can't open DRS files for CREATE, error code %d\n", ierr);
       exit(1);

    }

    ierr = Cluvdb();
    ierr = Setname(" ", "longitude", " ", "degrees", " ");
    ierr = Putvdim(20, ix, xg, &ig1, &ig2);

    ierr = Cluvdb();
    ierr = Setname(" ", "latitude", " ", "degrees", " ");
    ierr = Putvdim(20, iy, yg, &ig1, &ig2);

 
    for(i=0; i<ntime; i++){

        time = (float)(i+1);

        fgets(str, 100, fin);
        printf("%s\n", str);

        fread(field, dim*sizeof(float), 1, fin);


        ierr = Cluvdb();
        ierr = Setname(" ", "VAR", " ", " ", "R*4");
        ierr = Setvdim(1, " ", "longitude", " ", " ", *xg, *(xg + ix - 1));
        ierr = Setvdim(2, " ", "latitude", " ", " ", *yg, *(yg + iy - 1));
        ierr = Setdim(3, "time", " ", 1, (float)(i+1), (float)(i+1));

        ierr = Putdat(20, field);

        if(!fgets(str, 100, fin))break;


    }

   if(Cllun(20) != IDRS_SUCCESS){

      printf("***ERROR***, closing DRS files.\n");
      exit(1);

   }


    fclose(fin);

    free(xg);
    free(yg);
    free(field);

    return;

}

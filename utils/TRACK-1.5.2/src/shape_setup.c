#include <Stdio.h>
#include "st_fo.h"
#include "boundary.h"

#define MAXEX 20

/* fourier descriptor setup */


void shape_setup(struct boundary_cntl *bcntl, int fd)

{


    if(fd){

       printf("Do you want Fourier descriptor determination for boundary.\r\n" 
              "Input '1' for yes or '0' for no.                          \r\n"
              "Useful for smoothing object shapes.                       \n\n");
       scanf("%d", &(bcntl->fd));
       if(bcntl->fd < 0 || bcntl->fd > 1){

          printf("****ERROR****, option %d not recognised for fourier descriptors,\r\n"
                 "               defaulting to none.                              \n\n", bcntl->fd);
          bcntl->fd = 0;

       }

    }

    else bcntl->fd = 1;

    if(bcntl->fd > 0){

       bcntl->nexp = 0;
       bcntl->nadd = 0;

       printf("Do you want fixed number of frequencies for all objects, \r\n"
              "or variable number of frequencies (varies with boundary  \r\n"
              "size). For fixed number of frequencies the number chosen \r\n"
              "must be larger than the maximum object boundary number.  \r\n"
              "Input 'f' for fixed, 'v' for variable.                   \n\n");

       scanf("\n");
       if(getchar() == 'f'){

           printf("What is the factor nxp (nf=2^nxp), required for the number of frequencies nf?\n\n");
           scanf("%d", &(bcntl->nexp));
       }

       else{

           printf("What expansion factor is required for number of boundary points,\r\n"
                  "goes as power of 2.                                             \n\n");
           scanf("%d", &(bcntl->nadd));

       }

       printf("How many modes do you want to keep in addition to the zeroth mode?\n\n");
      scanf("%d", &(bcntl->nmode));

       printf("Do you want object fill, i.e. add points to object inside new boundary, 'y' or 'n'\n\n");
       scanf("\n");
       bcntl->ofill = getchar();

       if(bcntl->ofill == 'y'){

          printf("What object expansion is required, in grid points?\n\n");
          scanf("%d", &(bcntl->obex));

          if(bcntl->obex < 0 || bcntl->obex > MAXEX) {

             printf("****WARNING*****, object search expansion out of range,\r\n"
                    "                  defaulting to none.                  \n\n");
             bcntl->obex = 0;

          }

       }

    }

    return;

}

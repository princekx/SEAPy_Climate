#include <Stdio.h>

/* function to determine type of interpolation required */



int query_interp(int it)

{

    int ityp=0;

        
    printf("what kind of interpolation/smoothing is required?             \n\n");

    if(it)
      printf("***INFORMATION***, use option '0' if working on a projection  \n"
             "                   other than cylindrical lat-long.         \n\n");

    printf("input '0' to use smoopy (rect. mesh, no spherical continuity).\n"
           "          Mostely used for whole field or reduced field surface fit.\n"
           "input '1' to use sphery (rect. mesh, spherical continuity).   \n"
           "          Used exclusively for whole sphere fitting.      \n\n");

    scanf("%d", &ityp);



    return ityp;

}

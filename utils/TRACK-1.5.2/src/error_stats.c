#include <Stdio.h>
#include <stdlib.h>


void error_stats(int ier, int line)
{

   printf("***ERROR in STATS data*** for file 'run_means.c', exiting from program due to:-\n\n");

   switch(ier){
      case 1:
         printf("    Incompatable number of grid points, grids don't match.\n\n");
         break;
      case 2:
         printf("    Incompatable regions, data covers different domains.\n\n");
         break;
      case 3:
         printf("    Incompatable kernals, estimates produced with different kernal functions.\n\n");
         break;
      case 4:
         printf("    Incompatable global smoothing parameters.\n\n");
         break;
      case 5:
         printf("    Incompatable grids.\n\n");
         break;
      case 6:
         printf("    Density statistics are of different types.\n\n");
         break;

   }

   printf("    Error occured at line %d\n\n", line);

   exit(1);

   return;

}

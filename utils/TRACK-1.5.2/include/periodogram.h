/* header file for averaged periodogram weights */

#define  COUNT  1

typedef struct ws{
   int m;
   int istart;
   int iend;
   double wsum;
   double *w;
   float wt;
} WS;

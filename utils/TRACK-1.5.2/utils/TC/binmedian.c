#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Some sample C code for the binmedian algorithm. 
 * Note: I haven't written the code for n even yet. */

#define SWAP(a,b) temp=a;a=b;b=temp;

float binmedian(int n, float *x)
{
  int i, j;
  int bottomcount=0;
  int bincounts[1001];
  int bin=0;
  int k=0, r=0, count=0, medbin=0;
  int oldbin=0;
  int samepoints=1;
  
  float sum=0.0;
  float mu=0.0;
  float sigma=0.0;
  float scalefactor=0.0;
  float leftend=0.0;
  float rightend=0.0;
  float oldscalefactor=0.0, oldleftend=0.0;
  float temp=0.0;
  float a=0.0;


  if(n == 1) return x[0];
  
/* Compute the mean and standard deviation */
  sum = 0.0;
  for (i = 0; i < n; i++) {
    sum += x[i];
  }
  mu = sum/n;

  sum = 0.0;
  for (i = 0; i < n; i++) {
    sum += (x[i]-mu)*(x[i]-mu);
  }
  sigma = sqrt(sum/n);

/* Bin x across the interval [mu-sigma, mu+sigma] */

  for (i = 0; i < 1001; i++) {
    bincounts[i] = 0;
  }
  scalefactor = 1000/(2*sigma);
  leftend =  mu-sigma;
  rightend = mu+sigma;

  for (i = 0; i < n; i++) {
    if (x[i] < leftend) {
      bottomcount++;
    }
    else if (x[i] < rightend) {
      bin = (int)((x[i]-leftend) * scalefactor);
      bincounts[bin]++;
    }
  }

/* If n is odd */
  if (n & 1) {
    /* Recursive step */

    k = (n+1)/2;
    r = 0;

    while(1) {
      /* Find the bin that contains the median, and the order */
      /* of the median within that bin                        */
      count = bottomcount;
      for (i = 0; i < 1001; i++) {
        count += bincounts[i];

        if (count >= k) {
          medbin = i;
          k = k - (count-bincounts[i]);
          break;
        }
      }

      bottomcount = 0;
      for (i = 0; i < 1001; i++) {
        bincounts[i] = 0;
      }
      oldscalefactor = scalefactor;
      oldleftend = leftend;
      scalefactor = 1000*oldscalefactor;
      leftend = medbin/oldscalefactor + oldleftend;
      rightend = (medbin+1)/oldscalefactor + oldleftend;

      /* Determine which points map to medbin, and put */
      /* them in spots r,...n-1                        */
      i = r; r = n;
      while (i < r) {
        oldbin = (int)((x[i]-oldleftend) * oldscalefactor);

        if (oldbin == medbin) {
          r--;
          SWAP(x[i],x[r]);

          /* Re-bin on a finer scale */
          if (x[i] < leftend) {
            bottomcount++;
          }
          else if (x[i] < rightend) {
            bin = (int)((x[i]-leftend) * scalefactor);
            bincounts[bin]++;
          }
        }
        else {
          i++;
        }
      }

      /* Stop if all points in medbin are the same  */
      samepoints = 1;
      for (i = r+1; i < n; i++) {
        if (x[i] != x[r]) {
          samepoints = 0;
          break;
        }
      }
      if (samepoints) {
        return x[r];
      }

      /* Stop if there's <= 20 points left */
      if (n-r <= 20) {
        break;
      }
    }

    /* Perform insertion sort on the remaining points, */
    /* and then pick the kth smallest                  */

    for (i = r+1; i < n; i++) {
      a = x[i];
      for (j = i-1; j >= r; j--) {
        if (x[j] > a) {
          break;
        }
        x[j+1] = x[j];
      }
      x[j+1] = a;
    }
    return x[r-1+k];
  }

  /* If n is even (not implemented yet) */
  else {
    printf("****ERROR****, can't handle even number of data points at the moment.\n\n");
    exit(1);
  }
}

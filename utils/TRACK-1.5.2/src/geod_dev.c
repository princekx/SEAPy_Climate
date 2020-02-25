#include <Stdio.h>
#include <Math.h>
#include "st_fo.h"

#define  TOLPHI 0.00000001

struct pt{
    double x;
    double y;
    double z;
};

/* function to compute the geodesic deviation for three true feature points */

double geod_dist(struct feature_pts * , struct feature_pts * );
void geo_conv(struct feature_pts * );

extern int tf;
extern float w1, w2;

double geod_dev(struct feature_pts *fp0, struct feature_pts *fp1, struct feature_pts *fp2)

{

    double alpha1, alpha2;
    double salp1, salp2, dotp, phi;

    struct pt pt1, pt2, t1, t2;

    if(!fp0->gwky) geo_conv(fp0);
    if(!fp1->gwky) geo_conv(fp1);
    if(!fp2->gwky) geo_conv(fp2);


/* compute the angular deviation */

    alpha1 = geod_dist(fp0, fp1);
    alpha2 = geod_dist(fp1, fp2);

    if(alpha1 <= 0. && alpha2 <= 0.) phi = 0.;

    else if(alpha1 <= 0. || alpha2 <= 0.) phi = w2;

    else{

/* compute the norm of the ''tangent'' vectors */

       salp1 = sin(alpha1);
       salp2 = sin(alpha2);

/* compute tangent vectors at central point */

       pt1.x = fp0->gwk[0];
       pt1.y = fp0->gwk[1];
       pt1.z = fp0->gwk[2];

       pt2.x = fp1->gwk[0];
       pt2.y = fp1->gwk[1];
       pt2.z = fp1->gwk[2];

       dotp = pt1.x * pt2.x + pt1.y * pt2.y + pt1.z * pt2.z;

       t1.x = (pt1.x - dotp * pt2.x) / salp1;
       t1.y = (pt1.y - dotp * pt2.y) / salp1;
       t1.z = (pt1.z - dotp * pt2.z) / salp1;


       pt1.x = fp2->gwk[0];
       pt1.y = fp2->gwk[1];
       pt1.z = fp2->gwk[2];


       dotp = pt1.x * pt2.x + pt1.y * pt2.y + pt1.z * pt2.z;

       t2.x = (dotp * pt2.x - pt1.x) / salp2;
       t2.y = (dotp * pt2.y - pt1.y) / salp2;
       t2.z = (dotp * pt2.z - pt1.z) / salp2;

       dotp = t1.x * t2.x + t1.y * t2.y + t1.z * t2.z; 

       phi = w1 * (1.0 - dotp) + w2 * (1.0 - 2.0 * (sqrt(alpha1 * alpha2)/(alpha1 + alpha2)));


    }

    if(fabs(phi) < TOLPHI) phi = 0.;

    return phi;

}

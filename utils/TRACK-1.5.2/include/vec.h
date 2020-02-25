/* Include file for spherical quad tree structure definitions */

#define  TOLVEC    1.0e-6   /* tolerance for vector length check */

typedef struct vec {
    double x;
    double y;
    double z;
} VEC;

/* dot product */

#define  dotp(v1, v2) ((v1)->x * (v2)->x + (v1)->y * (v2)->y + (v1)->z * (v2)->z)

/* vector addition */

#define  addv(v1, v2, v3) {(v3)->x = (v1)->x + (v2)->x; \
                           (v3)->y = (v1)->y + (v2)->y; \
                           (v3)->z = (v1)->z + (v2)->z;}
			   
/* vector subtraction */

#define  subv(v1, v2, v3) {(v3)->x = (v1)->x - (v2)->x; \
                           (v3)->y = (v1)->y - (v2)->y; \
                           (v3)->z = (v1)->z - (v2)->z;}

/* vector normalization */

#define  normv(vv, nrm) {(vv)->x /= (nrm); (vv)->y /= (nrm); (vv)->z /= (nrm);}

/* vector cross product v3=v1Xv2 */ 

#define  crosp(v1, v2, v3) {(v3)->x = (v1)->y * (v2)->z - (v1)->z * (v2)->y; \
                            (v3)->y = (v1)->z * (v2)->x - (v1)->x * (v2)->z; \
                            (v3)->z = (v1)->x * (v2)->y - (v1)->y * (v2)->x;}

/* vector negation */

#define  negate(vv) {(vv)->x *= -1.0; (vv)->y *= -1.0; (vv)->z *= -1.0;}

/* mutiplication by a constant */

#define  mulv(vi, vv, mf) {(vv)->x = (vi)->x * mf; (vv)->y = (vi)->y * mf; (vv)->z = (vi)->z * mf;}

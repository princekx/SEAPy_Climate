/* definitions and macro's for complex arithmetic */

typedef struct comp{         /* definition for a complex number. */
    double real;             /* real part                        */
    double imag;             /* imaginary part                   */
} complex;            


/* complex addition */

#define  cadd(z1, z2, znew) {(znew)->real = (z1).real + (z2).real; \
                             (znew)->imag = (z1).imag + (z2).imag; }

/* complex subtraction */

#define  csub(z1, z2, znew) {(znew)->real = (z1).real - (z2).real; \
                              (znew)->imag = (z1).imag - (z2).imag;} 

/* complex assign */

#define  comp(x, y, znew) {(znew)->real=(x); (znew)->imag=(y);}



/* complex multiply */

#define  cmul(z1, z2, znew) \
       {(znew)->real = (z1).real * (z2).real - (z1).imag * (z2).imag; \
       (znew)->imag = (z1).real * (z2).imag + (z1).imag * (z2).real;}

/* complex conjugate */

#define conjg(z, znew) {(znew)->real = (z).real; (znew)->imag = -(z).imag;}

/* modulas */

#define cabs(z) ( sqrt((z).real * (z).real + (z).imag * (z).imag) )

/* modulas squared */

#define cabs2(z) ( (z).real * (z).real + (z).imag * (z).imag )

/* multiplication by a real */

#define cmx(x, zz, znew) {(znew)->real = (x) * (zz).real; \
                         (znew)->imag = (x) * (zz).imag;}
 
/* complex copy */

#define ccopy(z1, z2) {(z2)->real = (z1)->real; (z2)->imag = (z1)->imag;}


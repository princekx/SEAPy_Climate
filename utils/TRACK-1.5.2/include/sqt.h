/* Include file for spherical quad tree structure definitions */

/* include vector definitions */

#include "vec.h"

typedef struct sqt *TPTR;
typedef struct leaf *LPTR;

typedef struct sqt {
   char ud;
   TPTR parent;
   TPTR tt[4];
   int ivec[3];
} SQT;


typedef struct leaf {
   char ud;
   TPTR parent;
   LPTR lft;
   LPTR rgt;
   LPTR udw;
   int ingh;         /* neighbour flag */
   int ivec[3];
   int inrm[3];
   int ndata;
   int *ldata;
   int ng;
   int *lgrid;
   float area;
} LEAF;

/* Base vertices of polyhedra (octahedron and icosahedron) for NH
   SH obtained by reflection.                                      */

/* Octahedron */

#define XP  {  1.0,  0.0,  0.0 }
#define YP  {  0.0,  1.0,  0.0 }
#define ZP  {  0.0,  0.0,  1.0 }
#define XPM { -1.0,  0.0,  0.0 }
#define YPM {  0.0, -1.0,  0.0 }
#define ZPM {  0.0,  0.0, -1.0 }


/* Icosahedron */

#define NP { 0.0,       0.0,       1.0      }
#define NA { 0.0,       0.894427,  0.447214 }
#define NB {-0.850651,  0.276393,  0.447214 }
#define NC {-0.525731, -0.723607,  0.447214 }
#define ND { 0.525731, -0.723607,  0.447214 }
#define NE { 0.850651,  0.276393,  0.447214 }
#define SP { 0.0,       0.0,      -1.0      }
#define SA { 0.0,      -0.894427, -0.447214 }
#define SB { 0.850651, -0.276393, -0.447214 }
#define SC { 0.525731,  0.723607, -0.447214 }
#define SD {-0.525731,  0.723607, -0.447214 }
#define SE {-0.850651, -0.276393, -0.447214 }

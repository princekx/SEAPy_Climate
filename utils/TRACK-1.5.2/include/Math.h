#if defined(CYGWIN) || defined(MAC)
#include "/usr/include/math.h"
#else
#include <math.h>
#endif

#if defined(SUNOS5)

#include <sunmath.h>

#elif !defined(SUNOS4)

void sincos(double , double * , double * );

#endif

#include <Stdio.h>
#include <stdlib.h>

long int julian(long int tim)
{
    long int jt=0;
    
    long int year=0, mnth=0, day=0, hr;
    
    year = tim / 1000000;
    mnth = (tim - year * 1000000) / 10000;
    day =  (tim - year * 1000000 - mnth * 10000) / 100;
    hr = tim - year * 1000000 - mnth * 10000 - 100 * day;
    
    jt = (day-32075+1461*(year+4800+(mnth-14)/12)/4 + 367*(mnth-2-((mnth-14)/12)*12)/12 - 3*((year + 4900 + (mnth - 14)/12)/100)/4) * 24 + hr;

    return jt;
}


#include <Stdio.h>
#include <Math.h>

double hoskins_filt(float n, double cutf)

{

    double hfilt;
    double wn, wn1;

    wn = (double) n;
    wn1 = wn + 1.0;

    wn *= wn;
    wn1 *= wn1;

    hfilt = exp(-cutf * wn * wn1);

    return hfilt;

}

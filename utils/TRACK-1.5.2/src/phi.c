#include "zones.h"
  

/* function to compute values for phimax based on point seperation */


float phi(float d1, float d2, ADPT *add)

{

    float dd = 0.5 * (d1 + d2);

    if(dd < add->ct[0]) return add->phii[0];

    else if(dd >= add->ct[3]) return add->phii[3];

    else if(dd >= add->ct[0] && dd < add->ct[1])
         return add->psl[0] * dd + add->incp[0];
	 
    else if(dd >= add->ct[1] && dd < add->ct[2])
         return add->psl[1] * dd + add->incp[1];

    else return add->psl[2] * dd + add->incp[2];

}

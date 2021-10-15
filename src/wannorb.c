#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "constants.h"
#include "vector.h"
#include "wannorb.h"

//#define __DEBUG

int find_index_of_wannorb(wannorb * wann,int num_wann, vector site, int r, int l, int mr, int ms)
{
    int i,j,k;
    //num_wann = sizeof(wann)/sizeof(wannorb);
    for(i=0;i<num_wann;i++){
        //seq ms mr l r
        if( (wann+i)->ms == ms && (wann+i)->mr == mr && (wann+i)->l == l && (wann+i)->r == r){
            if( equal((wann+i)->site,site)){
                return i;
           }
        }

    }
    return -1;
}

void init_wannorb(wannorb * orb,vector site, int l, int mr, int ms, int r, vector axisz, vector axisx){
    orb->site    = site;
    orb->l       = l;
    orb->mr      = mr;
    orb->ms      = ms;
    orb->r       = r;
    orb->axis[2] = axisz;
    orb->axis[0] = axisx;
    orb->axis[1] = cross_product(axisz, axisx);
}

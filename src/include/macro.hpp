

#pragma once 


#define real double
#define ind int



#include <math.h>

#include <stdio.h>



#define GAMMA 1.4

enum BndType{
    PERIODIC,
    DIRICLET_SODL,
    DIRICLET_SODR,
};

#define LEFTT 0
#define RIGHT 1
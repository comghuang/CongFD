

#pragma once 


#define real double
#define ind int



#include <math.h>

#include <stdio.h>
#include <map>
#include <string>



#define GAMMA 1.4

enum BndType{
    TYPENULL,
    PERIODIC,
    DIRICLET_SODL,
    DIRICLET_SODR,
};
enum SpaceDisMethod{
    FIRSTORDER,
    MUSCL,
    WCNSJS5
};

enum FluxType{
    LINEARCONV1D,
    BURGERS1D,
    EULER1D
};
#define LEFTT 0
#define RIGHT 1




#pragma once 


#define real double
#define ind int



#include <math.h>

#include <stdio.h>
#include <map>
#include <string>
#include <memory>
#include <iostream>
#include <cstring>
#include <array>
#include <fstream>
#include <vector>
#include <cgnslib.h>
#include <map>
#include <algorithm>



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


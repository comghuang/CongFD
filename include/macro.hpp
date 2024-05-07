

#pragma once 


#define real double
#define ind int


#include <iostream>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <stdio.h>
#include <cstring>
#include <array>

#define GAMMA 1.4

enum BndType{
    PERIODIC,
    DIRICLET_SODL,
    DIRICLET_SODR,
};

#define LEFTT 0
#define RIGHT 1
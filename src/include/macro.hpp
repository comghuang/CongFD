

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
#include <format>



#define GAMMA 1.4

enum BndType{
    TYPENULL,
    PERIODIC1D,
    DIRICLET_SODL,
    DIRICLET_SODR,
    FLUXGHOST,
    SUPERSONICOUTLET
};
enum InterMethod{
    FIRSTORDER,
    MUSCL,
    WCNSJS5
};
enum DiffMethod{
    HDS6,
    TRAD6,
    TRAD2
};

enum EquationType{
    LINEARCONV1D,
    BURGERS1D,
    EULER1D
};

enum TimeMethod{
    RK3SSP,
    EulerFront};



#define LEFTT 0
#define RIGHT 1

int index(int,int,int,std::array<int,3>);
std::array<int,2> calOffset(int dim,int i,int j,std::array<int,3>);
std::array<int,2> calOffsetInverse(int idim,int i,int j,std::array<int,3> iMax);
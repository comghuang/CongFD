

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
    EULER,
};

enum FluxMethod{
    HLLC1D,
    ROE1D,
    HLLC2D
};

enum TimeMethod{
    RK3SSP,
    EulerFront};





// int index(int,int,int,std::array<int,3>);
// std::array<int,2> calOffset(int dim,int i,int j,std::array<int,3>);
// std::array<int,2> calOffsetInverse(int idim,int i,int j,std::array<int,3> iMax);

constexpr int index(int i,int j,int k,std::array<int,3> iMax)
{
    return i+j*iMax[0]+k*iMax[0]*iMax[1];
}

constexpr std::array<int,2> calOffset(int idim,int i,int j,std::array<int,3> iMax)
{
    std::array<int,3> offsets{1,iMax[0],iMax[0]*iMax[1]};
    std::array<int,2> res;
    if(idim==1)
    {
        res[0]=i*offsets[1]+j*offsets[2];//i0
        res[1]=offsets[0];//offset
    }
    else if(idim==2)
    {
        res[0]=i*offsets[0]+j*offsets[2];//i0
        res[1]=offsets[1];//offset
    }
    else if (idim==3)
    {
        res[0]=i*offsets[0]+j*offsets[1];//i0
        res[1]=offsets[2];//offset
    }
    return res;
}

constexpr std::array<int,2> calOffsetInverse(int idim,int i,int j,std::array<int,3> iMax)
{
    std::array<int,3> offsets{1,iMax[0],iMax[0]*iMax[1]};
    std::array<int,2> res;
    if(idim==1)
    {
        res[0]=i*offsets[1]+j*offsets[2]+(iMax[0]-1)*offsets[0];//i0
        res[1]=-offsets[0];//offset
    }
    else if(idim==2)
    {
        res[0]=i*offsets[0]+j*offsets[2]+(iMax[1]-1)*offsets[1];//i0
        res[1]=-offsets[1];//offset
    }
    else if (idim==3)
    {
        res[0]=i*offsets[0]+j*offsets[1]+(iMax[2]-1)*offsets[2];//i0
        res[1]=-offsets[2];//offset
    }
    return res;
}
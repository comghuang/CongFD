#include "macro.hpp"

int index(int i,int j,int k,std::array<int,3> iMax)
{
    return i+j*iMax[0]+k*iMax[0]*iMax[1];
}
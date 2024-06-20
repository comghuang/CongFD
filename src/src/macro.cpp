#include "macro.hpp"

int index(int i,int j,int k,std::array<int,3> iMax)
{
    return i+j*iMax[0]+k*iMax[0]*iMax[1];
}

std::array<int,2> calOffset(int idim,int i,int j,std::array<int,3> iMax)
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

std::array<int,2> calOffsetInverse(int idim,int i,int j,std::array<int,3> iMax)
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
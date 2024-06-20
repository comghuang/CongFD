#include"sp_distributor.hpp"

void SpDistributor::rhsSolve()
{
    if (dim<=0||dim>3)
    {
        std::cout<<"SpDistributor error: dim error";
        return;
    }
    
    if (dim<=1)
    {
        for (int i = 0; i < iMax[1]; i++)
        for (int j = 0; j < iMax[2]; j++)
        {
            auto oneDBnds=bnds->getOneDBnd(1,i,j);
            auto offsets=calOffset(1,i,j);
            SpaceDis spDis(iMax[0],prim,rhs,oneDBnds[0],oneDBnds[1],info);
            spDis.difference();
        }
    }
    
}

std::array<int,2> SpDistributor::calOffset(int idim,int i,int j)
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
    
    
    
}
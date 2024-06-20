#include"sp_distributor.hpp"

void SpDistributor::rhsSolve()
{
    if (dim<=0||dim>3)
    {
        std::cout<<"SpDistributor error: dim error";
        return;
    }
    
    if (dim>=1)
    {
        for (int i = 0; i < iMax[1]; i++)
        for (int j = 0; j < iMax[2]; j++)
        {
            auto oneDBnds=bnds->getOneDBnd(1,i,j);
            auto offsets=calOffset(1,i,j,iMax);
            SpaceDis spDis(iMax[0],prim,rhs,oneDBnds[0],oneDBnds[1],info);
            spDis.difference();
        }
    }

    if (dim>=2)
    {
        for (int i = 0; i < iMax[0]; i++)
        for (int j = 0; j < iMax[2]; j++)
        {
            auto oneDBnds=bnds->getOneDBnd(2,i,j);
            auto offsets=calOffset(2,i,j,iMax);
            SpaceDis spDis(iMax[1],prim,rhs,oneDBnds[0],oneDBnds[1],info);
            spDis.difference();
        }
    }

    if (dim>=3)
    {
        for (int i = 0; i < iMax[0]; i++)
        for (int j = 0; j < iMax[1]; j++)
        {
            auto oneDBnds=bnds->getOneDBnd(3,i,j);
            auto offsets=calOffset(3,i,j,iMax);
            SpaceDis spDis(iMax[2],prim,rhs,oneDBnds[0],oneDBnds[1],info);
            spDis.difference();
        }
    }
    
}

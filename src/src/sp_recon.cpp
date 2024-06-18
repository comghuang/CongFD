#include "SpaceDis.hpp"
#include "interScheme.hpp"

real SpaceDis::reconL(ind i,ind ivar)
{
    if(spDisMethod==MUSCL)
    {
        real delta;
        real deltam,deltap;
        deltam=(*data)(i-1,ivar)-(*data)(i-2,ivar);
        deltap=(*data)(i,ivar)-(*data)(i-1,ivar);
        
        //minmod
        real beta=1.0;
        if (deltap>0)
        {
            delta=std::max(0.0,std::max(std::min(beta*deltam,deltap),std::min(deltam,beta*deltap)));
        }
        else
        {
            delta=std::min(0.0,std::min(std::max(beta*deltam,deltap),std::max(deltam,beta*deltap)));
        }
        return (*data)(i-1,ivar)+delta*0.5;
    }
    else if(spDisMethod==WCNSJS5)
    {
        return weno5_JSchen((*data)(i-3,ivar),(*data)(i-2,ivar),(*data)(i-1,ivar),(*data)(i,ivar),(*data)(i+1,ivar));
    }
    else return (*data)(i-1,ivar);
}

real SpaceDis::reconR(ind i,ind ivar)
{
    if(spDisMethod==MUSCL)
    {
        real delta;
        real deltam,deltap;
        deltam=(*data)(i,ivar)-(*data)(i-1,ivar);
        deltap=(*data)(i-1,ivar)-(*data)(i,ivar);
        
        //minmod
        real beta=1.0;
        if (deltap>0)
        {
            delta=std::max(0.0,std::max(std::min(beta*deltam,deltap),std::min(deltam,beta*deltap)));
        }
        else
        {
            delta=std::min(0.0,std::min(std::max(beta*deltam,deltap),std::max(deltam,beta*deltap)));
        }
        return (*data)(i,ivar)-delta*0.5;
    }
    else if(spDisMethod==WCNSJS5)
    {
        /*WENO-JS 5 order*/
        return weno5_JSchen((*data)(i+2,ivar),(*data)(i+1,ivar),(*data)(i,ivar),(*data)(i-1,ivar),(*data)(i-2,ivar));
    }
    else return (*data)(i,ivar);

}

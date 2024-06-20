#include "SpaceDis.hpp"
#include "interScheme.hpp"

real SpaceDis::reconL(ind i,ind ivar)
{
    if(interMethod==MUSCL)
    {
        real delta;
        real deltam,deltap;
        deltam=at(i-1,ivar)-at(i-2,ivar);
        deltap=at(i,ivar)-at(i-1,ivar);
        
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
        return at(i-1,ivar)+delta*0.5;
    }
    else if(interMethod==WCNSJS5)
    {
        return weno5_JSchen(at(i-3,ivar),at(i-2,ivar),at(i-1,ivar),at(i,ivar),at(i+1,ivar));
    }
    else return at(i-1,ivar);
}

real SpaceDis::reconR(ind i,ind ivar)
{
    if(interMethod==MUSCL)
    {
        real delta;
        real deltam,deltap;
        deltam=at(i,ivar)-at(i-1,ivar);
        deltap=at(i-1,ivar)-at(i,ivar);
        
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
        return at(i,ivar)-delta*0.5;
    }
    else if(interMethod==WCNSJS5)
    {
        /*WENO-JS 5 order*/
        return weno5_JSchen(at(i+2,ivar),at(i+1,ivar),at(i,ivar),at(i-1,ivar),at(i-2,ivar));
    }
    else return at(i,ivar);

}

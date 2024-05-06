#include "SpaceDis.hpp"
#include <array>
#include "interScheme.hpp"

SpaceDis::SpaceDis(Data* data_,Data* coor_,ind n_,ind nVar_)
{
    coor=coor_;
    data=data_;
    n=n_;
    nVar=nVar_;
    nHalf=n;
    flux.init(nHalf*nVar,nVar);
}

void SpaceDis::calFluxConv()
{
    /*for u_t + a * u_x == 0*/
    real a=1.0;
    real ul,ur;

    for(ind i=0;i<nHalf;i++)
    {
        
        ul=reconL(i,0);
        ur=reconR(i,0);
        
        /*
        ul=data(i-1,0);
        ur=data(i,0);*/
        /*if (a>0)
        {
            flux(i,0)=a*ul;
        }
        else
        {
            flux(i,0)=a*ur;
        }*/

        flux(i,0)=0.5*(a*ul+a*ur-abs(a)*(ur-ul));
    }
}

void SpaceDis::calFluxBurgers()
{
    /*for u_t + a * u_x == 0*/
    
    real ul,ur,aLF;
    aLF=data->maxElement(0);

    for(ind i=0;i<nHalf;i++)
    {
        
        ul=reconL(i,0);
        ur=reconR(i,0);
        
        
        //ul=data(i-1,0);
        //ur=data(i,0);
        //L-F flux
        real al,ar;
        al=std::abs(ul);
        ar=std::abs(ur);
        //real a=std::max(ar,ar);
        //real a=(al+ar)/2;
        real a=aLF;
        flux(i,0)=0.5*(ul*ul/2+ur*ur/2-a*(ur-ul));
    }
}

std::vector<real> SpaceDis::difference()
{
    calFlux();
    std::vector<real> rhs;
    rhs.resize(n*nVar);
    for(ind i=0;i<n;i++)
    {
        real h;
        h=2.0/n;
        for(ind j=0;j<nVar;j++)
        {
            rhs[i*nVar+j]=75.0/64.0*(flux(i+1,j)-flux(i,j))/h
                         -25.0/128.0*(flux(i+2,j)-flux(i-1,j))/(3*h)
                         +3.0/128.0*(flux(i+3,j)-flux(i-2,j))/(5*h);
        }
    }
    return rhs;
}

void SpaceDis::calFlux()
{
    calFluxBurgers();
}

real SpaceDis::reconL(ind i,ind ivar)
{
    ind weight=1;
    if(weight==0)
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
    else if(weight==1)
    {
        return weno5_JSchen((*data)(i-3,ivar),(*data)(i-2,ivar),(*data)(i-1,ivar),(*data)(i,ivar),(*data)(i+1,ivar));
    }
    else return (*data)(i-1,ivar);
}

real SpaceDis::reconR(ind i,ind ivar)
{
    ind weight=1;
    if(weight==1)
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
    else if(weight==1)
    {
        /*WENO-JS 5 order*/
        return weno5_JSchen((*data)(i+2,ivar),(*data)(i+1,ivar),(*data)(i,ivar),(*data)(i-1,ivar),(*data)(i-2,ivar));
    }
    else return (*data)(i-1,ivar);

}


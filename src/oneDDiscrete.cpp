#include "oneDDiscrete.hpp"

oneDDiscrete::oneDDiscrete(std::vector<real>* data_,ind n_,ind nVar_)
{
    data.accept(data_,nVar_);
    n=n_;
    nVar=nVar_;
    nHalf=n+1;
    flux.allocate(nHalf*nVar,nVar);
}

void oneDDiscrete::calFluxConv()
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
        if (a>0)
        {
            flux(i,0)=a*ul;
        }
        else
        {
            flux(i,0)=a*ur;
        }
    }
}

void oneDDiscrete::calFluxBurgers()
{
    /*for u_t + a * u_x == 0*/
    
    real ul,ur;

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
        real a=(al+ar)/2;
        flux(i,0)=0.5*(ul*ul/2+ur*ur/2-a*(ur-ul));
    }
}

std::vector<real> oneDDiscrete::difference()
{
    calFlux();
    std::vector<real> rhs;
    rhs.resize(n*nVar);
    for(ind i=0;i<n;i++)
    for(ind j=0;j<nVar;j++)
    {
        rhs[i*nVar+j]=(flux(i+1,j)-flux(i,j))/(2.0/n);
    }
    return rhs;
}

void oneDDiscrete::calFlux()
{
    calFluxBurgers();
}

real oneDDiscrete::reconL(ind i,ind ivar)
{
    real delta;
    real deltam,deltap;
    deltam=data(i-1,ivar)-data(i-2,ivar);
    deltap=data(i,ivar)-data(i-1,ivar);
    
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
    return data(i-1,ivar)+delta*0.5;

}

real oneDDiscrete::reconR(ind i,ind ivar)
{
    real delta;
    real deltam,deltap;
    deltam=data(i,ivar)-data(i-1,ivar);
    deltap=data(i-1,ivar)-data(i,ivar);
    
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
    return data(i,ivar)-delta*0.5;

}
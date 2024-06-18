#include "SpaceDis.hpp"

void SpaceDis::calFluxConv(ind i)
{
    /*for u_t + a * u_x == 0*/
    real a=1.0;
    real ul,ur;
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

void SpaceDis::calFluxBurgers(ind i)
{
    /*for u_t + a * u_x == 0*/
    
    real ul,ur,aLF;
    aLF=data->maxElement(0);

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

void SpaceDis::calFluxEuler(ind i)
{
    /*for u_t + a * u_x == 0*/
    
    real ul,ur,aLF;
    aLF=data->maxElement(0);
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
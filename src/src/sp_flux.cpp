#include "SpaceDis.hpp"
#include "fluxScheme.hpp"

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

    fAt(i,0)=0.5*(a*ul+a*ur-abs(a)*(ur-ul));
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
    fAt(i,0)=0.5*(ul*ul/2+ur*ur/2-a*(ur-ul));
}
typedef std::array<real,2> arr2;
void SpaceDis::calFluxEuler(ind i)
{
    /*for u_t + a * u_x == 0*/
    arr2 r={reconL(i,0),reconR(i,0)};
    arr2 u={reconL(i,1),reconR(i,1)};
    arr2 p={reconL(i,2),reconR(i,2)};
    arr2 H={reconL(i,3),reconR(i,3)};
    arr2 RT={reconL(i,4),reconR(i,4)};
    std::vector<real> iflux=roeFlux1D2(r,u,p,H,RT);
    //std::vector<real> iflux=HLLCFlux1D(r,u,p,H,RT);
    //std::vector<real> iflux2=roeFlux1D(r,u,p,H,RT);
    for (ind ivar = 0; ivar < 3; ivar++)
    {
        fAt(i,ivar)=iflux[ivar];
    }
    
}
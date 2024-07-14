#include "SpaceDis.hpp"
#include "fluxScheme.hpp"
#include "interScheme.hpp"

void SpaceDis::calFluxConv(int i)
{
    /*for u_t + a * u_x == 0*/
    real a=1.0;
    real ul,ur;
    ul=(this->*reconLMethod)(i).at(0);
    ur=(this->*reconRMethod)(i).at(0);
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

    fluxAt(i,0)=0.5*(a*ul+a*ur-abs(a)*(ur-ul));
}
void SpaceDis::calFluxAccuracyTest(int i)
{
    /*for u_t + a * u_x == 0*/
    real a=1.0;
    real ul,ur;
    ul=(this->*reconLMethod)(i).at(0);

    fluxAt(i,0)=ul;
}

void SpaceDis::calFluxBurgers(int i)
{
    /*for u_t + a * u_x == 0*/
    
    real ul,ur,aLF;
    aLF=data->maxElement(0);

    ul=(this->*reconLMethod)(i).at(0);
    ur=(this->*reconRMethod)(i).at(0);
    
    
    //ul=data(i-1,0);
    //ur=data(i,0);
    //L-F flux
    real al,ar;
    al=std::abs(ul);
    ar=std::abs(ur);
    //real a=std::max(ar,ar);
    //real a=(al+ar)/2;
    real a=aLF;
    fluxAt(i,0)=0.5*(ul*ul/2+ur*ur/2-a*(ur-ul));
}


typedef std::array<real,2> arr2;
void SpaceDis::calFluxEuler1DHLLC(int i)
{
    /*for u_t + a * u_x == 0*/
    auto WL=(this->*reconLMethod)(i);
    auto WR=(this->*reconRMethod)(i);
    real cl=(WL[1]-WL[2])/2;
    real ul=(WL[1]+WL[2])/2;
    real rl=WL[0];
    real pl=cl*cl/GAMMA*rl;

    real cr=(WR[1]-WR[2])/2;
    real ur=(WR[1]+WR[2])/2;
    real rr=WR[0];
    real pr=cr*cr/GAMMA*rr;

    //auto iflux=roeFlux1D2(rl,rr,ul,ur,pl,pr);
    auto iflux=HLLCFlux1D(rl,rr,ul,ur,pl,pr);
    //std::vector<real> iflux2=roeFlux1D(r,u,p,H,RT);
    for (int ivar = 0; ivar < 3; ivar++)
    {
        fluxAt(i,ivar)=iflux[ivar];
    }
    
}

void SpaceDis::calFluxEuler1D(int i)
{
    /*for u_t + a * u_x == 0*/
    auto WL=(this->*reconLMethod)(i);
    auto WR=(this->*reconRMethod)(i);
    //std::vector<real> iflux=roeFlux1D2(WL[0],WR[0],WL[1],WR[1],WL[2],WR[2]);

    auto iflux=HLLCFlux1D(WL[0],WR[0],WL[1],WR[1],WL[2],WR[2]);
    //std::vector<real> iflux2=roeFlux1D(r,u,p,H,RT);
    for (int ivar = 0; ivar < 3; ivar++)
    {
        fluxAt(i,ivar)=iflux[ivar];
    }
}

void SpaceDis::calFluxEuler2D(int i)
{
    /*for u_t + a * u_x == 0*/
    auto WL=(this->*reconLMethod)(i);
    auto WR=(this->*reconRMethod)(i);
    if (WL[3]<0||WL[0]<0)
    {
        //std::cout<<std::format("SpaceDis error: negative pressure WL={} WR={} Wi={}\n",WL[3],WR[3],at(i-1,3));
        WL={at(i-1,0),at(i-1,1),at(i-1,2),at(i-1,3)};
    }
    if (WR[3]<0||WR[0]<0)
    {
        //std::cout<<std::format("SpaceDis error: negative pressure WL={} WR={} Wi={}\n",WL[3],WR[3],at(i,3));
        WR={at(i,0),at(i,1),at(i,2),at(i,3)};
    }
    auto iflux=roeFlux2D(WL[0],WR[0],WL[1],WR[1],WL[2],WR[2],WL[3],WR[3],norm);
    //auto iflux=HLLCFlux2D(WL[0],WR[0],WL[1],WR[1],WL[2],WR[2],WL[3],WR[3],norm);
    //std::vector<real> iflux2=roeFlux1D(r,u,p,H,RT);
    for (int ivar = 0; ivar < 4; ivar++)
    {
        fluxAt(i,ivar)=iflux[ivar];
    }
}

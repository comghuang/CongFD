#include "SpaceDis.hpp"
#include <array>


SpaceDis::SpaceDis(Data* data_,Data* coor_,ind n_,ind nVar_)
{
    coor=coor_;
    data=data_;
    n=n_;
    nVar=nVar_;
    nHalf=n+1;
    flux.init(nHalf*nVar,nVar);
    flux.setGhostVertex(2,LEFTT);
    flux.setGhostVertex(2,RIGHT);
    data->setGhostVertex(5,LEFTT);
    data->setGhostVertex(5,RIGHT);
}



std::vector<real> SpaceDis::difference()
{
    data->updateGhostVertex();
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
    std::array<ind,2> nGhost=flux.getNGhost();
    for(ind i=0-nGhost[LEFTT];i<nHalf+nGhost[RIGHT];i++)
    {
        switch (fluxType)
        {
        case LINEARCONV:
            calFluxConv(i);
            break;
        case BURGERS:
            calFluxBurgers(i);
            break;
        case EULER:
            calFluxEuler(i);
            break;
        
        default:
            break;
        }
    }
}



void SpaceDis::setMethod(SpaceDisMethod method_,FluxType type_)
{
    fluxType=type_;
    spDisMethod=method_;
}
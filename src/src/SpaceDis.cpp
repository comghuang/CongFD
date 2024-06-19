#include "SpaceDis.hpp"
#include <array>


SpaceDis::SpaceDis(){};


void SpaceDis::init(std::shared_ptr<Data> data_,std::shared_ptr<Block> grid_,ind n_,ind nVar_,ind nPrim_)
{
    grid=grid_;
    data=data_;
    nPrim=nPrim_;
    n=n_;
    nVar=nVar_;
    nHalf=n+1;
    flux.init(nHalf*nVar,nVar);
    OneDBnd* fluxl=new OneDBnd;
    OneDBnd* fluxr=new OneDBnd;
    OneDBnd* datal=new OneDBnd;
    OneDBnd* datar=new OneDBnd;
    fluxl->init(2,nVar,TYPENULL);
    fluxr->init(2,nVar,TYPENULL);
    datal->init(5,nPrim,DIRICLET_SODL);
    datar->init(5,nPrim,DIRICLET_SODR);
    flux.setGhostVertex(fluxl,fluxr);
    data->setGhostVertex(datal,datar);
}



std::vector<real> SpaceDis::difference()
{
    data->updateGhostVertex();
    calFlux();
    return (this->*difMethod)();
}

void SpaceDis::calFlux()
{
    std::array<ind,2> nGhost=flux.getNGhost();
    flux.setZeros();
    for(ind i=0-nGhost[LEFTT];i<nHalf+nGhost[RIGHT];i++)
    {
        (this->*calTypeFlux)(i);
    }
}



void SpaceDis::setMethod(SpaceDisMethod method_,EquationType type_)
{
    fluxType=type_;
    spDisMethod=method_;
    switch (fluxType)
    {
    case LINEARCONV1D:
        calTypeFlux=&SpaceDis::calFluxConv;
        break;
    case BURGERS1D:
        calTypeFlux=&SpaceDis::calFluxBurgers;
        break;
    case EULER1D:
        calTypeFlux=&SpaceDis::calFluxEuler;
        break;
    
    default:
        break;
    }

    switch (spDisMethod)
    {
    case FIRSTORDER:
        difMethod=&SpaceDis::dif2Order;
        break;
    case MUSCL:
        difMethod=&SpaceDis::dif2Order;
        break;
    case WCNSJS5:
        difMethod=&SpaceDis::difTraditional6;
        break;
    
    default:
        break;
    }
}
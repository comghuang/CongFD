#include "SpaceDis.hpp"


SpaceDis::SpaceDis(){};
SpaceDis::SpaceDis(int n_,std::shared_ptr<Data> data_,std::shared_ptr<Data> rhs_
            ,std::shared_ptr<OneDBnd> bndL_,std::shared_ptr<OneDBnd> bndR_,std::shared_ptr<Info> info)
{
    data=data_;
    nPrim=info->nPrim();
    n=n_;
    nVar=info->nCons();
    nHalf=n+1;
    flux=std::make_shared<Data>();
    fBndL=std::make_shared<OneDBnd>(info->nFluxPoint(),nVar,FLUXGHOST);
    fBndR=std::make_shared<OneDBnd>(info->nFluxPoint(),nVar,FLUXGHOST);

    fluxType=info->eqType;
    diffMethod=info->diffMethod;
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

    switch (diffMethod)
    {
    case TRAD2:
        difMethod=&SpaceDis::dif2Order;
        break;
    case TRAD6:
        difMethod=&SpaceDis::difTraditional6;
        break;
    case HDS6:
        difMethod=&SpaceDis::difHCS;
        break;
    
    default:
        break;
    }
}

void SpaceDis::setOffset(int i0_,int offset_)
{
    i0=i0_;offset=offset_;
}



std::vector<real> SpaceDis::difference()
{
    data->updateGhostVertex();
    calFlux();
    return (this->*difMethod)();
}

void SpaceDis::calFlux()
{
    std::array<ind,2> nGhost=flux->getNGhost();
    flux->setZeros();
    for(ind i=0-nGhost[LEFTT];i<nHalf+nGhost[RIGHT];i++)
    {
        (this->*calTypeFlux)(i);
    }
}



void SpaceDis::setMethod(EquationType type_,DiffMethod method_)
{
    fluxType=type_;
    diffMethod=method_;
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

    switch (diffMethod)
    {
    case TRAD2:
        difMethod=&SpaceDis::dif2Order;
        break;
    case TRAD6:
        difMethod=&SpaceDis::difTraditional6;
        break;
    case HDS6:
        difMethod=&SpaceDis::difHCS;
        break;
    
    default:
        break;
    }
}


real SpaceDis::at(int i,int ivar)
{
    //get the values in data
    if (i<0) 
    {
        return (*bndL)(-(i+1),ivar);
    }
    if (i>=n) 
    {
        return (*bndR)(i-n,ivar);
    }
    return (*data)[(i0+offset*i)*nVar+ivar];
}

real& SpaceDis::fAt(int i,int ivar)
{
    if (i<0) 
    {
        return (*fBndL)(-(i+1),ivar);
    }
    if (i>=n) 
    {
        return (*fBndR)(i-n,ivar);
    }
    return (*flux)[i*nVar+ivar];
}
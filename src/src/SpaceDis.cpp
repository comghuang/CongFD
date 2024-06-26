#include "SpaceDis.hpp"


SpaceDis::SpaceDis(){};
SpaceDis::SpaceDis(int n_,Data* data_,Data* rhs_
            ,std::shared_ptr<OneDBnd> bndL_,std::shared_ptr<OneDBnd> bndR_,Info* info_)
{
    data=data_;
    info=info_;
    nPrim=info->nPrim();
    n=n_;
    nVar=info->nCons();
    nHalf=n+1;
    bndL=bndL_;
    bndR=bndR_;
    flux=std::make_shared<Data>(nHalf,nVar);

    fBndL=std::make_shared<OneDBnd>(info->nFluxPoint(),nVar,FLUXGHOST);
    fBndR=std::make_shared<OneDBnd>(info->nFluxPoint(),nVar,FLUXGHOST);
    rhs=rhs_;

    fluxType=info->eqType;
    diffMethod=info->diffMethod;
    interMethod=info->interMethod;
    switch (fluxType)
    {
    case LINEARCONV1D:
        calTypeFlux=&SpaceDis::calFluxConv;
        break;
    case BURGERS1D:
        calTypeFlux=&SpaceDis::calFluxBurgers;
        break;
    case EULER:
        if (info->dim()==1)
        calTypeFlux=&SpaceDis::calFluxEuler1D;
        else if (info->dim()==2)calTypeFlux=&SpaceDis::calFluxEuler2D;
        
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



void SpaceDis::difference()
{
    calFlux();
    (this->*difMethod)();
}

void SpaceDis::calFlux()
{
    int fGhostL=fBndL->getN(),fGhostR=fBndR->getN();
    flux->setZeros();
    for(int i=0-fGhostL;i<nHalf+fGhostR;i++)
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
    case EULER:
        calTypeFlux=&SpaceDis::calFluxEuler1D;
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
    return (*data)[(i0+offset*i)*nPrim+ivar];
}

real& SpaceDis::fluxAt(int i,int ivar)
{
    if (i<0) 
    {
        return (*fBndL)(-(i+1),ivar);
    }
    if (i>=nHalf) 
    {
        return (*fBndR)(i-nHalf,ivar);
    }
    return (*flux)[i*nVar+ivar];
}

void SpaceDis::setConstNorm(std::array<real,3>&& norm_)
{
    norm=norm_;
}

void SpaceDis::setIDim(int idim_)
{
    //x:0,y:1,z:2 要直接用来作索引
    idim=idim_;
}
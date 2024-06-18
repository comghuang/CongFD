#include "zone.hpp"
static std::map<FluxType,std::string> fluxStr= {{LINEARCONV1D,"LINEARCONV1D"}
                                        ,{BURGERS1D,"BURGERS1D"}
                                        ,{EULER1D,"EULER1D"}};
static std::map<SpaceDisMethod,std::string> disStr={
    {FIRSTORDER,"FIRSTORDER"},
    {MUSCL,"MUSCL"},
    {WCNSJS5,"WCNSJS5"}
};
void Zone::init(Info info_,std::shared_ptr<Block> grid_)
{
    
    iMax=grid_->icMax;
    dim=grid_->dim;
    fluxType=info_.fType;
    spMethod=info_.spMethod;
    grid=grid_;
    switch (fluxType)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        nVar=1;
        nPrim=1;
        break;

    case EULER1D:
        nVar=3;
        nPrim=5;
        break;
    
    default:
        break;
    }
    len=iMax[0]*iMax[1]*iMax[2];
    
    
    cons->init(len,nVar);
    
    prim->init(len,nPrim);
    rhs.init(len,nVar);
    discrete.init(prim,grid,len,nVar,nPrim);
    discrete.setMethod(spMethod,fluxType);

    cons->solInit(len,nVar);

    grid->outputCgns();
    cons->oneDsolOutput(0,fluxStr[fluxType]+disStr[spMethod]);
    
}

void Zone::RK3(real dt)
{
        Data tempdata=(*cons);

        //third order RK
        //stage 1
        consToPrim();
        rhs.setZeros();
        rhs+=discrete.difference();
        for(ind i=0;i<len;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            (*cons)[i*nVar+ivar]=tempdata[i*nVar+ivar]-dt*rhs[i*nVar+ivar];
        }

        //stage 2
        consToPrim();
        rhs.setZeros();
        rhs+=discrete.difference();
        for(ind i=0;i<len;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            (*cons)[i*nVar+ivar]=0.75*tempdata[i*nVar+ivar]
                             -0.25*dt*rhs[i*nVar+ivar]
                             +0.25*(*cons)[i*nVar+ivar];
        }

        //stage 3
        consToPrim();
        rhs.setZeros();
        rhs+=discrete.difference();
        for(ind i=0;i<len;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            (*cons)[i*nVar+ivar]=1.0/3.0*tempdata[i*nVar+ivar]
                             -2.0/3.0*dt*rhs[i*nVar+ivar]
                             +2.0/3.0*(*cons)[i*nVar+ivar];
        }
}


void Zone::oneDsolOutput(real t)
{
    cons->oneDsolOutput(t,fluxStr[fluxType]+disStr[WCNSJS5]);
}

void Zone::consToPrim()
{
    switch (fluxType)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        std::copy(cons->begin(),cons->end(),prim->begin());
        break;
    case EULER1D:
        consToPrimEuler1D();
        break;
    
    default:
        break;
    }
}

void Zone::consToPrimEuler1D()
{
    if(iMax[1]!=1&&iMax[2]!=1)
    {
        std::cout<<"Euler 1d equation dimension error \n";
        return;
    }
    if(nVar!=3,nPrim!=5)
    {
        std::cout<<"Euler 1d equation variable number error \n";
    }

    for (ind i = 0; i < iMax[0]; i++)
    {
        real r=(*cons)(i,0);
        real ru=(*cons)(i,1);
        real rE=(*cons)(i,2);
        real u=ru/r;
        real E=rE/r;
        real e=E-u*u/2;
        real gamma=GAMMA;
        real RT=(gamma-1)*e;
        real p=r*RT;
        real H=gamma/(gamma-1)*RT+u*u/2;
        (*prim)(i,0)=r;
        (*prim)(i,1)=u;
        (*prim)(i,2)=p;
        (*prim)(i,3)=H;
        (*prim)(i,4)=RT;
    }
    
}

Zone::Zone()
{
    prim=std::make_shared<Data>();
    cons=std::make_shared<Data>();
}

std::shared_ptr<Data> Zone::getCons()
{
    return cons;
}
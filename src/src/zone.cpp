#include "zone.hpp"
std::map<FluxType,std::string> fluxStr= {{LINEARCONV1D,"LINEARCONV1D"}
                                        ,{BURGERS1D,"BURGERS1D"}
                                        ,{EULER1D,"EULER1D"}};
std::map<SpaceDisMethod,std::string> disStr={
    {FIRSTORDER,"FIRSTORDER"},
    {MUSCL,"MUSCL"},
    {WCNSJS5,"WCNSJS5"}
};
void Zone::init(ind I_,ind J_,ind K_,ind dim_,FluxType type)
{
    
    
    iMax[0]=I_;
    iMax[1]=J_;
    iMax[2]=K_;
    dim=dim_;
    switch (type)
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
    fluxType=type;
    len=I_*J_*K_;
    dat.init(len,nVar);
    prim.init(len,nPrim);
    coor.init(len,dim);
    rhs.init(len,nVar);
    discrete.init(&prim,&coor,len,nVar,nPrim);
    discrete.setMethod(WCNSJS5,type);


    coor.uniMesh();
    dat.solInit(len,nVar);

    coor.cgnsoutputInit1D();
    dat.oneDsolOutput(0,fluxStr[type]+disStr[WCNSJS5]);
    
}

void Zone::RK3(real dt)
{
        Data tempdata=dat;

        //third order RK
        //stage 1
        consToPrim();
        rhs.setZeros();
        rhs+=discrete.difference();
        for(ind i=0;i<len;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            dat[i*nVar+ivar]=tempdata[i*nVar+ivar]-dt*rhs[i*nVar+ivar];
        }

        //stage 2
        consToPrim();
        rhs.setZeros();
        rhs+=discrete.difference();
        for(ind i=0;i<len;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            dat[i*nVar+ivar]=0.75*tempdata[i*nVar+ivar]
                             -0.25*dt*rhs[i*nVar+ivar]
                             +0.25*dat[i*nVar+ivar];
        }

        //stage 3
        consToPrim();
        rhs.setZeros();
        rhs+=discrete.difference();
        for(ind i=0;i<len;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            dat[i*nVar+ivar]=1.0/3.0*tempdata[i*nVar+ivar]
                             -2.0/3.0*dt*rhs[i*nVar+ivar]
                             +2.0/3.0*dat[i*nVar+ivar];
        }
}


void Zone::oneDsolOutput(real t)
{
    dat.oneDsolOutput(t,fluxStr[fluxType]+disStr[WCNSJS5]);
}

void Zone::consToPrim()
{
    switch (fluxType)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        prim=dat;
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
        real r=dat(i,0);
        real ru=dat(i,1);
        real rE=dat(i,2);
        real u=ru/r;
        real E=rE/r;
        real e=E-u*u/2;
        real gamma=GAMMA;
        real RT=(gamma-1)*e;
        real p=r*RT;
        real H=gamma/(gamma-1)*RT+u*u/2;
        prim(i,0)=r;
        prim(i,1)=u;
        prim(i,2)=p;
        prim(i,3)=H;
        prim(i,4)=RT;
    }
    
}
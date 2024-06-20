#include "equation.hpp"

void Equation::consToPrim()
{
    switch (type)
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

void Equation::consToPrimEuler1D()
{
    if(nCons!=3,nPrim!=5)
    {
        std::cout<<"Equation error: Euler 1d equation variable number error \n";
    }

    for (ind i = 0; i < n; i++)
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

std::shared_ptr<Data> Equation::getPrim()
{
    return prim;
}

std::shared_ptr<Data> Equation::getCons()
{
    return cons;
}

std::shared_ptr<Data> Equation::getRhs()
{
    return rhs;
}
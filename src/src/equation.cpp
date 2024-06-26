#include "equation.hpp"

void Equation::consToPrim()
{
    switch (type)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        std::copy(cons->begin(),cons->end(),prim->begin());
        break;
    case EULER:
        {
            if (dim==1) consToPrimEuler1D();
            else if(dim==2) consToPrimEuler2D();
        }
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

    for (int i = 0; i < n; i++)
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
void Equation::consToPrimEuler2D()
{
    if(nCons!=3,nPrim!=5)
    {
        std::cout<<"Equation error: Euler 1d equation variable number error \n";
    }

    for (int i = 0; i < n; i++)
    {
        real r=(*cons)(i,0);
        real ru=(*cons)(i,1);
        real rv=(*cons)(i,2);
        real rE=(*cons)(i,3);
        real u=ru/r;
        real v=rv/r;
        real E=rE/r;
        real q2=(u*u+v*v)/2;
        real e=E-q2;
        real gamma=GAMMA;
        real RT=(gamma-1)*e;
        real p=r*RT;
        real H=gamma/(gamma-1)*RT+q2;
        (*prim)(i,0)=r;
        (*prim)(i,1)=u;
        (*prim)(i,2)=v;
        (*prim)(i,3)=p;
        (*prim)(i,4)=H;
    }
}

Data* Equation::getPrim()
{
    return prim;
}

Data* Equation::getCons()
{
    return cons;
}

Data* Equation::getRhs()
{
    return rhs;
}
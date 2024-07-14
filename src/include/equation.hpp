#pragma once

#include"Data.hpp"

class Equation
{
    public:
    Equation(){};
    void consToPrim();
    
    Data* getCons();
    Data* getPrim();
    Data* getRhs();

    protected:
    friend class Initializer;

    void consToPrimEuler1D();
    void consToPrimEuler2D();
    void consToPrimEuler1DHLL();//没啥用的东西
    void consToPrimEuler2DHLL();
    bool inited=false;
    int n,nCons,nPrim;
    int dim;
    Data* cons,*rhs,*prim;
    EquationType type;
};
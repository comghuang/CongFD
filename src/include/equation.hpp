#pragma once

#include"Data.hpp"

class Equation
{
    public:
    Equation(){};
    void consToPrim();
    
    std::shared_ptr<Data> getCons();
    std::shared_ptr<Data> getPrim();
    std::shared_ptr<Data> getRhs();

    protected:
    friend class Initializer;

    void consToPrimEuler1D();
    bool inited=false;
    int n,nCons,nPrim;
    std::shared_ptr<Data> cons,rhs,prim;
    EquationType type;
};
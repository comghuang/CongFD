#pragma once
#include "initializer.hpp"
#include "cgnsio.hpp"

class BlockSolver
{
    public:
    BlockSolver();
    BlockSolver(Info);
    void solve();
    ~BlockSolver();


    void outputGrid();
    void outputCons();
    void outputPrim();


    void stepsLoop();
    void Test();

    private:
    
    CgnsIO cgnsIO;
    Info* info;
    Block* block;
    Initializer* initer;
    Equation* eqn;
    Bnds* bnds;
    SpDistributor* spDis;
    Data* cons,*rhs;
    void RK3();
    TimeMethod timeMethod=RK3SSP;
};
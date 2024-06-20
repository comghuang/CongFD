#pragma once
#include "initializer.hpp"
class BlockSolver
{
    public:
    BlockSolver();
    BlockSolver(std::shared_ptr<Info>);
    void solve();
    private:
    std::shared_ptr<Info> info;
    std::shared_ptr<Block> block;
    std::shared_ptr<Initializer> initer;
    std::shared_ptr<Equation> eqn;
    std::shared_ptr<Bnds> bnds;
    std::shared_ptr<SpDistributor> spDis;
    std::shared_ptr<Data> cons,rhs;
    void RK3();
    TimeMethod timeMethod=RK3SSP;
    double dt=0.001;
};
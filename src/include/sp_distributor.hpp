#pragma once
#include "Bnds.hpp"
#include "SpaceDis.hpp"

class SpDistributor
{
    public:
    void rhsSolve();

    private:
    friend class Initializer;

    std::array<int,2> calOffset(int dim,int i,int j);
    int nCons,nPrim,dim;
    std::array<int,3> iMax;
    std::shared_ptr<Data> prim,cons;
    std::shared_ptr<Data> rhs;
    std::shared_ptr<Bnds> bnds;
    std::shared_ptr<Info> info;
};
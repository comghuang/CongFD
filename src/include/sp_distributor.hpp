#pragma once
#include "Bnds.hpp"
#include "SpaceDis.hpp"

class SpDistributor
{
    public:
    void consToPrim(std::shared_ptr<Data>,std::shared_ptr<Data>);


    private:
    int n,nCons,nPrim;
    std::array<ind,3> iMax;
    std::shared_ptr<Data> prim,cons;
    std::shared_ptr<Data> rhs;
    std::shared_ptr<Block> grid;
    std::shared_ptr<Bnds> bnds;
};
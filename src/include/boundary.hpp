#pragma once
#include "macro.hpp"
#include <vector>
#include <array>


class OneDBnd
{
    public:
    real& operator()(ind,ind);
    void init(ind,ind,BndType);
    void setGhostValue(real*);
    ind getN();
    BndType getType();

    private:
    BndType type=TYPENULL;
    std::vector<real> data;
    ind n=0;
    ind nVar=0;
};
enum bndPosition{
    IL,
    IR,
    JL,
    JR,
    KL,
    KR
};

class bndCache
{
    /*Cache for Periodical Boundary*/
    
    private:
    std::vector<real> data;


};

class Bnds
{
    public:
    void init(ind dim_,std::array<ind,3> iMax_);

    OneDBnd* getBnd(bndPosition,std::array<ind,2>);
    void setBnd(bndPosition,std::array<ind,4>);
    private:
    ind dim;
    std::array<ind,3> indexx;
    std::vector<OneDBnd> bnds;
};
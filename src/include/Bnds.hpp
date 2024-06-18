#pragma once
#include "block.hpp"



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
    void initFromCode(std::shared_ptr<Block>);

    OneDBnd* getBnd(bndPosition,std::array<ind,2>);
    void setBnd(bndPosition,std::array<ind,4>);
    private:
    ind dim;
    std::array<ind,3> iMax;
    //idim*
    std::vector<std::shared_ptr<OneDBnd>> bnds;
};
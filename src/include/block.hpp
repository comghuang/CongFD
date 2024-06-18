#pragma once
#include"data.hpp"

class Block
{
    public:
    void initUniform(std::array<int,3>,std::array<double,6>);
    void outputCgns();
    real operator()(int,int);

    Data coorVer;
    Data coorCel;
    int nVer,nCel;
    
    int dim;
    std::array<int,3> iMax,icMax;

};
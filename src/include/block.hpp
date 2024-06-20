#pragma once
#include"data.hpp"

class Block
{
    public:
    void outputCgns();
    real operator()(int,int);


    protected:
    friend class Initializer;
    Data coorVer;
    Data coorCel;
    int nVer,nCel;
    bool inited=false;
    int dim;
    std::array<int,3> iMax,icMax;

};
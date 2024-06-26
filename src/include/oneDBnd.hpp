#pragma once

#include"macro.hpp"
#include "data.hpp"
class OneDBnd
{
    public:
    OneDBnd(){};
    OneDBnd(int,int,BndType);
    real& operator()(int,int);

    int getN();
    BndType getType();
    void setUpdate(Data*,int,int);

    void update();

    private:
    void setValue(std::vector<real>);
    BndType type=TYPENULL;
    std::vector<real> data;
    Data* prim;
    int i0=0,offset=1;
    int n=0;
    int nVar=0;
};
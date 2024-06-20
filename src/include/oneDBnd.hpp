#pragma once

#include"macro.hpp"
class OneDBnd
{
    public:
    OneDBnd(){};
    OneDBnd(ind,ind,BndType);
    real& operator()(ind,ind);

    void setValue(std::vector<real>);
    ind getN();
    BndType getType();

    private:
    BndType type=TYPENULL;
    std::vector<real> data;
    ind n=0;
    ind nVar=0;
};
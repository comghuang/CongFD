#pragma once

#include"macro.hpp"
class OneDBnd
{
    public:
    real& operator()(ind,ind);
    void init(ind,ind,BndType);
    void setValue(real*);
    ind getN();
    BndType getType();

    private:
    BndType type=TYPENULL;
    std::vector<real> data;
    ind n=0;
    ind nVar=0;
};
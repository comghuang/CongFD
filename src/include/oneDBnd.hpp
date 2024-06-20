#pragma once

#include"macro.hpp"
#include "data.hpp"
class OneDBnd
{
    public:
    OneDBnd(){};
    OneDBnd(ind,ind,BndType);
    real& operator()(ind,ind);

    ind getN();
    BndType getType();
    void setUpdate(std::shared_ptr<Data>,int,int);

    void update();

    private:
    void setValue(std::vector<real>);
    BndType type=TYPENULL;
    std::vector<real> data;
    std::shared_ptr<Data> prim;
    int i0=0,offset=1;
    ind n=0;
    ind nVar=0;
};
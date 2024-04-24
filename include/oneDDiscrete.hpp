
#pragma once

#include <vector>
#include "vecs.hpp"
#include "macro.hpp"
class oneDDiscrete
{
    public:
    oneDDiscrete(std::vector<real>*,ind,ind);
    std::vector<real> difference();
    void calFlux();
    void calFluxConv();
    real reconL(ind,ind);
    real reconR(ind,ind);


    private:
    vecs data;
    ind n,nVar;
    ind nHalf;
    vecs flux;

};
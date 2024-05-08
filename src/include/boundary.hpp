#pragma once
#include "macro.hpp"
#include <vector>


class OneDGhost
{
    public:
    real& operator()(ind,ind);
    void init(ind,ind);
    void setGhostValue(std::vector<real>);
    ind getN();

    private:
    std::vector<real> data;
    ind n=5;
    ind nVar;
};
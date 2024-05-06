
#pragma once
#include "macro.hpp"
#include "data.hpp"
class SpaceDis
{
    public:
    SpaceDis(Data*,Data*,ind,ind);
    std::vector<real> difference();
    void calFlux();
    void calFluxConv();
    void calFluxBurgers();
    real reconL(ind,ind);
    real reconR(ind,ind);


    private:
    Data* data,*coor;
    ind n,nVar;
    ind nHalf;
    Data flux;

};
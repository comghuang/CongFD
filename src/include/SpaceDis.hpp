
#pragma once
#include "macro.hpp"
#include "data.hpp"
#include <algorithm>



class SpaceDis
{
    public:
    SpaceDis(Data*,Data*,ind,ind);
    SpaceDis();
    std::vector<real> difference();

    void setMethod(SpaceDisMethod,FluxType);
    void init(Data*,Data*,ind,ind,ind);

    void calFlux();
    void (SpaceDis::*calTypeFlux)(ind);
    void calFluxConv(ind);
    void calFluxBurgers(ind);
    void calFluxEuler(ind);

    real reconL(ind,ind);
    real reconR(ind,ind);

    std::vector<real> (SpaceDis::*difMethod)();
    std::vector<real> difHCS();
    std::vector<real> difTraditional6();
    std::vector<real> dif2Order();
    


    private:
    Data* data,*coor;
    ind n,nVar,nPrim;
    ind nHalf;
    Data flux;
    FluxType fluxType;
    SpaceDisMethod spDisMethod;


};
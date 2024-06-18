
#pragma once
#include "data.hpp"
#include "block.hpp"




class SpaceDis
{
    public:
    SpaceDis();
    std::vector<real> difference();

    void setMethod(SpaceDisMethod,FluxType);
    void init(std::shared_ptr<Data>,std::shared_ptr<Block>,ind,ind,ind);

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
    std::shared_ptr<Data> data;
    std::shared_ptr<Block> grid;
    ind n,nVar,nPrim;
    ind nHalf;
    Data flux;
    FluxType fluxType;
    SpaceDisMethod spDisMethod;


};
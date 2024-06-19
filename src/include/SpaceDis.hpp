
#pragma once
#include "data.hpp"
#include "block.hpp"




class SpaceDis
{
    public:
    SpaceDis();
    std::vector<real> difference();

    void setMethod(SpaceDisMethod,EquationType);
    void init(std::shared_ptr<Data>,std::shared_ptr<Block>,ind,ind,ind);

    void calFlux();
    void (SpaceDis::*calTypeFlux)(ind);
    void calFluxConv(ind);
    void calFluxBurgers(ind);
    void calFluxEuler(ind);

    real at(int,int);
    real reconL(ind,ind);
    real reconR(ind,ind);

    std::vector<real> (SpaceDis::*difMethod)();
    std::vector<real> difHCS();
    std::vector<real> difTraditional6();
    std::vector<real> dif2Order();
    


    private:
    std::shared_ptr<Data> data;
    std::shared_ptr<Block> grid;
    std::shared_ptr<Data> rhs;
    std::shared_ptr<OneDBnd> bndL,bndR;
    
    ind n,nVar,nPrim;
    int i0=0,offset=1;
    ind nHalf;
    Data flux;
    OneDBnd fBndL,fBndR;
    EquationType fluxType;
    SpaceDisMethod spDisMethod;


};
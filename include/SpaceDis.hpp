
#pragma once
#include "macro.hpp"
#include "data.hpp"

enum SpaceDisMethod{
    FIRSTORDER,
    MUSCL,
    WCNSJS5
};

enum FluxType{
    LINEARCONV,
    BURGERS,
    EULER
};

class SpaceDis
{
    public:
    SpaceDis(Data*,Data*,ind,ind);
    std::vector<real> difference();

    void setMethod(SpaceDisMethod,FluxType);

    void calFlux();
    void calFluxConv(ind);
    void calFluxBurgers(ind);
    void calFluxEuler(ind);

    std::array<real,2> getLR(ind,ind);
    real reconL(ind,ind);
    real reconR(ind,ind);
    


    private:
    Data* data,*coor;
    ind n,nVar;
    ind nHalf;
    Data flux;
    FluxType fluxType;
    SpaceDisMethod spDisMethod;


};
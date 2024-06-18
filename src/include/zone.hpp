#pragma once
#include "data.hpp"
#include "SpaceDis.hpp"

#include "info.hpp"


class Zone
{
    public:
    Zone();
    void init(Info,std::shared_ptr<Block>);
    void RK3(real);
    void oneDsolOutput(real);
    void consToPrim();
    void consToPrimEuler1D();
    std::shared_ptr<Data> getCons();


    private:
    FluxType fluxType;
    SpaceDisMethod spMethod;
    ind nVar,nPrim,dim,len;
    std::array<ind,3> iMax;
    ind offset,i0;
    std::shared_ptr<Data> prim,cons;
    SpaceDis discrete;
    Data rhs;
    std::shared_ptr<Block> grid;
    std::vector<OneDBnd> bnds;
};
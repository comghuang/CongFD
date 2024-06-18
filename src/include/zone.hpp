#pragma once
#include "data.hpp"
#include "SpaceDis.hpp"
#include <map>

class Zone
{
    public:
    void init(ind,ind,ind,ind,FluxType);
    void RK3(real);
    void oneDsolOutput(real);
    void consToPrim();
    void consToPrimEuler1D();


    private:
    FluxType fluxType;
    ind nVar,nPrim,dim,len;
    std::array<ind,3> iMax;
    ind offset,i0;
    Data dat;
    Data prim;
    SpaceDis discrete;
    Data rhs;
    Data coor;
    std::vector<OneDBnd> bnds;
};
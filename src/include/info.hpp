#pragma once
#include "macro.hpp"

class Info
{
    public:
    SpaceDisMethod spMethod=WCNSJS5;
    EquationType eqType=EULER1D;
    int nCase=0;
    DiffMethod diffMethod=TRAD6;

    std::array<int,3> iMax{100,2,2};
    std::array<double,6> calZone{-1,1,0,2,0,2};

    int nGhostCell();
    int nFluxPoint();
    int nPrim();
    int nCons();
    int dim();
    BndType defaultBndType();
    std::array<int,3> icMax();

};
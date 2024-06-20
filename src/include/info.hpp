#pragma once
#include "macro.hpp"

class Info
{
    public:
    InterMethod spMethod=WCNSJS5;
    EquationType eqType=EULER1D;
    int nCase=0;
    real t=0;
    DiffMethod diffMethod=TRAD6;
    InterMethod interMethod=WCNSJS5;

    std::array<int,3> iMax{201,2,2};
    std::array<double,6> calZone{-1,1,0,2,0,2};

    int nGhostCell();
    int nFluxPoint();
    int nPrim();
    int nCons();
    int dim();
    BndType defaultBndType();
    std::array<int,3> icMax();
    std::string filename();

};
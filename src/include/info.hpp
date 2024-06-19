#pragma once
#include "macro.hpp"

class Info
{
    public:
    SpaceDisMethod spMethod=WCNSJS5;
    EquationType eqType=EULER1D;
    int nCase=0;
    DiffMethod diffMethod=TRAD6;

    int nGhostCell();
    int nFluxPoint();
    int nPrim();
    int nCons();
};
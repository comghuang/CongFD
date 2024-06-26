#pragma once
#include "macro.hpp"


class Info
{
    public:
    InterMethod spMethod=FIRSTORDER;
    EquationType eqType=EULER;
    int nCase=0;
    real t=0;
    real dt=0.0001;
    int step=0;
    int endStep=10000;
    int outputInterval=100;
    DiffMethod diffMethod=TRAD2;
    InterMethod interMethod=FIRSTORDER;

    std::array<int,3> iMax{501,501,2};
    std::array<double,6> calZone{-0.5,0.5,-0.5,0.5,0,2};

    int nGhostCell();
    int nFluxPoint();
    int nPrim();
    int nCons();
    int dim();
    real geth(int);
    BndType defaultBndType();
    std::array<int,3> icMax();
    std::string filename();

    std::vector<std::string> getVarNameListCons();
    std::vector<std::string> getVarNameListPrim();
    std::vector<std::string> getVarNameListRhs();
};
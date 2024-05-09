#pragma once

#include "macro.hpp"
#include <fstream>
#include <cgnslib.h>
#include "boundary.hpp"
#include <vector>
#include <iostream>
#include <string>
#include <cstring>
#include <array>

class Data
{
    public:
    void solInit(ind,ind);
    void init(ind,ind);
    void setValue(real*,ind);
    //void setDim(ind,std::vector<ind>);
    real& operator() (ind,ind);
    real& operator[] (ind);
    
    real maxElement(ind);

    void uniMesh();
    void output(std::fstream*,ind);


    void cgnsoutputInit();
    void oneDsolOutput(real);

    //for ghost vertex
    void setGhostVertex(ind,ind);
    void updateGhostVertex();
    std::array<ind,2> getNGhost();

    private:
    std::vector<real> data;
    std::vector<real> x;
    ind n=200;
    ind nVar=1;
    OneDGhost ghVertex[2];
};
#pragma once

#include "macro.hpp"
#include <fstream>
#include <cgnslib.h>

class Data
{
    public:
    void init(ind,ind);
    //void setDim(ind,std::vector<ind>);
    real& operator() (ind,ind);
    real& operator[] (ind);
    
    real maxElement(ind);

    void uniMesh();
    void output(std::fstream*,ind);


    void cgnsoutputInit();
    void oneDsolOutput(real);

    private:
    std::vector<real> data;
    std::vector<real> x;
    ind n=200;
    ind nVar=1;
};
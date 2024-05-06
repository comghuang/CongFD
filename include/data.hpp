#pragma once

#include "macro.hpp"
#include <fstream>

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

    private:
    std::vector<real> data;
    std::vector<real> x;
    ind n=200;
    ind nVar=1;
};
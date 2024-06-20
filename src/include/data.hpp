#pragma once
#include "macro.hpp"

class Data
{
    public:

    Data(){};
    Data(int,int);


    void solInit(ind,ind);
    void init(ind,ind);
    void setValue(std::vector<real>);
    //void setDim(ind,std::vector<ind>);
    real& operator() (ind,ind);
    real& operator[] (ind);
    
    void operator= (Data&);
    void operator+= (std::vector<real>);
    void setZeros();
    real size();
    std::vector<real>::iterator begin();
    std::vector<real>::iterator end();
    
    //for global LF flux in burgers equation
    real maxElement(ind);

    private:
    std::vector<real> data;
    ind n=200;
    ind nVar=1;
};
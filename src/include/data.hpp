#pragma once
#include <fstream>
#include <cgnslib.h>
#include "boundary.hpp"
#include <iostream>
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
    void operator= (Data&);
    void operator+= (std::vector<real>);
    void setZeros();
    
    //for global LF flux in burgers equation
    real maxElement(ind);

    //get uniform mesh in 1D workbench
    void uniMesh();
    //output sol to a file;
    void output(std::fstream*,ind);

    //output to cgns file 1D 
    void cgnsoutputInit1D();
    void cgnsoutputInit2D();
    void oneDsolOutput(real,std::string);
    

    //for ghost vertex
    void setGhostVertex(OneDBnd*,OneDBnd*);
    void updateGhostVertex();
    std::array<ind,2> getNGhost();

    private:
    std::vector<real> data;
    std::vector<real> x;
    ind n=200;
    ind nVar=1;
    ind i0=0,offset=1;
    OneDBnd* ghVertex[2];
};
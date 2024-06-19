#pragma once
#include "oneDBnd.hpp"


class Data
{
    public:
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

    //get uniform mesh in 1D workbench
    void uniMesh();
    //output sol to a file;
    void output(std::fstream*,ind);

    //output to cgns file 1D 
    void cgnsoutputInit1D();
    void cgnsoutputInit2D();
    void oneDsolOutput(real,std::string);
    

    //for ghost vertex //准备废弃！不应该放在这里面
    void setGhostVertex(OneDBnd*,OneDBnd*);
    void updateGhostVertex();
    std::array<ind,2> getNGhost();

    private:
    std::vector<real> data;
    ind n=200;
    ind nVar=1;
    ind i0=0,offset=1;
    OneDBnd* ghVertex[2];
};
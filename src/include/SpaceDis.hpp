
#pragma once
#include "block.hpp"
#include "oneDBnd.hpp"
#include "info.hpp"




class SpaceDis
{
    public:
    SpaceDis(int n_,Data* data_,Data* rhs_
            ,std::shared_ptr<OneDBnd> bndL_,std::shared_ptr<OneDBnd> bndR_,Info* info);
    SpaceDis();
    void difference();
    void setOffset(int,int);
    void setMethod(EquationType ,DiffMethod);
    void setIDim(int);

    void setConstNorm(std::array<real,3>&&);
    


    private:
    Data* data;
    Data* rhs;
    Info* info;
    std::shared_ptr<OneDBnd> bndL,bndR;

    void calFlux();
    void calFluxConv(int);
    void calFluxBurgers(int);
    void calFluxEuler1D(int);
    void calFluxEuler2D(int);
    void (SpaceDis::*calTypeFlux)(int);

    real at(int,int);
    real& fluxAt(int,int);
    real reconL(int,int);
    real reconR(int,int);


    void difHCS();
    void difTraditional6();
    void dif2Order();
    void (SpaceDis::*difMethod)();
    
    int n,nVar,nPrim;
    int idim;
    int i0=0,offset=1;
    int nHalf;

    std::array<real,3> norm;

    std::shared_ptr<Data> flux;
    std::shared_ptr<OneDBnd> fBndL,fBndR;
    EquationType fluxType;
    DiffMethod diffMethod=TRAD2;
    InterMethod interMethod=FIRSTORDER;


};
#pragma once

#include "reconstructorFaceCenter.hpp"

class ReconerBVD : public Recon5OrderFaceCenter<nullptr> {
public:
    ReconerBVD() {};
    ReconerBVD(DataReader varsR_, std::shared_ptr<OneDBnd> bndL_,
        std::shared_ptr<OneDBnd> bndR_)
        : Recon5OrderFaceCenter<nullptr>(varsR_, bndL_, bndR_) {};
    void solve() final
    {
        initVarsR();
        iter = data->begin();
        candidateIter = candidateValues.begin();
        candidate2Iter = candidateValues2.begin();
        // candidate3Iter = candidateValues3.begin();
        leftBnd();
        internal();
        rightBnd();
        assert(iter == data->end());
        boundaryDimishing();
    }

private:
    Data candidateValues;
    std::vector<real>::iterator candidateIter;
    Data candidateValues2;
    std::vector<real>::iterator candidate2Iter;
    // Data candidateValues3;
    // std::vector<real>::iterator candidate3Iter;

    void initVarsR() final;
    void boundaryDimishing();
    void reconI(std::array<std::vector<real>::iterator, 6>) final;
};
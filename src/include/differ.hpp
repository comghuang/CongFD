#pragma once

#include "data.hpp"

#include "info.hpp"
#include <functional>

class Differ {
public:
    Differ() {};
    Differ(std::shared_ptr<Data> fluxesHalf_, std::shared_ptr<Data> fluxesNode_,
        DataReader rhsR_)
        : fluxesHalf(fluxesHalf_)
        , fluxesNode(fluxesNode_)
        , rhsR(rhsR_) {};
    void init(std::shared_ptr<Data> fluxesHalf_,
        std::shared_ptr<Data> fluxesNode_, DataReader rhsR_)
    {
        fluxesHalf = fluxesHalf_;
        fluxesNode = fluxesNode_;
        rhsR = rhsR_;
    }
    void setConstantH(real constantH_) { constantH = constantH_; }
    void check() { std::cout << "initialized successfully Differ\n"; }
    virtual void solve();

protected:
    std::shared_ptr<Data> fluxesHalf;
    std::shared_ptr<Data> fluxesNode;
    DataReader rhsR;
    real constantH;
    int idim;
};

class ExplicitDif6 : public Differ {
public:
    void solve() final;
};

class MidNodeAndNodeDif6 : public Differ {
public:
    void solve() final;
};
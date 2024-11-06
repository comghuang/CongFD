#pragma once
#include "dataManipulater.hpp"

class Reconstuctor : public DataManipulater {
public:
  Reconstuctor(DataReader varsR_, std::shared_ptr<OneDBnd> bndL_,
               std::shared_ptr<OneDBnd> bndR_, Info *info_)
      : DataManipulater(varsR_, bndL_, bndR_, info_) {};
  void solve() override {
    initVarsR();
    iter = data->begin();
    leftBnd();
    internal();
    rightBnd();
    assert(iter == data->end());
  }
};

class Recon1Order : public Reconstuctor {
public:
protected:
  void initVarsR() override;
  void leftBnd() override;
  void internal() override;
  void rightBnd() override;
};

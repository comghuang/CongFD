#pragma once
#include "data.hpp"
#include "macro.hpp"
#include "oneDBnd.hpp"
class DataManipulater {
public:
  DataManipulater(DataReader varsReader_, std::shared_ptr<OneDBnd> bndL_,
                  std::shared_ptr<OneDBnd> bndR_, Info *info_);
  void setConstNorm(std::array<real, 3> norm_);
  virtual void solve() = 0;

  std::shared_ptr<Data> data;

protected:
  std::vector<real>::iterator iter;
  virtual void initVarsR() = 0;
  virtual void leftBnd() = 0;
  virtual void internal() = 0;
  virtual void rightBnd() = 0;

  DataReader varsReader;
  std::shared_ptr<OneDBnd> bndL, bndR;
  std::array<real, 3> norm;
  Info *info;
  int n, nvar;
};
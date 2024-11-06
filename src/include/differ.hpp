#pragma once

#include "data.hpp"
#include "differenceScheme.hpp"
#include "info.hpp"
#include <functional>

class Differ {
public:
  Differ(std::shared_ptr<std::vector<real>> fluxesHalf_,
         std::shared_ptr<std::vector<real>> fluxesNode_, DataReader rhsR_,
         Info *info_)
      : fluxesHalf(fluxesHalf_), fluxesNode(fluxesNode_), rhsR(rhsR_),
        info(info_) {};
  virtual void solve();
  std::function<real(std::vector<real>::iterator, std::vector<real>::iterator,
                     int)>
      diffFunc = secondOrder;

protected:
  std::shared_ptr<std::vector<real>> fluxesHalf;
  std::shared_ptr<std::vector<real>> fluxesNode;
  DataReader rhsR;
  Info *info;
  int idim;
};

class ExplicitDif6 : public Differ {
public:
  void solve() override;
};
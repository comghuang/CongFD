#pragma once
#include "dataManipulater.hpp"
#include <concepts>

/*------------------concepts bigin---------------------------*/
template <typename T, std::size_t NVar>
concept Flux =
    requires(T f, std::span<real, NVar> arr, std::span<real, NVar> arr2,
             std::array<real, 3> norm) { T(arr, arr2, norm); };

template <std::size_t NVar>
using FluxN = void (*)(const std::span<real, NVar> &,
                       const std::span<real, NVar> &,
                       const std::array<real, 3> &);
/*-----------------concepts end--------------------------------*/

template <std::size_t NVar> class SolvePointSolver : public DataManipulater {
public:
  SolvePointSolver(DataReader varsReader_, std::shared_ptr<OneDBnd> bndL_,
                   std::shared_ptr<OneDBnd> bndR_, Info *info_,
                   FluxN<NVar> fluxSolver_, int nl_, int nr_)
      : DataManipulater(varsReader_, bndL_, bndR_, info_),
        fluxSolver(fluxSolver_), nl(nl_), nr(nr_) {};

  void solve() override {
    initVarsR();
    iter = data->begin();
    leftBnd();
    internal();
    rightBnd();
    assert(iter == data->end());
  };
  void initVarsR() override;
  void leftBnd() override;
  void internal() override;
  void rightBnd() override;
  void solveI(std::vector<real>::iterator);
  FluxN<NVar> fluxSolver;

protected:
  int nl = 0, nr = 0;
};
template <std::size_t NVar> void SolvePointSolver<NVar>::initVarsR() {
  if (!data) {
    data = std::make_shared<Data>(n + nl + nr, NVar);
  }
  assert(NVar == nvar);
}

template <std::size_t NVar> void SolvePointSolver<NVar>::leftBnd() {
  for (int i = nl - 1; i >= 0; i--) {
    auto iterBnd = (*bndL)(i);
    solveI(iterBnd);
  }
}

template <std::size_t NVar> void SolvePointSolver<NVar>::internal() {
  auto endIter = varsReader(n);
  for (auto iterVar = varsReader(0); iterVar != endIter; iterVar += NVar) {
    solveI(iterVar);
    assert(iterVar <= varsReader(n));
  }
}

template <std::size_t NVar> void SolvePointSolver<NVar>::rightBnd() {
  for (int i = 0; i < nr; i++) {
    auto iterBnd = (*bndR)(i);
    solveI(iterBnd);
  }
}

template <std::size_t NVar>
void SolvePointSolver<NVar>::solveI(std::vector<real>::iterator iterBnd) {
  std::span<real, NVar> input(iterBnd, NVar), output(iter, NVar);
  (*fluxSolver)(input, output, norm);
  iter += NVar;
}
#pragma once

#include "data.hpp"
#include "info.hpp"
#include <concepts>
#include <span>
#include <type_traits>
/*------------------concepts bigin---------------------------*/

template <std::size_t NVar>
using RiemannSolverN = void (*)(const std::span<real, NVar * 2> &,
                                const std::span<real, NVar> &,
                                const std::array<real, 3> &);

template <typename T, std::size_t NVar>
concept RiemannSolver =
    requires(T f, std::span<real, NVar> arr, std::array<real, 3> norm) {
      { T(arr, norm) } -> std::same_as<std::array<real, NVar>>;
    };

/*-----------------concepts end--------------------------------*/

template <std::size_t NVar, RiemannSolver<NVar> RSolver> class FluxPointSolver {
public:
  FluxPointSolver(std::shared_ptr<std::vector<real>> valsR_, Info *info_,
                  int idim_, int nvar_)
      : valsR(valsR_), info(info_), idim(idim_), nvar(nvar_) {};

  void setConstNorm(std::array<real, 3> norm_);
  virtual void solve() = 0;
  std::shared_ptr<Data> fluxes;

protected:
  std::shared_ptr<std::vector<real>> valsR;
  Info *info;
  int idim, nvar;
  std::array<real, 3> norm;
};

template <std::size_t NVar, RiemannSolver<NVar> RSolver>
void FluxPointSolver<NVar, RSolver>::setConstNorm(std::array<real, 3> norm_) {
  norm = norm_;
}

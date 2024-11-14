#pragma once

#include <concepts>
#include <span>
#include <type_traits>

#include "data.hpp"
#include "info.hpp"
/*------------------concepts bigin---------------------------*/

template <std::size_t NVar>
using RiemannSolverN = void (*)(const std::span<real, NVar * 2>&,
    const std::span<real, NVar>&,
    const std::array<real, 3>&);

template <typename T, std::size_t NVar>
concept RiemannSolver = requires(T f, std::span<real, NVar> arr, std::array<real, 3> norm) {
    { T(arr, norm) } -> std::same_as<std::array<real, NVar>>;
};

/*-----------------concepts end--------------------------------*/

template <std::size_t NVar, RiemannSolverN<NVar> solver>
class FluxPointSolver {
public:
    FluxPointSolver() {};
    FluxPointSolver(std::shared_ptr<Data> valsR_, int idim_, int nvar_)
        : valsR(valsR_)
        , idim(idim_)
        , nvar(nvar_) {};

    void init(std::shared_ptr<Data> valsR_, int idim_, int nvar_)
    {
        valsR = valsR_;
        idim = idim_;
        nvar = nvar_;
    }
    std::shared_ptr<Data> getData() { return fluxes; }
    void setConstNorm(std::array<real, 3> norm_);
    void check() { std::cout << "initialized successfully FluxPointSolver\n"; }
    void solve();
    std::shared_ptr<Data> fluxes;

protected:
    std::shared_ptr<Data> valsR;
    int idim, nvar;
    std::array<real, 3> norm;
};

template <std::size_t NVar, RiemannSolverN<NVar> solver>
void FluxPointSolver<NVar, solver>::setConstNorm(std::array<real, 3> norm_)
{
    norm = norm_;
}

template <std::size_t NVar, RiemannSolverN<NVar> solver>
void FluxPointSolver<NVar, solver>::solve()
{
    auto end = valsR->end();
    assert(nvar == NVar);
    assert(valsR->size() % 2 == 0);
    if (!fluxes)
        fluxes = std::make_shared<Data>(valsR->getN(), NVar);
    auto ivar = valsR->begin();
    auto iflux = fluxes->begin();
    for (; ivar != end; ivar += 2 * NVar, iflux += NVar) {
        std::span<real, 2 * NVar> input(ivar, 2 * NVar);
        std::span<real, NVar> output(iflux, NVar);
        solver(input, output, norm);
    }
    assert(ivar == end);
    assert(iflux == fluxes->end());
}

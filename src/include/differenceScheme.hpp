#pragma once
#include "macro.hpp"

class DifferenceScheme {
public:
    const int nFluxHalf = 0;
    const int nFluxNode = 0;
    virtual real operator()(std::vector<real>::iterator fluxHalf,
        std::vector<real>::iterator fluxNode, int nvar)
        = 0;
};

class DifSecondOrder : public DifferenceScheme {
public:
    const int nFluxHalf = 0;
    const int nFluxNode = 0;
    real operator()(std::vector<real>::iterator fluxHalf,
        std::vector<real>::iterator fluxNode, int nvar) override
    {
        return fluxHalf[1] - fluxHalf[0];
    };
};
real secondOrder(std::vector<real>::iterator fluxHalf,
    std::vector<real>::iterator fluxNode, int nvar)
{
    return fluxHalf[1] - fluxHalf[0];
}

real explicit6(std::vector<real>::iterator fluxHalf,
    std::vector<real>::iterator fluxNode, int nvar)
{
    constexpr std::array<real, 3> w = { 75.0 / 64.0, 25.0 / (128.0 * 3.0),
        3.0 / (128.0 * 50) };
    return w[0] * (*(fluxHalf - nvar * 3) - *(fluxHalf - nvar * 2)) - w[1] * (*(fluxHalf - nvar * 4) - *(fluxHalf - nvar * 1)) + w[2] * (*(fluxHalf - nvar * 5) - *(fluxHalf - nvar * 0));
}

// real MND6(std::vector<real>::iterator fluxHalf,
//           std::vector<real>::iterator fluxNode, int nvar) {
//   constexpr std::array<real, 3> w = {3.0 / 2.0, -3.0 / 10.0, 1.0 / 30.0};
//   return w[1] * (*(fluxNode - nvar * 1) - *(fluxNode - nvar * 0)) +
//          w[0] * (*(fluxHalf - nvar * 2) - *(fluxHalf - nvar * 1)) +
//          w[2] * (*(fluxHalf - nvar * 3) - *(fluxHalf - nvar * 0));
// }
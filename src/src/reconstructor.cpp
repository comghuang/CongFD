#include "reconstructor.hpp"
#include "eigenSystem.hpp"
#include "interScheme.hpp"

void Recon1Order::initVarsR()
{
    assert(bndL->getN() >= 1);
    assert(bndR->getN() >= 1);
    int nPointR = n + 1;
    data = std::make_shared<Data>(n * 2, nvar);
}
void Recon1Order::leftBnd()
{
    auto iterVar = varsReader(0);
    std::copy((*bndL)(0), (*bndL)(0) + nvar, iter);
    std::copy(iterVar, iterVar + nvar, iter + nvar);
    iter += 2 * nvar;
}
void Recon1Order::internal()
{
    for (int i = 0; i < n - 2; i++) {
        auto iterVar = varsReader(i);
        std::copy(iterVar, iterVar + nvar, iter);
        std::copy(iterVar + nvar, iterVar + 2 * nvar, iter + nvar);
    }
    iter += 2 * nvar;
}
void Recon1Order::rightBnd()
{
    auto iterVar = varsReader(n - 1);
    std::copy(iterVar, iterVar + nvar, iter);
    std::copy((*bndR)(0), (*bndR)(0) + nvar, iter + nvar);
    iter += 2 * nvar;
}
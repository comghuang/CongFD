#include "differ.hpp"
void Differ::solve() {
  assert(rhsR.getN() + 1 == fluxesHalf->size());
  int n = rhsR.getN(), nvar = rhsR.getNVar();
  real h = info->geth(rhsR.idim);

  for (auto iterRhs = rhsR(0); iterRhs < rhsR(n); iterRhs += rhsR.getOffset()) {
    auto iterFluxHalf = fluxesHalf->begin();
    auto iterFluxNode = fluxesNode->begin();
    for (int ivar = 0; ivar < nvar; ivar++) {
      iterRhs[ivar] =
          diffFunc(iterFluxHalf + ivar, iterFluxNode + ivar, nvar) / h;
    }
  }
}
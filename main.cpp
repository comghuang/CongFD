#include "solvePointFlux.hpp"

#include "fluxSchemes.hpp"
#include "oneDBnd.hpp"
#include "reconstructor5order.hpp"

int main() {
  auto bndl = std::make_shared<OneDBnd>(10, 3, SUPERSONICOUTLET);
  auto bndr = std::make_shared<OneDBnd>(10, 3, SUPERSONICOUTLET);
  auto data = new Data(50, 3);
  std::vector<real> initvalues(150);
  for (int i = 0; i < 150; i++) {
    initvalues[i] = i;
  }
  data->setValue(initvalues);

  DataReader dataReader(50, 0, 1, 0, data);
  // = std::make_shared<DataReader>(50, 0, 1, 0, data);
  Info *info = new Info;

  Recon5OrderFaceCenter Reconer(dataReader, bndl, bndr, info, weno5_JSchen);
  SolvePointSolver<3> solvePointFlux(dataReader, bndl, bndr, info,
                                     fluxSolveEuler1D, 1, 1);

  Reconer.solve();

  return 0;
}
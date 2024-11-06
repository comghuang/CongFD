#include "dataManipulater.hpp"
DataManipulater::DataManipulater(DataReader varsR_,
                                 std::shared_ptr<OneDBnd> bndL_,
                                 std::shared_ptr<OneDBnd> bndR_, Info *info_)
    : varsReader(varsR_), info(info_), bndL(bndL_), bndR(bndR_) {
  n = varsReader.getN();
  nvar = varsReader.getNVar();
}

void DataManipulater::setConstNorm(std::array<real, 3> norm_) { norm = norm_; }
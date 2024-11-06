#pragma once
#include "interScheme.hpp"
#include "reconstructor.hpp"
#include <functional>

// template <typename T>
// concept InterpolationScheme5Order = std::is_same()
// requires(T F, real a) {
//   { F(a, a, a, a, a) } -> std::same_as<real>;
// };

typedef real (*InterpolationScheme5Order)(std::array<real, 5>);

class Recon5OrderFaceCenter : public Reconstuctor {
public:
  Recon5OrderFaceCenter(DataReader varsR_, std::shared_ptr<OneDBnd> bndL_,
                        std::shared_ptr<OneDBnd> bndR_, Info *info_,
                        InterpolationScheme5Order inter5_)
      : Reconstuctor(varsR_, bndL_, bndR_, info_), inter5(inter5_) {};

protected:
  void initVarsR() override;
  void leftBnd() override;
  void internal() override;
  void rightBnd() override;
  InterpolationScheme5Order inter5 = weno5_JSchen;
  virtual void reconI(std::array<std::vector<real>::iterator, 6>);
  int tt = 0;
};

class Recon5Order1DEulerEig : public Recon5OrderFaceCenter {

public:
  Recon5Order1DEulerEig(DataReader varsR_, std::shared_ptr<OneDBnd> bndL_,
                        std::shared_ptr<OneDBnd> bndR_, Info *info_,
                        InterpolationScheme5Order inter5_)
      : Recon5OrderFaceCenter(varsR_, bndL_, bndR_, info_, inter5_) {};

protected:
  const int NVAR = 3;
  void reconI(std::array<std::vector<real>::iterator, 6>) final;
};

class Recon5Order2DEulerEig : public Recon5OrderFaceCenter {
public:
  Recon5Order2DEulerEig(DataReader varsR_, std::shared_ptr<OneDBnd> bndL_,
                        std::shared_ptr<OneDBnd> bndR_, Info *info_,
                        InterpolationScheme5Order inter5_)
      : Recon5OrderFaceCenter(varsR_, bndL_, bndR_, info_, inter5_) {};

protected:
  const int NVAR = 4;
  void reconI(std::array<std::vector<real>::iterator, 6>) final;
};
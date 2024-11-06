
#include "fluxScheme.hpp"
#include "macro.hpp"
#include <algorithm>

enum { L, R };

std::array<real, 3> roeFlux1D(real rl, real ul, real pl, real rr, real ur,
                              real pr) {
  std::array<real, 3> res;

  real gamma = GAMMA;
  arr2 H = {pl / rl * GAMMA / (GAMMA - 1) + (ul * ul) / 2,
            pr / rr * GAMMA / (GAMMA - 1) + (ur * ur) / 2};

  real rBar = std::sqrt(rl * rr);
  real uBar = (ul * std::sqrt(rl) + ur * std::sqrt(rr)) /
              (std::sqrt(rl) + std::sqrt(rr));
  real HBar = (H[L] * std::sqrt(rl) + H[R] * std::sqrt(rr)) /
              (std::sqrt(rl) + std::sqrt(rr));
  real cBar = std::sqrt((gamma - 1) * (HBar - uBar * uBar / 2));
  real cBar2 = (gamma - 1) * (HBar - uBar * uBar / 2);

  real dr = rr - rl;
  real du = ur - ul;
  real dp = pr - pl;

  real K1[3] = {1, uBar, uBar * uBar / 2};
  real K2[3] = {1, uBar - cBar, HBar - uBar * cBar};
  real K3[3] = {1, uBar + cBar, HBar + uBar * cBar};

  real alpha[3] = {dr - dp / cBar2,
                   1.0 / (2.0 * cBar2) * (dp - rBar * cBar * du),
                   1.0 / (2.0 * cBar2) * (dp + rBar * cBar * du)};
  real lambda[3] = {std::abs(uBar), std::abs(uBar - cBar),
                    std::abs(uBar + cBar)};
  /*
  real eps=0.1*(std::abs(uBar)+std::abs(uBar)+cBar);
  for(int i=0;i<3;i++)
  {
      if(lambda[i]<eps)
      {
          lambda[i]=(lambda[i]*lambda[i]+eps*eps)/(2.0*eps);
      }
  }*/

  real FL[3] = {rl * ul, rl * ul * ul + pl, rl * H[L] * ul};
  real FR[3] = {rr * ur, rr * ur * ur + pr, rr * H[R] * ur};

  for (int i = 0; i < 3; i++) {
    res[i] = (FL[i] + FR[i] - lambda[0] * alpha[0] * K1[i] -
              lambda[1] * alpha[1] * K2[i] - lambda[2] * alpha[2] * K3[i]) /
             2.0;
  }
  return res;
}

std::array<real, 3> roeFlux1D2(real rl, real ul, real pl, real rr, real ur,
                               real pr) {
  // reference: https://blog.csdn.net/Tankrun1997/article/details/132743487
  std::array<real, 3> res;

  real gamma = GAMMA;

  arr2 H = {pl / rl * GAMMA / (GAMMA - 1) + (ul * ul) / 2,
            pr / rr * GAMMA / (GAMMA - 1) + (ur * ur) / 2};

  real rBar = std::sqrt(rl * rr);
  real uBar = (ul * std::sqrt(rl) + ur * std::sqrt(rr)) /
              (std::sqrt(rl) + std::sqrt(rr));
  real HBar = (H[L] * std::sqrt(rl) + H[R] * std::sqrt(rr)) /
              (std::sqrt(rl) + std::sqrt(rr));
  real cBar = std::sqrt((gamma - 1) * (HBar - uBar * uBar / 2));
  real cBar2 = (gamma - 1) * (HBar - uBar * uBar / 2);

  real dr = rr - rl;
  real du = ur - ul;
  real dp = pr - pl;

  real K1[3] = {1, uBar, uBar * uBar / 2};
  real K2[3] = {1, uBar - cBar, HBar - uBar * cBar};
  real K3[3] = {1, uBar + cBar, HBar + uBar * cBar};

  real lambda[3] = {std::abs(uBar), std::abs(uBar + cBar),
                    std::abs(uBar - cBar)};

  real eps = 0.1 * (std::abs(uBar) + cBar);
  for (int i = 0; i < 3; i++) {
    if (lambda[i] < eps) {
      lambda[i] = (lambda[i] * lambda[i] + eps * eps) / (2.0 * eps);
    }
  }

  real alpha[5] = {lambda[0] * (dr - dp / cBar2),
                   lambda[1] / (2.0 * cBar2) * (dp + rBar * cBar * du),
                   lambda[2] / (2.0 * cBar2) * (dp - rBar * cBar * du)};
  alpha[3] = alpha[0] + alpha[1] + alpha[2];
  alpha[4] = cBar * (alpha[1] - alpha[2]);

  real FL[3] = {rl * ul, rl * ul * ul + pl, rl * H[L] * ul};
  real FR[3] = {rr * ur, rr * ur * ur + pr, rr * H[R] * ur};

  for (int i = 0; i < 3; i++) {
    res[i] = (FL[i] + FR[i]) / 2.0;
  }
  res[0] -= 0.5 * alpha[3];
  res[1] -= 0.5 * (alpha[3] * uBar + alpha[4]);
  res[2] -= 0.5 * (HBar * alpha[3] + uBar * alpha[4] -
                   cBar2 * alpha[0] / (gamma - 1));

  return res;
}

std::array<real, 3> roeFlux1D2(std::array<real, 6> prims) {
  // reference: https://blog.csdn.net/Tankrun1997/article/details/132743487

  real rl = prims[0];
  real ul = prims[1];
  real pl = prims[2];
  real rr = prims[3];
  real ur = prims[4];
  real pr = prims[5];
  std::array<real, 3> res;

  real gamma = GAMMA;

  arr2 H = {pl / rl * GAMMA / (GAMMA - 1) + (ul * ul) / 2,
            pr / rr * GAMMA / (GAMMA - 1) + (ur * ur) / 2};

  real rBar = std::sqrt(rl * rr);
  real uBar = (ul * std::sqrt(rl) + ur * std::sqrt(rr)) /
              (std::sqrt(rl) + std::sqrt(rr));
  real HBar = (H[L] * std::sqrt(rl) + H[R] * std::sqrt(rr)) /
              (std::sqrt(rl) + std::sqrt(rr));
  real cBar = std::sqrt((gamma - 1) * (HBar - uBar * uBar / 2));
  real cBar2 = (gamma - 1) * (HBar - uBar * uBar / 2);

  real dr = rr - rl;
  real du = ur - ul;
  real dp = pr - pl;

  real K1[3] = {1, uBar, uBar * uBar / 2};
  real K2[3] = {1, uBar - cBar, HBar - uBar * cBar};
  real K3[3] = {1, uBar + cBar, HBar + uBar * cBar};

  real lambda[3] = {std::abs(uBar), std::abs(uBar + cBar),
                    std::abs(uBar - cBar)};

  real eps = 0.1 * (std::abs(uBar) + cBar);
  for (int i = 0; i < 3; i++) {
    if (lambda[i] < eps) {
      lambda[i] = (lambda[i] * lambda[i] + eps * eps) / (2.0 * eps);
    }
  }

  real alpha[5] = {lambda[0] * (dr - dp / cBar2),
                   lambda[1] / (2.0 * cBar2) * (dp + rBar * cBar * du),
                   lambda[2] / (2.0 * cBar2) * (dp - rBar * cBar * du)};
  alpha[3] = alpha[0] + alpha[1] + alpha[2];
  alpha[4] = cBar * (alpha[1] - alpha[2]);

  real FL[3] = {rl * ul, rl * ul * ul + pl, rl * H[L] * ul};
  real FR[3] = {rr * ur, rr * ur * ur + pr, rr * H[R] * ur};

  for (int i = 0; i < 3; i++) {
    res[i] = (FL[i] + FR[i]) / 2.0;
  }
  res[0] -= 0.5 * alpha[3];
  res[1] -= 0.5 * (alpha[3] * uBar + alpha[4]);
  res[2] -= 0.5 * (HBar * alpha[3] + uBar * alpha[4] -
                   cBar2 * alpha[0] / (gamma - 1));

  return res;
}
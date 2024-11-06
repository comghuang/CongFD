#include "reconstructor5order.hpp"
#include "eigenSystem.hpp"

void Recon5OrderFaceCenter::initVarsR() {
  int nPointR = bndL->getN() + n + bndR->getN() - 5;
  if (!data) {
    data = std::make_shared<Data>(nPointR * 2, nvar);
  }
}

void Recon5OrderFaceCenter::leftBnd() {
  int nBnd = bndL->getN();
  std::vector<std::vector<real>::iterator> tempIters(nBnd + 5);
  auto tempItersIter = tempIters.begin();
  for (int i = 0; i < nBnd; i++) {
    int iInver = nBnd - 1 - i;
    auto tempppp = (*bndL)(iInver);
    (*tempItersIter++) = (*bndL)(iInver);
  }
  for (int i = 0; i < 5; i++) {
    (*tempItersIter++) = varsReader(i);
  }
  assert(tempItersIter == tempIters.end());

  for (int i = 0; i < nBnd; i++) {
    reconI({tempIters[i + 0], tempIters[i + 1], tempIters[i + 2],
            tempIters[i + 3], tempIters[i + 4], tempIters[i + 5]});
  }
}

void Recon5OrderFaceCenter::internal() {

  for (int i = 0; i < n - 5; i++) {
    auto iter = varsReader(i);
    reconI({
        varsReader(i + 0),
        varsReader(i + 1),
        varsReader(i + 2),
        varsReader(i + 3),
        varsReader(i + 4),
        varsReader(i + 5),
    });
  }
}

void Recon5OrderFaceCenter::rightBnd() {
  int nBnd = bndR->getN();
  std::vector<std::vector<real>::iterator> tempIters(nBnd + 5);
  auto tempItersIter = tempIters.begin();
  for (int i = 0; i < 5; i++) {
    (*tempItersIter++) = varsReader(n - 5 + i);
  }
  for (int i = 0; i < nBnd; i++) {
    (*tempItersIter++) = (*bndR)(i);
  }

  assert(tempItersIter == tempIters.end());
  for (int i = 0; i < nBnd; i++) {
    reconI({tempIters[i + 0], tempIters[i + 1], tempIters[i + 2],
            tempIters[i + 3], tempIters[i + 4], tempIters[i + 5]});
  }
  assert(iter == data->end());
}

void Recon5OrderFaceCenter::reconI(
    std::array<std::vector<real>::iterator, 6> input) {
  // 每一次运行iter都应该加2*nVar哟
  // 这里实现了原始变量重构
  for (int i = 0; i < nvar; i++) {
    *iter = inter5({
        input[0][i],
        input[1][i],
        input[2][i],
        input[3][i],
        input[4][i],
    });
    *(iter + nvar) = inter5({
        input[5][i],
        input[4][i],
        input[3][i],
        input[2][i],
        input[1][i],
    });
    iter++;
  }
  iter += nvar;
  int temp = data->end() - iter;
  std::cout << temp << '\n';
}

void Recon5Order1DEulerEig::reconI(
    std::array<std::vector<real>::iterator, 6> input) {
  // 每一次运行iter都应该加2*nVar哟
  // 这里实现了原始变量重构
  enum { R, U, P };
  auto primL = input[2], primR = input[3];
  eigensystemEuler1D eig = eigensystemEuler1D({primL[R], primL[U], primL[P]},
                                              {primR[R], primR[U], primR[P]});
  std::array<real, 5> q1L, q2L, q3L, q1R, q2R, q3R;
  for (int j = 0; j < 6; j++) {

    auto charTemp = eig.primToChar({input[j][R], input[j][U], input[j][P]});

    if (j < 5) {
      q1L[j] = charTemp[0];
      q2L[j] = charTemp[1];
      q3L[j] = charTemp[2];
    }

    int jj = 5 - j;
    if (jj < 5) {
      q1R[jj] = charTemp[0];
      q2R[jj] = charTemp[1];
      q3R[jj] = charTemp[2];
    }
  }
  // auto start = std::chrono::steady_clock::now();
  auto Q1LL = inter5(q1L);
  auto Q1RR = inter5(q1R);
  auto Q2LL = inter5(q2L);
  auto Q2RR = inter5(q2R);
  auto Q3LL = inter5(q3L);
  auto Q3RR = inter5(q3R);
  // auto stop = std::chrono::steady_clock::now();
  // auto duration =
  // std::chrono::duration_cast<std::chrono::nanoseconds>(stop
  // - start).count(); timep+=duration;

  auto resTempL = eig.charToPrim({Q1LL, Q2LL, Q3LL});
  auto resTempR = eig.charToPrim({Q1RR, Q2RR, Q3RR});
  std::copy(resTempL.begin(), resTempL.end(), iter);
  std::copy(resTempR.begin(), resTempR.end(), iter);
  this->iter += 2 * this->nvar;
}

void Recon5Order2DEulerEig::reconI(
    std::array<std::vector<real>::iterator, 6> input) {
  // 每一次运行iter都应该加2*nVar哟
  // 这里实现了原始变量重构
  assert(this->nvar == 4);
  enum { R, U, V, P };
  auto primL = input[2], primR = input[3];
  eigensystemEuler2D eig =
      eigensystemEuler2D({primL[R], primL[U], primL[V], primL[P]},
                         {primR[R], primR[U], primR[V], primR[P]}, this->norm);
  std::array<real, 5> q1L, q2L, q3L, q4L, q1R, q2R, q3R, q4R;
  for (int j = 0; j < 6; j++) {
    auto charTemp =
        eig.primToChar({input[j][R], input[j][U], input[j][V], input[j][P]});

    if (j < 5) {
      q1L[j] = charTemp[R];
      q2L[j] = charTemp[U];
      q3L[j] = charTemp[V];
      q4L[j] = charTemp[P];
    }

    int jj = 5 - j;
    if (jj < 5) {
      q1R[jj] = charTemp[R];
      q2R[jj] = charTemp[U];
      q3R[jj] = charTemp[V];
      q4R[jj] = charTemp[P];
    }
  }
  // auto start = std::chrono::steady_clock::now();
  auto Q1LL = inter5(q1L);
  auto Q1RR = inter5(q1R);
  auto Q2LL = inter5(q2L);
  auto Q2RR = inter5(q2R);
  auto Q3LL = inter5(q3L);
  auto Q3RR = inter5(q3R);
  auto Q4LL = inter5(q4L);
  auto Q4RR = inter5(q4R);
  // auto stop = std::chrono::steady_clock::now();
  // auto duration =
  // std::chrono::duration_cast<std::chrono::nanoseconds>(stop
  // - start).count(); timep+=duration;

  auto resTempL = eig.charToPrim({Q1LL, Q2LL, Q3LL, Q4LL});
  auto resTempR = eig.charToPrim({Q1RR, Q2RR, Q3RR, Q4RR});
  memcpy(&this->iter[R], &resTempL[0], this->nvar * sizeof(real));
  memcpy(&this->iter[P + 1], &resTempR[0], this->nvar * sizeof(real));
  this->iter += 2 * this->nvar;
}
#pragma once
#include "eigenSystem.hpp"
#include "reconstructor.hpp"

typedef real (*InterpolationScheme5Order)(std::array<real, 5>);

template <InterpolationScheme5Order inter5>
class Recon5OrderFaceCenter : public Reconstuctor {
public:
    Recon5OrderFaceCenter() {};

    Recon5OrderFaceCenter(DataReader varsR_, std::shared_ptr<OneDBnd> bndL_,
        std::shared_ptr<OneDBnd> bndR_)
        : Reconstuctor(varsR_, bndL_, bndR_) {};

protected:
    void initVarsR() override;
    void leftBnd() override;
    void internal() override;
    void rightBnd() override;
    virtual void reconI(std::array<std::vector<real>::iterator, 6>);
    int tt = 0;
};

template <InterpolationScheme5Order inter5>
class Recon5Order1DEulerEig : public Recon5OrderFaceCenter<inter5> {

public:
    Recon5Order1DEulerEig() {};
    Recon5Order1DEulerEig(DataReader varsR_, std::shared_ptr<OneDBnd> bndL_,
        std::shared_ptr<OneDBnd> bndR_)
        : Recon5OrderFaceCenter<inter5>(varsR_, bndL_, bndR_) {};

protected:
    const int NVAR = 3;
    void reconI(std::array<std::vector<real>::iterator, 6>) final;
};

template <InterpolationScheme5Order inter5>
class Recon5Order2DEulerEig : public Recon5OrderFaceCenter<inter5> {
public:
    Recon5Order2DEulerEig() {};
    Recon5Order2DEulerEig(DataReader varsR_, std::shared_ptr<OneDBnd> bndL_,
        std::shared_ptr<OneDBnd> bndR_)
        : Recon5OrderFaceCenter<inter5>(varsR_, bndL_, bndR_) {};

protected:
    const int NVAR = 4;
    void reconI(std::array<std::vector<real>::iterator, 6>) final;
};

template <InterpolationScheme5Order inter5>
void Recon5OrderFaceCenter<inter5>::initVarsR()
{
    if (!data) {
        int nPointR = bndL->getN() + n + bndR->getN() - 5;
        data = std::make_shared<Data>(nPointR, nvar * 2);
    }
}

template <InterpolationScheme5Order inter5>
void Recon5OrderFaceCenter<inter5>::leftBnd()
{
    // std::cout << (iter - data->end()) / 2 / nvar << '\n';
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
        reconI({ tempIters[i + 0], tempIters[i + 1], tempIters[i + 2],
            tempIters[i + 3], tempIters[i + 4], tempIters[i + 5] });
    }
}

template <InterpolationScheme5Order inter5>
void Recon5OrderFaceCenter<inter5>::internal()
{

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

template <InterpolationScheme5Order inter5>
void Recon5OrderFaceCenter<inter5>::rightBnd()
{
    // std::cout << (iter - data->end()) / nvar / 2 << '\n';
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
        reconI({ tempIters[i + 0], tempIters[i + 1], tempIters[i + 2],
            tempIters[i + 3], tempIters[i + 4], tempIters[i + 5] });
    }

    assert(iter == data->end());
}

template <InterpolationScheme5Order inter5>
void Recon5OrderFaceCenter<inter5>::reconI(
    std::array<std::vector<real>::iterator, 6> input)
{
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
}

template <InterpolationScheme5Order inter5>
void Recon5Order1DEulerEig<inter5>::reconI(
    std::array<std::vector<real>::iterator, 6> input)
{
    // 每一次运行iter都应该加2*nVar哟
    enum { R,
        U,
        P };
    auto primL = std::span<real, 3>(input[2], 3), primR = std::span<real, 3>(input[3], 3);
    eigensystemEuler1D eig = eigensystemEuler1D(primL,
        primR, { 1., 0., 0. });
    std::array<real, 5> q1L, q2L, q3L, q1R, q2R, q3R;
    for (int j = 0; j < 6; j++) {

        auto primj = std::span<real, 3>(input[j], 3);
        auto charTemp = eig.primToChar(primj);

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
    std::array<real, 3> charl = { Q1LL, Q2LL, Q3LL };
    std::array<real, 3> charr = { Q1RR, Q2RR, Q3RR };

    auto resTempL = eig.charToPrim(std::span<real, 3>(charl.begin(), 3));
    auto resTempR = eig.charToPrim(std::span<real, 3>(charl.begin(), 3));
    std::copy(resTempL.begin(), resTempL.end(), this->iter);
    std::copy(resTempR.begin(), resTempR.end(), this->iter + this->nvar);
    this->iter += 2 * this->nvar;
}

template <InterpolationScheme5Order inter5>
void Recon5Order2DEulerEig<inter5>::reconI(
    std::array<std::vector<real>::iterator, 6> input)
{
    // 每一次运行iter都应该加2*nVar哟
    assert(this->nvar == 4);
    enum { R,
        U,
        V,
        P };
    auto primL = std::span<real, 4>(input[2], 4), primR = std::span<real, 4>(input[3], 4);
    eigensystemEuler2D eig = eigensystemEuler2D(primL,
        primR, this->norm);
    std::array<real, 5> q1L, q2L, q3L, q4L, q1R, q2R, q3R, q4R;
    for (int j = 0; j < 6; j++) {
        auto primj = std::span<real, 4>(input[j], 4);
        auto charTemp = eig.primToChar(primj);

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

    std::array<real, 4> charl = { Q1LL, Q2LL, Q3LL, Q4LL };
    std::array<real, 4> charr = { Q1RR, Q2RR, Q3RR, Q4RR };

    auto resTempL = eig.charToPrim(std::span<real, 4>(charl.begin(), 4));
    auto resTempR = eig.charToPrim(std::span<real, 4>(charr.begin(), 4));
    std::copy(resTempL.begin(), resTempL.end(), this->iter);
    std::copy(resTempR.begin(), resTempR.end(), this->iter + this->nvar);
    this->iter += 2 * this->nvar;
}
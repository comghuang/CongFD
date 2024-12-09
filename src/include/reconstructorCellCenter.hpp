#pragma once
#include "eigenSystem.hpp"
#include "reconstructor.hpp"

using InterpolationScheme5OrderCell = std::array<real, 2> (*)(std::span<real, 5>);

template <InterpolationScheme5OrderCell inter5>
class Recon5OrderCellCenter : public Reconstuctor {
public:
    Recon5OrderCellCenter() {};

    Recon5OrderCellCenter(DataReader varsR_, std::shared_ptr<OneDBnd> bndL_,
        std::shared_ptr<OneDBnd> bndR_)
        : Reconstuctor(varsR_, bndL_, bndR_) {};

protected:
    void initVarsR() override;
    void leftBnd() override;
    void internal() override;
    void rightBnd() override;
    virtual void reconI(std::array<std::vector<real>::iterator, 5>);
    virtual void reconFirst(std::array<std::vector<real>::iterator, 5>);
    virtual void reconLast(std::array<std::vector<real>::iterator, 5>);
};

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------*/

template <InterpolationScheme5OrderCell inter5>
void Recon5OrderCellCenter<inter5>::initVarsR()
{
    if (!data) {
        int nPointR = bndL->getN() + n + bndR->getN() - 5;
        data = std::make_shared<Data>(nPointR, nvar * 2);
    }
}

template <InterpolationScheme5OrderCell inter5>
void Recon5OrderCellCenter<inter5>::leftBnd()
{

    int nBnd = bndL->getN();
    std::vector<std::vector<real>::iterator> tempIters(nBnd + 4);
    auto tempItersIter = tempIters.begin();
    for (int i = 0; i < nBnd; i++) {
        int iInver = nBnd - 1 - i;
        auto tempppp = (*bndL)(iInver);
        (*tempItersIter++) = (*bndL)(iInver);
    }
    for (int i = 0; i < 4; i++) {
        (*tempItersIter++) = varsReader(i);
    }
    assert(tempItersIter == tempIters.end());

    int i = 0;
    reconFirst({ tempIters[i + 0], tempIters[i + 1], tempIters[i + 2],
        tempIters[i + 3], tempIters[i + 4] });

    for (i = 1; i < nBnd; i++) {
        reconI({ tempIters[i + 0], tempIters[i + 1], tempIters[i + 2],
            tempIters[i + 3], tempIters[i + 4] });
    }
    // std::cout << iter - data->begin() << '\n';
}

template <InterpolationScheme5OrderCell inter5>
void Recon5OrderCellCenter<inter5>::internal()
{
    for (int i = 0; i < n - 4; i++) {
        auto iter = varsReader(i);
        reconI({ varsReader(i + 0),
            varsReader(i + 1),
            varsReader(i + 2),
            varsReader(i + 3),
            varsReader(i + 4) });
    }
    // std::cout << iter - data->end() << '\n';
}

template <InterpolationScheme5OrderCell inter5>
void Recon5OrderCellCenter<inter5>::rightBnd()
{
    int nBnd = bndR->getN();
    std::vector<std::vector<real>::iterator> tempIters(nBnd + 4);
    auto tempItersIter = tempIters.begin();
    for (int i = 0; i < 4; i++) {
        (*tempItersIter++) = varsReader(n - 4 + i);
    }
    for (int i = 0; i < nBnd; i++) {
        (*tempItersIter++) = (*bndR)(i);
    }

    assert(tempItersIter == tempIters.end());
    int i = 0;
    for (; i < nBnd - 1; i++) {
        reconI({ tempIters[i + 0], tempIters[i + 1], tempIters[i + 2],
            tempIters[i + 3], tempIters[i + 4] });
    }
    i = nBnd - 1;
    reconLast({ tempIters[i + 0], tempIters[i + 1], tempIters[i + 2],
        tempIters[i + 3], tempIters[i + 4] });
    // std::cout << iter - data->end() << '\n';
    assert(iter == data->end());
}

template <InterpolationScheme5OrderCell inter5>
void Recon5OrderCellCenter<inter5>::reconI(
    std::array<std::vector<real>::iterator, 5> input)
{
    // 每一次运行iter都应该加2*nVar哟
    // 这里实现了原始变量重构
    for (int i = 0; i < nvar; i++) {
        std::array<real, 5> primTemp = { input[0][i],
            input[1][i],
            input[2][i],
            input[3][i],
            input[4][i] };
        auto valueLR = inter5(primTemp);
        *iter = valueLR[0];
        *(iter + nvar) = valueLR[1];
        iter++;
    }

    iter += nvar;
}
template <InterpolationScheme5OrderCell inter5>
void Recon5OrderCellCenter<inter5>::reconFirst(
    std::array<std::vector<real>::iterator, 5> input)
{
    // 每一次运行iter都应该加2*nVar哟
    // 这里实现了原始变量重构
    for (int i = 0; i < nvar; i++) {
        std::array<real, 5> primTemp = { input[0][i],
            input[1][i],
            input[2][i],
            input[3][i],
            input[4][i] };
        auto valueLR = inter5(primTemp);
        *iter = valueLR[1];
        iter++;
    }
}

template <InterpolationScheme5OrderCell inter5>
void Recon5OrderCellCenter<inter5>::reconLast(
    std::array<std::vector<real>::iterator, 5> input)
{
    // 每一次运行iter都应该加2*nVar哟
    // 这里实现了原始变量重构
    for (int i = 0; i < nvar; i++) {
        std::array<real, 5> primTemp = { input[0][i],
            input[1][i],
            input[2][i],
            input[3][i],
            input[4][i] };
        auto valueLR = inter5(primTemp);
        *iter = valueLR[0];
        iter++;
    }
}

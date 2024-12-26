#pragma once
#include "reconstructorCellCenter.hpp"
#include <boost/circular_buffer.hpp>
#include <fstream>

template <std::size_t nSten>
using TenoDetector = uint_fast8_t (*)(std::span<real, nSten>);
inline std::array<real, 2> THINC(real q1, real q2, real q3);
inline std::array<real, 2> TENO_N_THINC(uint_fast8_t flag, std::span<real, 5> q);

template <std::size_t NVar, EigenSystem<NVar> eigenSystem, std::size_t nSten = 5, TenoDetector<nSten> tenoDetector = Teno5_Congf>
class ReconElement {
public:
    ReconElement() {};
    ReconElement(std::span<real, NVar> prim_, const std::array<real, 3>& norm_, const std::array<std::vector<real>::iterator, nSten>& stencils);
    inline void solve(std::span<real> dest);
    std::array<real, NVar> flags;

private:
    eigenSystem eigen;
    std::array<real, NVar * nSten> chars;
};

template <std::size_t NVar, EigenSystem<NVar> eigenSystem, std::size_t nSten = 5, TenoDetector<nSten> tenoDetector = Teno5_Congf>
class ReconerTenoThinc : public Recon5OrderCellCenter<nullptr> {
public:
    ReconerTenoThinc() {};
    ReconerTenoThinc(DataReader varsR_, std::shared_ptr<OneDBnd> bndL_,
        std::shared_ptr<OneDBnd> bndR_)
        : Recon5OrderCellCenter(varsR_, bndL_, bndR_) {};
    // void solve() final;

protected:
    void initVarsR() final;
    void leftBnd() final;
    void internal() final;
    void rightBnd() final;
    void reconI(std::array<std::vector<real>::iterator, 5>);

    void outputFlags(const std::array<real, NVar> flags)
    {
        return;
        std::fstream file("flags.dat", std::ios::app);
        file << std::format("{} {} {} {}\n", count++, flags[0], flags[1], flags[2]);
        file.close();
    }
    int count = 0;

    using Element
        = ReconElement<NVar, eigenSystem, nSten, tenoDetector>;
    boost::circular_buffer<Element> buffer;
};

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------*/
template <std::size_t NVar, EigenSystem<NVar> eigenSystem, std::size_t nSten, TenoDetector<nSten> tenoDetector>
ReconElement<NVar, eigenSystem, nSten, tenoDetector>::ReconElement(std::span<real, NVar> prim_, const std::array<real, 3>& norm_, const std::array<std::vector<real>::iterator, nSten>& stencils)
    : eigen(prim_, norm_)
{
    for (auto isten : std::ranges::views::iota(0ul, nSten)) {
        auto primi = std::span<real, NVar>(stencils[isten], NVar);
        std::array<real, NVar> chartemp = eigen.primToChar(primi);
        for (auto ivar : std::ranges::views::iota(0ul, NVar)) {
            chars[ivar * nSten + isten] = chartemp[ivar];
        }
    }
    for (auto ivar : std::ranges::views::iota(0ul, NVar)) {
        auto interVars = std::span<real, nSten>(chars.begin() + ivar * nSten, nSten);
        flags[ivar] = tenoDetector(interVars);
    }
}
template <std::size_t NVar, EigenSystem<NVar> eigenSystem, std::size_t nSten, TenoDetector<nSten> tenoDetector>
inline void ReconElement<NVar, eigenSystem, nSten, tenoDetector>::solve(std::span<real> dest)
{
    assert(dest.size() == NVar * 2);
    std::array<real, NVar> charL, charR;
    for (auto ivar : std::ranges::views::iota(0ul, NVar)) {
        // for (auto ivar = 0; ivar < NVar; ivar++) {
        auto interVars = std::span<real, nSten>(chars.begin() + ivar * nSten, nSten);
        std::array<real, 2> res = TENO_N_THINC(flags[ivar], interVars);
        charL[ivar] = res[0];
        charR[ivar] = res[1];
    }
    std::array<real, NVar> primL = eigen.charToPrim(charL);
    std::array<real, NVar> primR = eigen.charToPrim(charR);
    std::copy(primL.begin(), primL.end(), dest.begin());
    std::copy(primR.begin(), primR.end(), dest.begin() + NVar);
}

template <std::size_t NVar, EigenSystem<NVar> eigenSystem, std::size_t nSten, TenoDetector<nSten> tenoDetector>
void ReconerTenoThinc<NVar, eigenSystem, nSten, tenoDetector>::initVarsR()
{
    if (!data) {
        int nPointR = bndL->getN() + n + bndR->getN() - 5;
        data = std::make_shared<Data>(nPointR, nvar * 2);
    }
    buffer = boost::circular_buffer<Element>(7);
    buffer.clear();
    count = 0;
    // std::fstream file("flags.dat", std::ios::out);
    // file.close();
}

template <std::size_t NVar, EigenSystem<NVar> eigenSystem, std::size_t nSten, TenoDetector<nSten> tenoDetector>
void ReconerTenoThinc<NVar, eigenSystem, nSten, tenoDetector>::leftBnd()
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

    for (int i = 0; i < nBnd; i++) {
        reconI({ tempIters[i + 0], tempIters[i + 1], tempIters[i + 2],
            tempIters[i + 3], tempIters[i + 4] });
    }
    // std::cout << iter - data->begin() << '\n';
}

template <std::size_t NVar, EigenSystem<NVar> eigenSystem, std::size_t nSten, TenoDetector<nSten> tenoDetector>
void ReconerTenoThinc<NVar, eigenSystem, nSten, tenoDetector>::internal()
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

template <std::size_t NVar, EigenSystem<NVar> eigenSystem, std::size_t nSten, TenoDetector<nSten> tenoDetector>
void ReconerTenoThinc<NVar, eigenSystem, nSten, tenoDetector>::rightBnd()
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
    for (int i = 0; i < nBnd; i++) {
        reconI({ tempIters[i + 0], tempIters[i + 1], tempIters[i + 2],
            tempIters[i + 3], tempIters[i + 4] });
    }

    // 对buffer中的所有元素执行solve，并抛弃最后一个重构单元的第二个返回值
    for (auto iEle = buffer.begin(); iEle != buffer.end() - 1; iEle++) {
        iEle->solve(std::span<real, NVar * 2>(iter, NVar * 2));
        outputFlags(iEle->flags); // output
        iter += NVar * 2;
    }
    std::array<real, NVar * 2> temp;
    buffer.back().solve(temp);
    outputFlags(buffer.back().flags); // output
    std::copy(temp.begin(), temp.begin() + NVar, iter);
    iter += NVar;
    assert(iter == data->end());
}

template <std::size_t NVar, EigenSystem<NVar> eigenSystem, std::size_t nSten, TenoDetector<nSten> tenoDetector>
void ReconerTenoThinc<NVar, eigenSystem, nSten, tenoDetector>::reconI(
    std::array<std::vector<real>::iterator, 5> input)
{
    // step1: if buffer.full()
    if (buffer.full()) {
        if (iter == data->begin()) {
            //      1.2 若iter在起始位置，抛弃第一个返回值
            std::array<real, NVar * 2> temp;
            buffer.front().solve(temp);
            outputFlags(buffer.front().flags); // output
            std::copy(temp.begin() + NVar, temp.end(), iter);
            iter += NVar;
        } else {
            //      1.1 若iter不在起始位置，对buffer的第一个元素执行solve。
            buffer.front().solve(std::span<real, NVar * 2>(iter, NVar * 2));
            outputFlags(buffer.front().flags); // output
            iter += NVar * 2;
        }
    }

    // step2: push_back新的重构单元，并执行TENO-THINC算法
    auto primi = std::span<real, NVar>(input[2], NVar);
    buffer.push_back(Element(primi, norm, input));
    //      2.1 for ivar : 0-NVar
    auto leftSign = [](uint_fast8_t flag) -> bool { uint_fast8_t sign = 0b100;return !((sign&flag)==0); }; // 判断最高位是否为1
    auto midSign = [](uint_fast8_t flag) -> bool { uint_fast8_t sign = 0b010;return !((sign&flag)==0); }; // 判断最高位是否为1
    auto rightSign = [](uint_fast8_t flag) -> bool { uint_fast8_t sign = 0b001;return !((sign&flag)==0); }; // 判断最低位是否为1
    for (auto ivar : std::ranges::views::iota(0ul, NVar)) {
        //          if buffer[end].flag[ivar]的第一位为0
        if (rightSign(buffer.back().flags[ivar])) {
            //           for iEle=end-1:-1(反向遍历)
            for (auto iEle = buffer.rbegin() + 1; iEle != buffer.rend(); iEle++) {
                //            if iEle的第三位为0
                if (leftSign(iEle->flags[ivar])) {
                    for (auto iEle2 = buffer.rbegin() + 1; iEle2 != iEle; iEle2++) {
                        {
                            // 令iEle的与end之间的重构单元flag[ivar]==0的flag[ivar]=7,即进入THINC重构，并退出循环
                            if (iEle2->flags[ivar] == 0) {
                                // iEle2->flags[ivar] = 7;
                                // (iEle2 + 1)->flags[ivar] = 7;
                                // (iEle2 - 1)->flags[ivar] = 7;
                            }
                        }
                    }
                    break;
                }
            }
        }
    }
}

// template <std::size_t NVar, EigenSystem<NVar> eigenSystem, std::size_t nSten, TenoDetector<nSten> tenoDetector>
// void ReconerTenoThinc<NVar, eigenSystem, nSten, tenoDetector>::solve()
// {
//     initVarsR();
//     for (ivar = 0; ivar < varsReader.getNVar(); ivar++) {
//         iterFlags = flags.begin();
//         leftBnd();
//         internal();
//         rightBnd();
//         finalReconstruction();
//     }
// }

// template <std::size_t NVar, EigenSystem<NVar> eigenSystem, std::size_t nSten, TenoDetector<nSten> tenoDetector>
// void ReconerTenoThinc<NVar, eigenSystem, nSten, tenoDetector>::initVarsR()
// {
//     int nPointR = bndL->getN() + n + bndR->getN() - 5;
//     if (!data) {
//         data = std::make_shared<Data>(nPointR, nvar * 2);
//     }

//     if (flags.size() != nPointR) {
//         flags.resize(nPointR);
//     }
// }

// template <std::size_t NVar, EigenSystem<NVar> eigenSystem, std::size_t nSten, TenoDetector<nSten> tenoDetector>
// void ReconerTenoThinc<NVar, eigenSystem, nSten, tenoDetector>::leftBnd()
// {
//     int nBnd = bndL->getN();
//     std::vector<std::vector<real>::iterator> tempIters(nBnd + 5);
//     auto tempItersIter = tempIters.begin();
//     for (int i = 0; i < nBnd; i++) {
//         int iInver = nBnd - 1 - i;
//         auto tempppp = (*bndL)(iInver);
//         (*tempItersIter++) = (*bndL)(iInver) + ivar;
//     }
//     for (int i = 0; i < 5; i++) {
//         (*tempItersIter++) = varsReader(i) + ivar;
//     }
//     assert(tempItersIter == tempIters.end());

//     for (i = 0; i < nBnd; i++) {
//         (*iterFlags) = Teno5_Congf({ *(tempIters[i + 0]),
//             *(tempIters[i + 1]),
//             *(tempIters[i + 2]),
//             *(tempIters[i + 3]),
//             *(tempIters[i + 4]) });
//         iterFlags++;
//     }
// }

// template <std::size_t NVar, EigenSystem<NVar> eigenSystem, std::size_t nSten, TenoDetector<nSten> tenoDetector>
// void ReconerTenoThinc<NVar, eigenSystem, nSten, tenoDetector>::internal()
// {
//     for (int i = 0; i < n - 4; i++) {
//         (*iterFlags) = Teno5_Congf({ *(varsReader(i + 0) + ivar),
//             *(varsReader(i + 1) + ivar),
//             *(varsReader(i + 2) + ivar),
//             *(varsReader(i + 3) + ivar),
//             *(varsReader(i + 4) + ivar) });
//         iterFlags++;
//     }
// }

// template <std::size_t NVar, EigenSystem<NVar> eigenSystem, std::size_t nSten, TenoDetector<nSten> tenoDetector>
// void ReconerTenoThinc<NVar, eigenSystem, nSten, tenoDetector>::rightBnd()
// {
//     int nBnd = bndR->getN();
//     std::vector<std::vector<real>::iterator> tempIters(nBnd + 4);
//     auto tempItersIter = tempIters.begin();
//     for (int i = 0; i < 4; i++) {
//         (*tempItersIter++) = varsReader(n - 4 + i) + ivar;
//     }
//     for (int i = 0; i < nBnd; i++) {
//         (*tempItersIter++) = (*bndR)(i) + ivar;
//     }
//     assert(tempItersIter == tempIters.end());

//     for (int i = 0; i < nBnd; i++) {
//         (*iterFlags) = Teno5_Congf({ *(tempIters[i + 0]),
//             *(tempIters[i + 1]),
//             *(tempIters[i + 2]),
//             *(tempIters[i + 3]),
//             *(tempIters[i + 4]) });
//         iterFlags++;
//     }
//     assert(iterFlags == flags.end());
//     // assert(iter == data->end());
// }

inline std::array<real, 2> TENO_N_THINC(uint_fast8_t flag, std::span<real, 5> q)
{
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return { 3.0 / 128.0 * q[4] - 5.0 / 32.0 * q[3] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[1] - 5.0 / 128.0 * q[0],
            3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4] };
        break;
    case 1:
        /* 0,1,1 */
        return { 1.0 / 16.0 * q[4] - 5.0 / 16.0 * q[3] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[1],
            -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4] };
        break;
    case 2:
        /* 1,0,1 */
        return { 3.0 / 8.0 * q[4] - 5.0 / 4.0 * q[3] + 15.0 / 8.0 * q[2],
            3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4] };
        break;
    case 4:
        /* 1,1,0 */
        return { -1.0 / 16.0 * q[3] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[1] - 1.0 / 16.0 * q[0],
            1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3] };
        break;
    case 3:
        /* 0,0,1 */
        return { 3.0 / 8.0 * q[4] - 5.0 / 4.0 * q[3] + 15.0 / 8.0 * q[2],
            3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4] };
        break;
    case 5:
        /* 0,1,0 */
        return { -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3],
            -1.0 / 8.0 * q[3] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[1] };
        break;
    case 6:
        /* 1,0,0 */
        return { 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[1] - 1.0 / 8.0 * q[0],
            3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2] };
        break;
    case 7:
        return THINC(q[1], q[2], q[3]);
        break;
    default:
        /* 0,0,0 */
        return { q[2], q[2] };
        break;
    }
}
inline std::array<real, 2> THINC(real q1, real q2, real q3)
{

    if ((q1 - q2) * (q2 - q3) < 1e-20)
        return { q2, q2 };

    real qmax, qmin;
    if (q1 > q3) {
        qmax = q1;
        qmin = q3;
    } else {
        qmax = q3;
        qmin = q1;
    }
    constexpr real beta = 0.70;
    constexpr real T1 = std::exp(2.0 * beta), T3 = exp(-2.0 * beta);

    constexpr real eps = 1e-20;
    real qb = (q1 + q3) / 2, dq = (q3 - q1) / 2;
    real pp = (q2 - qb) / (dq);
    assert(!std::isnan(pp));
    return { qb + dq * (pp - 1.0 + (pp + 1.0) * T3) / (1.0 - pp + (pp + 1.0) * T3),
        qb + dq * (pp - 1.0 + (pp + 1.0) * T1) / (1.0 - pp + (pp + 1.0) * T1) };
}
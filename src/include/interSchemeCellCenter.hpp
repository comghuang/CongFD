#pragma once
#include "macro.hpp"
#include <cstdint>

inline std::array<real, 2> weno5_JSCell(std::span<real, 5> q);
inline std::array<real, 2> weno5_ZCell(std::span<real, 5> q);
inline std::array<real, 2> teno5_CongCell(std::span<real, 5> q);
inline std::array<real, 2> teno5_ZCell(std::span<real, 5> q);

inline uint_fast8_t Teno5_Congf(std::span<real, 5> q);
inline uint_fast8_t Teno5_Zf(std::span<real, 5> q);
inline std::array<real, 2> TENO_N(uint_fast8_t flag, std::span<real, 5> q);
inline std::array<real, 2> THINC(real q1, real q2, real q3);
/*-----------------------------------------------------------------------------*/

inline std::array<real, 2> teno5_CongCell(std::span<real, 5> q)
{
    return TENO_N(Teno5_Congf(q), q);
}
inline std::array<real, 2> teno5_ZCell(std::span<real, 5> q)
{
    return TENO_N(Teno5_Congf(q), q);
}
inline std::array<real, 2> weno5_JSCell(std::span<real, 5> q)
{
    //   leftface|    cell i     |rightface
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    auto func = [](real q0, real q1, real q2, real q3, real q4, std::array<real, 3> beta) -> real {
        std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
        real res = 0;
        real sumbeta = 0;
        real eps = 1e-6;
        std::array<real, 3> alpha;
        std::array<real, 3> u;
        u[0] = 3.0 / 8.0 * q0 - 5.0 / 4.0 * q1 + 15.0 / 8.0 * q2;
        u[1] = -1.0 / 8.0 * q1 + 3.0 / 4.0 * q2 + 3.0 / 8.0 * q3;
        u[2] = 3.0 / 8.0 * q2 + 3.0 / 4.0 * q3 - 1.0 / 8.0 * q4;
        for (int i = 0; i < 3; i++) {
            alpha[i] = gamma[i] / std::pow(eps + beta[i], 2);
            sumbeta += alpha[i];
        }
        for (int i = 0; i < 3; i++)
            res += alpha[i] * u[i];
        res /= sumbeta;
        return res;
    };
    real resL = func(q[4], q[3], q[2], q[1], q[0], { beta[2], beta[1], beta[0] });
    real resR = func(q[0], q[1], q[2], q[3], q[4], { beta[0], beta[1], beta[2] });

    return { resL, resR };
}

inline std::array<real, 2> weno5_ZCell(std::span<real, 5> q)
{
    //   leftface|    cell i     |rightface
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    auto func = [](real q0, real q1, real q2, real q3, real q4, std::array<real, 3> beta) -> real {
        std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
        real res = 0;
        real sumbeta = 0;
        real eps = 1e-6;
        std::array<real, 3> alpha;
        std::array<real, 3> u;
        u[0] = 3.0 / 8.0 * q0 - 5.0 / 4.0 * q1 + 15.0 / 8.0 * q2;
        u[1] = -1.0 / 8.0 * q1 + 3.0 / 4.0 * q2 + 3.0 / 8.0 * q3;
        u[2] = 3.0 / 8.0 * q2 + 3.0 / 4.0 * q3 - 1.0 / 8.0 * q4;
        real C = 1, qq = 2, tau = std::abs(beta[2] - beta[0]);
        for (int i = 0; i < 3; i++) {
            alpha[i] = gamma[i] * (C + pow(tau / (beta[i] + eps), qq));
            sumbeta += alpha[i];
        }
        for (int i = 0; i < 3; i++)
            res += alpha[i] * u[i];
        res /= sumbeta;
        return res;
    };
    real resL = func(q[4], q[3], q[2], q[1], q[0], { beta[2], beta[1], beta[0] });
    real resR = func(q[0], q[1], q[2], q[3], q[4], { beta[0], beta[1], beta[2] });

    return { resL, resR };
}

inline uint_fast8_t Teno5_Congf(std::span<real, 5> q)
{
    //   leftface|    cell i     |rightface
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    uint_fast8_t minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    constexpr real CT = std::pow(1.5 * 1e-5, 1.0 / 6.0);
    constexpr real CT_1 = 1 - CT;
    real tau = std::abs(beta[2] - beta[0]);
    real rr = CT * tau - CT_1 * beta[minBeta];
    real ll = tau * beta[minBeta];
    uint_fast8_t flag = 0;
    if (ll < rr * beta[0])
        flag += 1;
    if (ll < rr * beta[1])
        flag += 2;
    if (ll < rr * beta[2])
        flag += 4;

    return flag;
}

inline uint_fast8_t Teno5_Zf(std::span<real, 5> q)
{
    real eps = 1e-40;
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    real sumbeta = 0;
    real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
    for (int i = 0; i < 3; i++) {
        real tempp = C + tau / (beta[i] + eps);
        tempp *= tempp;
        beta[i] = tempp * tempp * tempp;
        sumbeta += beta[i];
    }
    real CT = 1e-5 * sumbeta;
    // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);
    uint_fast8_t flag = 0;
    if (beta[0] < CT)
        flag += 1;
    if (beta[1] < CT)
        flag += 2;
    if (beta[2] < CT)
        flag += 4;
    return flag;
}

inline std::array<real, 2> TENO_N(uint_fast8_t flag, std::span<real, 5> q)
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
    default:
        /* 0,0,0 */
        return { q[2], q[2] };
        break;
    }
}
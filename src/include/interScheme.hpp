#pragma once
#include "macro.hpp"
#include <array>

inline real u1(real q1, real q2, real q3);
inline real u2(real q1, real q2, real q3);
inline real u3(real q1, real q2, real q3);
inline real weno5_JSchen(std::array<real, 5>);
inline real weno5_Z(std::array<real, 5> q);

inline real Teno5_ZCT4(std::array<real, 5> q);
inline real Teno5_ZCT7(std::array<real, 5> q);
inline real Teno5_Z(std::array<real, 5> q);
inline real Teno5_ZConvex(std::array<real, 5> q);
inline real Teno5_CongZ(std::array<real, 5> q);
inline real Teno5_CongZCT4(std::array<real, 5> q);
inline real Teno5_CongZCT7(std::array<real, 5> q);
inline real Teno5_CongA(std::array<real, 5> q);
inline real musclInterpolation(real q1, real q2, real q3);

inline real THINC1(real q1, real q2, real q3);
inline real THINC2(real q1, real q2, real q3);

inline real Linear5(std::array<real, 5> q);
/*-------------------------------------------------------------------------*/

inline real u1(real q1, real q2, real q3)
{
    return 3.0 / 8.0 * q1 - 5.0 / 4.0 * q2 + 15.0 / 8.0 * q3;
}
inline real u2(real q1, real q2, real q3)
{
    return -1.0 / 8.0 * q1 + 3.0 / 4.0 * q2 + 3.0 / 8.0 * q3;
}
inline real u3(real q1, real q2, real q3)
{
    return 3.0 / 8.0 * q1 + 3.0 / 4.0 * q2 - 1.0 / 8.0 * q3;
}

inline real weno5_JSchen(std::array<real, 5> q)
{
    real eps = 1e-6;
    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    std::array<real, 3> u;
    u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

    // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};

    real sumbeta = 0, result = 0;
    for (int i = 0; i < 3; i++) {
        beta[i] = gamma[i] / pow(eps + beta[i], 2);
        sumbeta += beta[i];
    }
    for (int i = 0; i < 3; i++)
        result += beta[i] * u[i];
    return result / sumbeta;
}

inline real weno5_Z(std::array<real, 5> q)
{
    real eps = 1e-40;
    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    std::array<real, 3> u;
    u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

    // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};

    real sumbeta = 0, result = 0;
    real C = 1, qq = 2, tau = std::abs(beta[2] - beta[0]);
    for (int i = 0; i < 3; i++) {
        beta[i] = gamma[i] * (C + pow(tau / (beta[i] + eps), qq));
        sumbeta += beta[i];
    }
    for (int i = 0; i < 3; i++)
        result += beta[i] * u[i];
    return result / sumbeta;
}

inline real Teno5_ZCT4(std::array<real, 5> q)
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
    real CT = 1e-4 * sumbeta;
    // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);
    unsigned short flag = 0;
    if (beta[0] < CT)
        flag += 1;
    if (beta[1] < CT)
        flag += 2;
    if (beta[2] < CT)
        flag += 4;
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}

inline real Teno5_ZCT7(std::array<real, 5> q)
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
    real CT = 1e-7 * sumbeta;
    // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);
    unsigned short flag = 0;
    if (beta[0] < CT)
        flag += 1;
    if (beta[1] < CT)
        flag += 2;
    if (beta[2] < CT)
        flag += 4;
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}

inline real Teno5_Z(std::array<real, 5> q)
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
    unsigned short flag = 0;
    if (beta[0] < CT)
        flag += 1;
    if (beta[1] < CT)
        flag += 2;
    if (beta[2] < CT)
        flag += 4;
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}
inline real Teno5_ZConvex(std::array<real, 5> q)
{
    real eps = 1e-40;
    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};
    std::array<real, 3> u;
    u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

    real sumbeta = 0, result = 0;
    real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
    for (int i = 0; i < 3; i++) {
        real tempp = C + tau / (beta[i] + eps);
        tempp *= tempp;
        beta[i] = tempp * tempp * tempp;
        sumbeta += beta[i];
    }

    real CT = 1e-5 * sumbeta, sumGamma = 0;
    for (int i = 0; i < 3; i++) {
        if (beta[i] > CT) {
            // result+=gamma[i]*(*u[i])(q[i],q[i+1],q[i+2]);
            result += gamma[i] * u[i];
            sumGamma += gamma[i];
        }
    }
    result /= sumGamma;
    return result;
}

const static real CTi = pow(1.5 * 1e-10, 1.0 / 6.0);
inline real Teno5_CongZ(std::array<real, 5> q)
{
    real eps = 1e-40;
    std::array<real, 3> beta = {
        1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),
        1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),
        1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)
    };
    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    constexpr real CT = std::pow(1.5 * 1e-5, 1.0 / 6.0);
    constexpr real CT_1 = 1 - CT;
    real tau = std::abs(beta[2] - beta[0]);
    real rr = CT * tau - CT_1 * beta[minBeta];
    real ll = tau * beta[minBeta];
    unsigned short flag = 0;
    if (ll < rr * beta[0])
        flag += 1;
    if (ll < rr * beta[1])
        flag += 2;
    if (ll < rr * beta[2])
        flag += 4;
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}

inline real Teno5_CongZCT4(std::array<real, 5> q)
{
    real eps = 1e-40;
    std::array<real, 3> beta = {
        1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

        1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

        1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)
    };

    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    constexpr real CT = std::pow(1.5 * 1e-4, 1.0 / 6.0);
    constexpr real CT_1 = 1 - CT;
    real tau = std::abs(beta[2] - beta[0]);
    real rr = CT * tau - CT_1 * beta[minBeta];
    real ll = tau * beta[minBeta];
    unsigned short flag = 0;
    if (ll < rr * beta[0])
        flag += 1;
    if (ll < rr * beta[1])
        flag += 2;
    if (ll < rr * beta[2])
        flag += 4;
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}
inline real Teno5_CongZCT7(std::array<real, 5> q)
{
    real eps = 1e-40;
    std::array<real, 3> beta = {
        1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),
        1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),
        1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)
    };
    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    constexpr real CT = std::pow(1.5 * 1e-7, 1.0 / 6.0);
    constexpr real CT_1 = 1 - CT;
    real tau = std::abs(beta[2] - beta[0]);
    real rr = CT * tau - CT_1 * beta[minBeta];
    real ll = tau * beta[minBeta];
    unsigned short flag = 0;
    if (ll < rr * beta[0])
        flag += 1;
    if (ll < rr * beta[1])
        flag += 2;
    if (ll < rr * beta[2])
        flag += 4;
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}

inline real Teno5_CongA(std::array<real, 5> q)
{
    real eps = 1e-40; // 1e-10;
    std::array<real, 3> beta = {
        1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

        1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

        1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)
    };
    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    real tau = std::abs(beta[2] - beta[0]);
    constexpr real CT = 0.1; // 6
    constexpr real CT_1 = 1 - CT;

    real mulbeta = beta[0] * beta[1] * beta[2];
    real rr = CT * tau * (beta[0] * beta[1] + beta[1] * beta[2] + beta[0] * beta[2]) - CT_1 * mulbeta;
    real ll = tau * mulbeta;
    unsigned short flag = 0;
    if (ll < rr * beta[0])
        flag += 1;
    if (ll < rr * beta[1])
        flag += 2;
    if (ll < rr * beta[2])
        flag += 4;
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}

inline real musclInterpolation(real q1, real q2, real q3)
{

    real delta;
    real deltam, deltap;
    deltam = q2 - q1;
    deltap = q3 - q2;

    // minmod
    real beta = 1.0;
    if (deltap > 0) {
        delta = std::max(0.0, std::max(std::min(beta * deltam, deltap), std::min(deltam, beta * deltap)));
    } else {
        delta = std::min(0.0, std::min(std::max(beta * deltam, deltap), std::max(deltam, beta * deltap)));
    }
    return q2 + delta * 0.5;
}

inline real THINC1(real q1, real q2, real q3)
{
    if ((q1 - q2) * (q2 - q3) < 1e-20)
        return q2;

    real qmax, qmin;
    if (q1 > q3) {
        qmax = q1;
        qmin = q3;
    } else {
        qmax = q3;
        qmin = q1;
    }
    constexpr real beta = 1.0;
    constexpr real T = std::exp(2.0 * beta);
    constexpr real eps = 1e-20;
    real qb
        = (q1 + q3) / 2,
        dq = (q3 - q1) / 2;
    real pp = (q2 - qb) / (dq);
    assert(!std::isnan(pp));
    real res = qb + dq * (pp - 1.0 + (pp + 1.0) * T) / (1.0 - pp + (pp + 1.0) * T);
    return res;
}

inline real THINC2(real q1, real q2, real q3)
{
    if ((q1 - q2) * (q2 - q3) < 1e-20)
        return q2;

    real qmax, qmin;
    if (q1 > q3) {
        qmax = q1;
        qmin = q3;
    } else {
        qmax = q3;
        qmin = q1;
    }
    constexpr real beta = 0.7;
    constexpr real T = std::exp(2.0 * beta);
    constexpr real eps = 1e-20;
    real qb
        = (q1 + q3) / 2,
        dq = (q3 - q1) / 2;
    real pp = (q2 - qb) / (dq);
    assert(!std::isnan(pp));
    real res = qb + dq * (pp - 1.0 + (pp + 1.0) * T) / (1.0 - pp + (pp + 1.0) * T);
    return res;
}

inline real Linear5(std::array<real, 5> q)
{
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
}

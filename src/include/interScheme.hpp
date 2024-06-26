#pragma once
#include <array>
#include "macro.hpp"
constexpr real weno5_JSchen(real,real,real,real,real);

constexpr real weno5_JSchen(real q1,real q2,real q3,real q4,real q5)
{
    real eps=1e-6;
    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta,u;
    beta[0]= 1.0/1.0 *pow(1.0*q1-2.0*q2+1.0*q3,2)
            + 1.0/4.0 *pow(1.0*q1-4.0*q2+3.0*q3,2);

    beta[1]= 1.0/1.0  *pow(1.0*q2-2.0*q3+1.0*q4,2)
            + 1.0/4.0 *pow(1.0*q2+0.0*q3-1.0*q4,2);

    beta[2]= 1.0/1.0 *pow(1.0*q3-2.0*q4+1.0*q5,2)
            + 1.0/4.0*pow(3.0*q3-4.0*q4+1.0*q5,2);
    
    u[0]= 3.0/8.0*q1-5.0/4.0*q2+15.0/8.0*q3;
    u[1]=-1.0/8.0*q2+3.0/4.0*q3+3.0 /8.0*q4;
    u[2]= 3.0/8.0*q3+3.0/4.0*q4-1.0 /8.0*q5;
    
    real sumbeta=0,result=0;
    for(int i=0;i<3;i++)
    {
        beta[i]=gamma[i]/pow(eps+beta[i],2.0);
        sumbeta+=beta[i];
    }
    for(int i=0;i<3;i++) result+=beta[i]/sumbeta*u[i];
    return result;
}
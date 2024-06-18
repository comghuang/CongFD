#include "SpaceDis.hpp"

std::vector<real> SpaceDis::difHCS()
{
    std::vector<real> rhs;
    rhs.resize(n*nVar,0);
    for(ind i=0;i<n;i++)
    {
        real h;
        h=2.0/n;
        for(ind j=0;j<nVar;j++)
        {
            rhs[i*nVar+j]=75.0/64.0*(flux(i+1,j)-flux(i,j))/h
                         -25.0/128.0*(flux(i+2,j)-flux(i-1,j))/(3*h)
                         +3.0/128.0*(flux(i+3,j)-flux(i-2,j))/(5*h);
            // rhs[i*nVar+j]=(flux(i+1,j)-flux(i,j))/h;
        }
    }
    return rhs;
}
std::vector<real> SpaceDis::difTraditional6()
{
    std::vector<real> rhs;
    rhs.resize(n*nVar,0);
    for(ind i=0;i<n;i++)
    {
        real h;
        h=2.0/n;
        for(ind j=0;j<nVar;j++)
        {
            rhs[i*nVar+j]=75.0/64.0*(flux(i+1,j)-flux(i,j))/h
                         -25.0/128.0*(flux(i+2,j)-flux(i-1,j))/(3*h)
                         +3.0/128.0*(flux(i+3,j)-flux(i-2,j))/(5*h);
            // rhs[i*nVar+j]=(flux(i+1,j)-flux(i,j))/h;
        }
    }
    return rhs;
}

std::vector<real> SpaceDis::dif2Order()
{
    std::vector<real> rhs;
    rhs.resize(n*nVar,0);
    for(ind i=0;i<n;i++)
    {
        real h;
        h=2.0/n;
        for(ind j=0;j<nVar;j++)
        {
            rhs[i*nVar+j]=(flux(i+1,j)-flux(i,j))/h;
        }
    }
    return rhs;
}
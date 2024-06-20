#include "SpaceDis.hpp"

void SpaceDis::difHCS()
{
    for(ind i=0;i<n;i++)
    {
        real h;
        h=2.0/n;
        for(ind j=0;j<nVar;j++)
        {
            (*rhs)(i0+i*offset,j)=75.0/64.0*(fAt(i+1,j)-fAt(i,j))/h
                         -25.0/128.0*(fAt(i+2,j)-fAt(i-1,j))/(3*h)
                         +3.0/128.0*(fAt(i+3,j)-fAt(i-2,j))/(5*h);
            // rhs[i*nVar+j]=(flux(i+1,j)-flux(i,j))/h;
        }
    }
}
void SpaceDis::difTraditional6()
{
    for(ind i=0;i<n;i++)
    {
        real h;
        h=2.0/n;
        for(ind j=0;j<nVar;j++)
        {
            (*rhs)(i0+i*offset,j)+=75.0/64.0*(fAt(i+1,j)-fAt(i,j))/h
                         -25.0/128.0*(fAt(i+2,j)-fAt(i-1,j))/(3*h)
                         +3.0/128.0*(fAt(i+3,j)-fAt(i-2,j))/(5*h);
            // rhs[i*nVar+j]=(flux(i+1,j)-flux(i,j))/h;
        }
    }
}

void SpaceDis::dif2Order()
{
    for(ind i=0;i<n;i++)
    {
        real h;
        h=2.0/n;
        for(ind j=0;j<nVar;j++)
        {
            (*rhs)(i0+i*offset,j)=(fAt(i+1,j)-fAt(i,j))/h;
        }
    }
}
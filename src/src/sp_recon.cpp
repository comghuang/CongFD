#include "SpaceDis.hpp"
#include "interScheme.hpp"
std::vector<real> SpaceDis::reconR(int i)
{
    std::vector<real> res(nPrim);
    if(interMethod==MUSCL)
    {
        for (int ivar = 0; ivar < nPrim; ivar++)
        {
            res[ivar]=musclInterpolation(at(i+1,ivar),at(i,ivar),at(i-1,ivar));
        }
        return res;
    }
    else 
    {
        for (int ivar = 0; ivar < nPrim; ivar++)
        {
            res[ivar]=at(i,ivar);
        }
        return res;
    }

}


std::vector<real> SpaceDis::reconL(int i)
{
    std::vector<real> res(nPrim);
    if(interMethod==MUSCL)
    {
        for (int ivar = 0; ivar < nPrim; ivar++)
        {
            res[ivar]=musclInterpolation(at(i-2,ivar),at(i-1,ivar),at(i,ivar));
        }
        return res;
    }
    else 
    {
        for (int ivar = 0; ivar < nPrim; ivar++)
        {
            res[ivar]=at(i-1,ivar);
        }
        return res;
    }
}

std::vector<real> SpaceDis::reconLChar1D(int i)
{
    assert(info->dim==1);
    assert(info->eqType==EULER);
    enum{R,U,P};
    real rRef=at(i-1,R);
    real pRef=at(i-1,P);
    real cRef=sqrt(GAMMA*(pRef/rRef));
    std::array<real,5> q1,q2,q3;
    for(int j=i-3;j<i+2;j++)
    {
        int iLocal=j-i+3;
        q1[iLocal]=at(j,R)-at(j,P)/(cRef*cRef);
        q2[iLocal]=(at(j,P)+at(j,U)*cRef*rRef)/2;
        q3[iLocal]=(at(j,P)-at(j,U)*cRef*rRef)/2;
    }
    real Q1=(*this->inter5Positive)(q1);
    real Q2=(*this->inter5)(q2);
    real Q3=(*this->inter5)(q3);
    return {Q1+(Q2+Q3)/(cRef*cRef),(Q2-Q3)/cRef/rRef,Q2+Q3};
}
std::vector<real> SpaceDis::reconRChar1D(int i)
{
    enum{R,U,P};
    real rRef=at(i,R);
    real pRef=at(i,P);
    real cRef=sqrt(GAMMA*(pRef/rRef));
    std::array<real,5> q1,q2,q3;
    for(int j=i+2;j>i-3;j--)
    {
        int iLocal=i+2-j;
        q1[iLocal]=at(j,R)-at(j,P)/(cRef*cRef);
        q2[iLocal]=(at(j,P)+at(j,U)*cRef*rRef)/2;
        q3[iLocal]=(at(j,P)-at(j,U)*cRef*rRef)/2;
    }
    real Q1=(*this->inter5Positive)(q1);
    real Q2=(*this->inter5)(q2);
    real Q3=(*this->inter5)(q3);
    return {Q1+(Q2+Q3)/(cRef*cRef),(Q2-Q3)/cRef/rRef,Q2+Q3};
}


std::vector<real> SpaceDis::reconLChar2D(int i)
{
    assert(info->dim==2);
    assert(info->eqType==EULER);
    enum{R,U,V,P};
    real rRef=at(i-1,R);
    real pRef=at(i-1,P);
    real cRef=sqrt(GAMMA*(pRef/rRef));
    // real rRef=(at(i-1,R)+at(i,R))/2;
    // real pRef=(at(i-1,P)+at(i,P))/2;
    // real cRef=sqrt(GAMMA*(pRef/rRef));
    std::array<real,5> q1,q2,q3,q4;
    for(int j=i-3;j<i+2;j++)
    {
        int iLocal=j-i+3;
        real Vn=norm[0]*at(j,U)+norm[1]*at(j,V);
        q1[iLocal]=norm[0]*at(j,V)-norm[1]*at(j,U);
        q2[iLocal]=at(j,R)-at(j,P)/(cRef*cRef);
        q3[iLocal]=(at(j,P)+Vn*cRef*rRef)/2;
        q4[iLocal]=(at(j,P)-Vn*cRef*rRef)/2;
    }
    real Q1=(*this->inter5)(q1);
    real Q2=(*this->inter5Positive)(q2);
    real Q3=(*this->inter5)(q3);
    real Q4=(*this->inter5)(q4);
    std::vector<real> res={Q2+(Q3+Q4)/(cRef*cRef)
            ,-norm[1]*Q1+norm[0]*(Q3-Q4)/cRef/rRef
            ,norm[0]*Q1+norm[1]*(Q3-Q4)/cRef/rRef
            ,Q3+Q4};
    return res;
}

std::vector<real> SpaceDis::reconRChar2D(int i)
{
    assert(info->dim==2);
    enum{R,U,V,P};
    real rRef=at(i,R);
    real pRef=at(i,P);
    real cRef=sqrt(GAMMA*(pRef/rRef));
    std::array<real,5> q1,q2,q3,q4;
    for(int j=i+2;j>i-3;j--)
    {
        int iLocal=i+2-j;
        real Vn=norm[0]*at(j,U)+norm[1]*at(j,V);
        q1[iLocal]=norm[0]*at(j,V)-norm[1]*at(j,U);
        q2[iLocal]=at(j,R)-at(j,P)/(cRef*cRef);
        q3[iLocal]=(at(j,P)+Vn*cRef*rRef)/2;
        q4[iLocal]=(at(j,P)-Vn*cRef*rRef)/2;
    }
    real Q1=(*this->inter5)(q1);
    real Q2=(*this->inter5Positive)(q2);
    real Q3=(*this->inter5)(q3);
    real Q4=(*this->inter5)(q4);

    std::vector<real> res={Q2+(Q3+Q4)/(cRef*cRef)
            ,-norm[1]*Q1+norm[0]*(Q3-Q4)/cRef/rRef
            ,norm[0]*Q1+norm[1]*(Q3-Q4)/cRef/rRef
            ,Q3+Q4};
    return res;
}

std::vector<real> SpaceDis::reconRprim(int i)
{
    std::vector<real> res(nPrim);
    for (int ivar = 0; ivar < nPrim; ivar++)
    {
        std::array iprims={at(i+2,ivar),at(i+1,ivar),at(i,ivar),at(i-1,ivar),at(i-2,ivar)};
        res[ivar]=(*this->inter5)(iprims);
    }
    return res;
}

std::vector<real> SpaceDis::reconLprim(int i)
{
    std::vector<real> res(nPrim);
    for (int ivar = 0; ivar < nPrim; ivar++)
    {
        std::array iprims={at(i-3,ivar),at(i-2,ivar),at(i-1,ivar),at(i,ivar),at(i+1,ivar)};
        res[ivar]=(*this->inter5)(iprims);
    }
    return res;
}

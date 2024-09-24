#include "SpaceDis.hpp"
#include "fluxScheme.hpp"

void SpaceDis::difHCS()
{
    constexpr std::array<real,3> w={64.0/45.0,-2.0/9.0,1.0/180.0};
    real h=info->geth(idim);
    std::vector<real> (*fFunction)(const std::vector<real>& ,std::array<real,3>);
    if(info->eqType==EULER)
    {if(info->dim==1) fFunction=&fEuler1D;
    else if (info->dim==2) fFunction=&fEuler2D;}
    else{fFunction=&fDefault;}
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<nVar;j++)
        {
            (*rhs)(i0+i*offset,j)+=w[0]*(fluxAt(i+1,j)-fluxAt(i,j))/h;
        }
    }
    //i==-2
    std::vector<real> qs(nVar);
    int iNode=-2;
    memcpy(&qs[0],&at(iNode,0),nVar);
    auto fluxNode=fFunction(qs,norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+2)*offset,j)+=-w[2]*fluxNode[j]/h;
    iNode=-1;
    memcpy(&qs[0],&at(iNode,0),nVar);
    fluxNode=fFunction(qs,norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+2)*offset,j)+=-w[2]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+1)*offset,j)+=-w[1]*fluxNode[j]/h;
    iNode=0;
    memcpy(&qs[0],&at(iNode,0),nVar);
    fluxNode=fFunction(qs,norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+2)*offset,j)+=-w[2]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+1)*offset,j)+=-w[1]*fluxNode[j]/h;
    iNode=1;
    memcpy(&qs[0],&at(iNode,0),nVar);
    fluxNode=fFunction(qs,norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+2)*offset,j)+=-w[2]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+1)*offset,j)+=-w[1]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-1)*offset,j)+= w[1]*fluxNode[j]/h;

    for (iNode = 2; iNode < n-2; iNode++)
    {
        memcpy(&qs[0],&at(iNode,0),nVar);
        fluxNode=fFunction(qs,norm);
        for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+2)*offset,j)+=-w[2]*fluxNode[j]/h;
        for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+1)*offset,j)+=-w[1]*fluxNode[j]/h;
        for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-1)*offset,j)+= w[1]*fluxNode[j]/h;
        for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-2)*offset,j)+= w[2]*fluxNode[j]/h;
    }
    

    iNode=n-2;
    memcpy(&qs[0],&at(iNode,0),nVar);
    fluxNode=fFunction(qs,norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+1)*offset,j)+=-w[1]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-1)*offset,j)+= w[1]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-2)*offset,j)+= w[2]*fluxNode[j]/h;
    iNode=n-1;
    memcpy(&qs[0],&at(iNode,0),nVar);
    fluxNode=fFunction(qs,norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-1)*offset,j)+= w[1]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-2)*offset,j)+= w[2]*fluxNode[j]/h;
    iNode=n;
    memcpy(&qs[0],&at(iNode,0),nVar);
    fluxNode=fFunction(qs,norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-1)*offset,j)+= w[1]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-2)*offset,j)+= w[2]*fluxNode[j]/h;
    iNode=n+1;
    memcpy(&qs[0],&at(iNode,0),nVar);
    fluxNode=fFunction(qs,norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-2)*offset,j)+= w[2]*fluxNode[j]/h;
}

void SpaceDis::difMND6()
{
    constexpr std::array<real,3> w={3.0/2.0,-3.0/10.0,1.0/30.0};
    std::vector<real> (*fFunction)(const std::vector<real>& ,std::array<real,3>);
    if(info->eqType==EULER){
    if(info->dim==1) fFunction=&fEuler1D;
    else if (info->dim==2) fFunction=&fEuler2D;}
    else{fFunction=&fDefault;}
    real h=info->geth(idim);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<nVar;j++)
        {
            (*rhs)(i0+i*offset,j)+=w[0]*(fluxAt(i+1,j)-fluxAt(i,j))/h
                                   +w[2]*(fluxAt(i+2,j)-fluxAt(i-1,j))/h;
        }
    }
    //i==-2
    std::vector<real> qs(nVar);
    memcpy(&qs[0],&at(-1,0),nVar*sizeof(real));
    auto fluxNodeM=fFunction(qs,norm);
    memcpy(&qs[0],&at(0,0),nVar*sizeof(real));
    auto fluxNodeR=fFunction(qs,norm);
    for (int iNode = 0; iNode < n; iNode++)
    {
        auto fluxNodeL=fluxNodeM;
        fluxNodeM=fluxNodeR;

        memcpy(&qs[0],&at(iNode+1,0),nVar*sizeof(real));
        fluxNodeR=fFunction(qs,norm);
        for(int j=0;j<nVar;j++)  (*rhs)(i0+iNode*offset,j)+=w[1]*(fluxNodeR[j]-fluxNodeL[j])/h;
    }
}

void SpaceDis::difTraditional6()
{
    real h=info->geth(idim);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<nVar;j++)
        {
            (*rhs)(i0+i*offset,j)+=75.0/64.0*(fluxAt(i+1,j)-fluxAt(i,j))/h
                         -25.0/128.0*(fluxAt(i+2,j)-fluxAt(i-1,j))/(3*h)
                         +3.0/128.0*(fluxAt(i+3,j)-fluxAt(i-2,j))/(5*h);
        }
    }
}

void SpaceDis::dif2Order()
{
    real h=info->geth(idim);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<nVar;j++)
        {
            (*rhs)(i0+i*offset,j)+=(fluxAt(i+1,j)-fluxAt(i,j))/h;
        }
    }
}
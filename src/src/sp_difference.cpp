#include "SpaceDis.hpp"
#include "fluxScheme.hpp"

void SpaceDis::difHCS()
{
    constexpr std::array<real,3> w={64/45,-2/9,1/180};
    real h=info->geth(idim);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<nVar;j++)
        {
            (*rhs)(i0+i*offset,j)+=w[0]*(fluxAt(i+1,j)-fluxAt(i,j))/h;
        }
    }
    //i==-2
    int iNode=-2;
    auto fluxNode=fEuler2D(at(iNode,0),at(iNode,1),at(iNode,2),at(iNode,3),at(iNode,4),norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+2)*offset,j)+=-w[2]*fluxNode[j]/h;
    iNode=-1;
    fluxNode=fEuler2D(at(iNode,0),at(iNode,1),at(iNode,2),at(iNode,3),at(iNode,4),norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+2)*offset,j)+=-w[2]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+1)*offset,j)+=-w[1]*fluxNode[j]/h;
    iNode=0;
    fluxNode=fEuler2D(at(iNode,0),at(iNode,1),at(iNode,2),at(iNode,3),at(iNode,4),norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+2)*offset,j)+=-w[2]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+1)*offset,j)+=-w[1]*fluxNode[j]/h;
    iNode=1;
    fluxNode=fEuler2D(at(iNode,0),at(iNode,1),at(iNode,2),at(iNode,3),at(iNode,4),norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+2)*offset,j)+=-w[2]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+1)*offset,j)+=-w[1]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-1)*offset,j)+= w[1]*fluxNode[j]/h;

    for (iNode = 2; iNode < n-2; iNode++)
    {
        fluxNode=fEuler2D(at(iNode,0),at(iNode,1),at(iNode,2),at(iNode,3),at(iNode,4),norm);
        for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+2)*offset,j)+=-w[2]*fluxNode[j]/h;
        for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+1)*offset,j)+=-w[1]*fluxNode[j]/h;
        for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-1)*offset,j)+= w[1]*fluxNode[j]/h;
        for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-2)*offset,j)+= w[2]*fluxNode[j]/h;
    }
    

    iNode=n-2;
    fluxNode=fEuler2D(at(iNode,0),at(iNode,1),at(iNode,2),at(iNode,3),at(iNode,4),norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode+1)*offset,j)+=-w[1]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-1)*offset,j)+= w[1]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-2)*offset,j)+= w[2]*fluxNode[j]/h;
    iNode=n-1;
    fluxNode=fEuler2D(at(iNode,0),at(iNode,1),at(iNode,2),at(iNode,3),at(iNode,4),norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-1)*offset,j)+= w[1]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-2)*offset,j)+= w[2]*fluxNode[j]/h;
    iNode=n;
    fluxNode=fEuler2D(at(iNode,0),at(iNode,1),at(iNode,2),at(iNode,3),at(iNode,4),norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-1)*offset,j)+= w[1]*fluxNode[j]/h;
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-2)*offset,j)+= w[2]*fluxNode[j]/h;
    iNode=n+1;
    fluxNode=fEuler2D(at(iNode,0),at(iNode,1),at(iNode,2),at(iNode,3),at(iNode,4),norm);
    for(int j=0;j<nVar;j++)  (*rhs)(i0+(iNode-2)*offset,j)+= w[2]*fluxNode[j]/h;
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
            // rhs[i*nVar+j]=(flux(i+1,j)-flux(i,j))/h;
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
#pragma once

#include <array>
#include <vector>

typedef std::array<real,2> arr2;
std::vector<real> roeFlux1D(arr2 r, arr2 u, arr2 p, arr2 H, arr2 RT);
std::vector<real> roeFlux1D2(arr2 r, arr2 u, arr2 p, arr2 H, arr2 RT);
//std::vector<real> HLLCFlux1D(arr2 r, arr2 u, arr2 p, arr2 H, arr2 RT);
constexpr std::vector<real> HLLCFlux1D(arr2 r, arr2 u, arr2 p, arr2 H, arr2 RT)
{
    //reference:https://zhuanlan.zhihu.com/p/583555029
    enum{
    L,
    R
    };
    std::vector<real> res;
    res.resize(3);
    real gamma=GAMMA;
    real cl=sqrt(gamma*RT[L]);
    real cr=sqrt(gamma*RT[R]);

    real SL=std::min(u[L]-cl,u[R]-cr);
    real SR=std::max(u[L]+cl,u[R]+cr);

    real Sstar=(r[L]*u[L]*(SL-u[L])-p[L]
                +p[R]-r[R]*u[R]*(SR-u[R]))/
               (r[L]*(SL-u[L])-r[R]*(SR-u[R]));
    if (SL>=0)
    {
        res[0]=r[L]*u[L];
        res[1]=r[L]*u[L]*u[L]+p[L];
        res[2]=r[L]*H[L]*u[L];
    }
    else if (Sstar>=0)
    {
        real pStar=p[L]+r[L]*(SL-u[L])*(Sstar-u[L]);
        real U[3]={r[L],r[L]*u[L],p[L]/(gamma-1)+r[L]*u[L]*u[L]/2};
        real F[3]={r[L]*u[L],r[L]*u[L]*u[L]+p[L],r[L]*H[L]*u[L]};
        real D[3]={0,1,Sstar};
        for (int i = 0; i < 3; i++)
        {
            res[i]=(Sstar*(SL*U[i]-F[i])+SL*pStar*D[i])/(SL-Sstar);
        }
        
    }
    else if (SR>=0)
    {
        real pStar=p[R]+r[R]*(SR-u[R])*(Sstar-u[R]);
        real U[3]={r[R],r[R]*u[R],p[R]/(gamma-1)+r[R]*u[R]*u[R]/2};
        real F[3]={r[R]*u[R],r[R]*u[R]*u[R]+p[R],r[R]*H[R]*u[R]};
        real D[3]={0,1,Sstar};
        for (int i = 0; i < 3; i++)
        {
            res[i]=(Sstar*(SR*U[i]-F[i])+SR*pStar*D[i])/(SR-Sstar);
        }
    }
    else
    {
        res[0]=r[R]*u[R];
        res[1]=r[R]*u[R]*u[R]+p[R];
        res[2]=r[R]*H[R]*u[R];
    }
    return res;
}

constexpr std::array<real,4> HLLCFlux2D(arr2 r, arr2 u,arr2 v, arr2 p, arr2 H,std::array<real,3> norm)
{
    enum{
    L,
    R
    };
    //reference:https://zhuanlan.zhihu.com/p/583555029
    std::array<real,4> res;
    arr2 Vn={u[L]*norm[0]+v[L]*norm[1],u[R]*norm[0]+v[R]*norm[1]};
    real gamma=GAMMA;
    real q2L=(u[L]*u[L]+v[L]*v[L])/2;
    real q2R=(u[R]*u[R]+v[R]*v[R])/2;
    real cl=sqrt((gamma-1)*(H[L]-q2L));
    real cr=sqrt((gamma-1)*(H[R]-q2R));

    real SL=std::min(Vn[L]-cl,Vn[R]-cr);
    real SR=std::max(Vn[L]+cl,Vn[R]+cr);

    real Sstar=(r[L]*Vn[L]*(SL-Vn[L])-p[L]
                +p[R]-r[R]*Vn[R]*(SR-Vn[R]))/
               (r[L]*(SL-Vn[L])-r[R]*(SR-Vn[R]));
    if (SL>=0)
    {
        res[0]=r[L]         *Vn[L];
        res[1]=r[L]*u[L]*Vn[L]+p[L]*norm[0];
        res[2]=r[L]*v[L]*Vn[L]+p[L]*norm[1];
        res[3]=r[L]*H[L]*Vn[L];
    }
    else if (Sstar>=0)
    {
        real pStar=p[L]+r[L]*(SL-Vn[L])*(Sstar-Vn[L]);
        real U[4]={r[L],
                   r[L]*u[L],
                   r[L]*v[L],
                   r[L]*H[L]-p[L]};
        real F[4]={r[L]     *Vn[L]
                  ,r[L]*u[L]*Vn[L]+p[L]*norm[0]
                  ,r[L]*v[L]*Vn[L]+p[L]*norm[1]
                  ,r[L]*H[L]*Vn[L]};
        real D[4]={0,norm[0],norm[1],Sstar};
        for (int i = 0; i < 4; i++)
        {
            res[i]=(Sstar*(SL*U[i]-F[i])+SL*pStar*D[i])/(SL-Sstar);
        }
        
    }
    else if (SR>=0)
    {
        real pStar=p[R]+r[R]*(SR-Vn[R])*(Sstar-Vn[R]);

        real U[4]={r[R]
                  ,r[R]*u[R]
                  ,r[R]*v[R]
                  ,r[L]*H[R]-p[R]};

        real F[4]={r[R]     *Vn[R]
                  ,r[R]*u[R]*Vn[R]+p[R]*norm[0]
                  ,r[R]*v[R]*Vn[R]+p[R]*norm[1]
                  ,r[R]*H[R]*Vn[R]};

        real D[4]={0,norm[0],norm[1],Sstar};
        for (int i = 0; i < 4; i++)
        {
            res[i]=(Sstar*(SR*U[i]-F[i])+SR*pStar*D[i])/(SR-Sstar);
        }
    }
    else
    {
        res[0]=r[R]     *Vn[R];
        res[1]=r[R]*u[R]*Vn[R]+p[R]*norm[0];
        res[2]=r[R]*v[R]*Vn[R]+p[R]*norm[1];
        res[3]=r[R]*H[R]*Vn[R];
    }
    return res;
}

constexpr std::array<real,4> rhoFlux2D(arr2 r, arr2 u,arr2 v, arr2 p, arr2 H,std::array<real,3> norm)
{
    enum{
    L,
    R
    };
    std::array<real,4> FcL,FcR;
    std::array<real,4> result;
    // arr2 H;
    // H[L]=p[L]/r[L]*GAMMA/(GAMMA-1)+(u[L]*u[L]+v[L]*v[L])/2;
    // H[R]=p[R]/r[R]*GAMMA/(GAMMA-1)+(u[R]*u[R]+v[R]*v[R])/2;
    arr2 Vn={u[L]*norm[0]+v[L]*norm[1],u[R]*norm[0]+v[R]*norm[1]};
    FcL[0]=r[L]*Vn[L];
    FcL[1]=r[L]*u[L]*Vn[L]+p[L]*norm[0];
    FcL[2]=r[L]*v[L]*Vn[L]+p[L]*norm[1];
    FcL[3]=r[L]*H[L]*Vn[L];

    FcR[0]=r[R]*Vn[R];
    FcR[1]=r[R]*u[R]*Vn[R]+p[R]*norm[0];
    FcR[2]=r[R]*v[R]*Vn[R]+p[R]*norm[1];
    FcR[3]=r[R]*H[R]*Vn[R];

    double rhoAvg,uAvg,vAvg,HAvg,cAvg,VnAvg,q_2Avg,coef1,coef2;
    coef1=sqrt(r[L])/(sqrt(r[L])+sqrt(r[R]));
    coef2=1-coef1;
    rhoAvg=sqrt(r[L]*r[R]);
    uAvg=coef1*u[L]+coef2*u[R];
    vAvg=coef1*v[L]+coef2*v[R];
    HAvg=coef1*H[L]+coef2*H[R];
    q_2Avg=(uAvg*uAvg+vAvg*vAvg)*0.5;

    cAvg=sqrt((GAMMA-1)*(HAvg-q_2Avg));
    if((HAvg-q_2Avg)<0)
    {
        cAvg=sqrt(0.1);
    }
    VnAvg=uAvg*norm[0]+vAvg*norm[1];

    double lambda[3]={abs(VnAvg-cAvg),abs(VnAvg),abs(VnAvg+cAvg)};
    double eps=0.2*(abs(uAvg)+abs(vAvg)+cAvg);
    for(int i=0;i<3;i++)
    {
        if(lambda[i]<eps)
        {
            lambda[i]=(lambda[i]*lambda[i]+eps*eps)/(2.0*eps);
        }
    }
    double deltaP=p[R]-p[L],deltaVn=Vn[R]-Vn[L],deltaU=u[R]-u[L],deltaV=v[R]-v[L],
           deltaRho=r[R]-r[L],coef;
    double FDispassion[4];
    coef1=(deltaP-rhoAvg*cAvg*deltaVn)/(2*cAvg*cAvg);
    FDispassion[0]=lambda[0]*coef1*1;
    FDispassion[1]=lambda[0]*coef1*(uAvg-cAvg*norm[0]);
    FDispassion[2]=lambda[0]*coef1*(vAvg-cAvg*norm[1]);
    FDispassion[3]=lambda[0]*coef1*(HAvg-cAvg*VnAvg);

    coef2=deltaRho-deltaP/(cAvg*cAvg);
    FDispassion[0] +=lambda[1]*(coef2*1.0       +  rhoAvg*0.0);
    FDispassion[1]+=lambda[1]*(coef2*uAvg      +  rhoAvg*(deltaU-deltaVn*norm[0]));
    FDispassion[2]+=lambda[1]*(coef2*vAvg      +  rhoAvg*(deltaV-deltaVn*norm[1]));
    FDispassion[3]+=lambda[1]*(coef2*q_2Avg+  rhoAvg*(uAvg*deltaU+vAvg*deltaV-VnAvg*deltaVn));

    coef=(deltaP+rhoAvg*cAvg*deltaVn)/(2.0*cAvg*cAvg);
    FDispassion[0] +=lambda[2]*(coef*1);
    FDispassion[1]+=lambda[2]*(coef*(uAvg+cAvg*norm[0]));
    FDispassion[2]+=lambda[2]*(coef*(vAvg+cAvg*norm[1]));
    FDispassion[3]+=lambda[2]*(coef*(HAvg+cAvg*VnAvg));

    for(int i=0;i<4;i++)
    {
        result[i]=0.5*(FcL[i]+FcR[i]-FDispassion[i]);
    }
    return result;
}

constexpr std::array<real,4> fEuler2D(real r,real u ,real v,real p,real H,std::array<real,3> norm)
{
    real Vn=(u*norm[0]+v*norm[1]);
    std::array<real,4> res={
        r*Vn,
        r*u*Vn,
        r*v*Vn,
        r*H*Vn
    };
    return res;
}
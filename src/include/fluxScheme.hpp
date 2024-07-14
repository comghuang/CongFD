#pragma once

#include <array>
#include <vector>

typedef std::array<real,2> arr2;
std::vector<real> roeFlux1D(real rl,real rr,real ul,real ur,real pl,real pr);
std::vector<real> roeFlux1D2(real rl,real rr,real ul,real ur,real pl,real pr);
//std::vector<real> HLLCFlux1D(arr2 r, arr2 u, arr2 p, arr2 H, arr2 RT);
constexpr std::array<real,3> HLLCFlux1D(real rl,real rr,real ul,real ur,real pl,real pr)
{
    //reference:https://zhuanlan.zhihu.com/p/583555029
    enum{
    L,
    R
    };
    arr2 H={pl/rl*GAMMA/(GAMMA-1)+(ul*ul)/2
           ,pr/rr*GAMMA/(GAMMA-1)+(ur*ur)/2};
    arr2 RT={pl/rl
           ,pr/rr};
    std::array<real,3> res;
    real gamma=GAMMA;
    real cl=sqrt(gamma*RT[L]);
    real cr=sqrt(gamma*RT[R]);

    real SL=std::min(ul-cl,ur-cr);
    real SR=std::max(ul+cl,ur+cr);

    real Sstar=(rl*ul*(SL-ul)-pl
                +pr-rr*ur*(SR-ur))/
               (rl*(SL-ul)-rr*(SR-ur));
    if (SL>=0)
    {
        res[0]=rl*ul;
        res[1]=rl*ul*ul+pl;
        res[2]=rl*H[L]*ul;
    }
    else if (Sstar>=0)
    {
        real pStar=pl+rl*(SL-ul)*(Sstar-ul);
        real U[3]={rl,rl*ul,pl/(gamma-1)+rl*ul*ul/2};
        real F[3]={rl*ul,rl*ul*ul+pl,rl*H[L]*ul};
        real D[3]={0,1,Sstar};
        for (int i = 0; i < 3; i++)
        {
            res[i]=(Sstar*(SL*U[i]-F[i])+SL*pStar*D[i])/(SL-Sstar);
        }
        
    }
    else if (SR>=0)
    {
        real pStar=pr+rr*(SR-ur)*(Sstar-ur);
        real U[3]={rr,rr*ur,pr/(gamma-1)+rr*ur*ur/2};
        real F[3]={rr*ur,rr*ur*ur+pr,rr*H[R]*ur};
        real D[3]={0,1,Sstar};
        for (int i = 0; i < 3; i++)
        {
            res[i]=(Sstar*(SR*U[i]-F[i])+SR*pStar*D[i])/(SR-Sstar);
        }
    }
    else
    {
        res[0]=rr*ur;
        res[1]=rr*ur*ur+pr;
        res[2]=rr*H[R]*ur;
    }
    return res;
}

constexpr std::array<real,4> HLLCFlux2D(real rl,real rr,real ul,real ur,real vl,real vr,real pl,real pr,std::array<real,3> norm)
{
    enum{
    L,
    R
    };
    //reference:https://zhuanlan.zhihu.com/p/583555029
    std::array<real,4> res;
    // if(pl<0 || pr<0 || rr<0 || rl<0)
    // {
    //     std::cout<<"fluxScheme error: positivity break\n";
    // }

     arr2 H;
     H[L]=pl/rl*GAMMA/(GAMMA-1)+(ul*ul+vl*vl)/2;
     H[R]=pr/rr*GAMMA/(GAMMA-1)+(ur*ur+vr*vr)/2;

    arr2 Vn={ul*norm[0]+vl*norm[1],ur*norm[0]+vr*norm[1]};
    real gamma=GAMMA;
    real q2L=(ul*ul+vl*vl)/2;
    real q2R=(ur*ur+vr*vr)/2;
    real cl=sqrt((gamma-1)*(H[L]-q2L));
    real cr=sqrt((gamma-1)*(H[R]-q2R));

    real SL=std::min(Vn[L]-cl,Vn[R]-cr);
    real SR=std::max(Vn[L]+cl,Vn[R]+cr);

    real Sstar=(rl*Vn[L]*(SL-Vn[L])-pl
                +pr-rr*Vn[R]*(SR-Vn[R]))/
               (rl*(SL-Vn[L])-rr*(SR-Vn[R]));
    if (SL>=0)
    {
        res[0]=rl         *Vn[L];
        res[1]=rl*ul*Vn[L]+pl*norm[0];
        res[2]=rl*vl*Vn[L]+pl*norm[1];
        res[3]=rl*H[L]*Vn[L];
    }
    else if (Sstar>=0)
    {
        real pStar=pl+rl*(SL-Vn[L])*(Sstar-Vn[L]);
        real U[4]={rl,
                   rl*ul,
                   rl*vl,
                   rl*H[L]-pl};
        real F[4]={rl     *Vn[L]
                  ,rl*ul*Vn[L]+pl*norm[0]
                  ,rl*vl*Vn[L]+pl*norm[1]
                  ,rl*H[L]*Vn[L]};
        real D[4]={0,norm[0],norm[1],Sstar};
        for (int i = 0; i < 4; i++)
        {
            res[i]=(Sstar*(SL*U[i]-F[i])+SL*pStar*D[i])/(SL-Sstar);
        }
        
    }
    else if (SR>=0)
    {
        real pStar=pr+rr*(SR-Vn[R])*(Sstar-Vn[R]);

        real U[4]={rr
                  ,rr*ur
                  ,rr*vr
                  ,rl*H[R]-pr};

        real F[4]={rr     *Vn[R]
                  ,rr*ur*Vn[R]+pr*norm[0]
                  ,rr*vr*Vn[R]+pr*norm[1]
                  ,rr*H[R]*Vn[R]};

        real D[4]={0,norm[0],norm[1],Sstar};
        for (int i = 0; i < 4; i++)
        {
            res[i]=(Sstar*(SR*U[i]-F[i])+SR*pStar*D[i])/(SR-Sstar);
        }
    }
    else
    {
        res[0]=rr     *Vn[R];
        res[1]=rr*ur*Vn[R]+pr*norm[0];
        res[2]=rr*vr*Vn[R]+pr*norm[1];
        res[3]=rr*H[R]*Vn[R];
    }
    return res;
}

constexpr std::array<real,4> roeFlux2D(real rl,real rr,real ul,real ur,real vl,real vr,real pl,real pr,std::array<real,3> norm)
{
    enum{
    L,
    R
    };
    std::array<real,4> FcL,FcR;
    std::array<real,4> result;
    if(pl<0 || pr<0 || rr<0 || rl<0)
    {
        std::cout<<"fluxScheme error: positivity break\n";
    }


    arr2 H;
    H[L]=pl/rl*GAMMA/(GAMMA-1)+(ul*ul+vl*vl)/2;
    H[R]=pr/rr*GAMMA/(GAMMA-1)+(ur*ur+vr*vr)/2;
    
    arr2 Vn={ul*norm[0]+vl*norm[1],ur*norm[0]+vr*norm[1]};
    FcL[0]=rl*Vn[L];
    FcL[1]=rl*ul*Vn[L]+pl*norm[0];
    FcL[2]=rl*vl*Vn[L]+pl*norm[1];
    FcL[3]=rl*H[L]*Vn[L];

    FcR[0]=rr*Vn[R];
    FcR[1]=rr*ur*Vn[R]+pr*norm[0];
    FcR[2]=rr*vr*Vn[R]+pr*norm[1];
    FcR[3]=rr*H[R]*Vn[R];

    double rhoAvg,uAvg,vAvg,HAvg,cAvg,VnAvg,q_2Avg,coef1,coef2;
    coef1=sqrt(rl)/(sqrt(rl)+sqrt(rr));
    coef2=1-coef1;
    rhoAvg=sqrt(rl*rr);
    uAvg=coef1*ul+coef2*ur;
    vAvg=coef1*vl+coef2*vr;
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
    double deltaP=pr-pl,deltaVn=Vn[R]-Vn[L],deltaU=ur-ul,deltaV=vr-vl,
           deltaRho=rr-rl,coef;
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

constexpr std::vector<real> fEuler2D(const std::vector<real>& q ,std::array<real,3> norm)
{
    real r=q[0],u=q[1],v=q[2],p=q[3];
    real Vn=(u*norm[0]+v*norm[1]);
    //real H=GAMMA/(GAMMA-1)*p/r+(u*u+v*v)/2;
    std::vector<real> res={
        r*Vn,
        r*u*Vn+p*norm[0],
        r*v*Vn+p*norm[1],
        (GAMMA/(GAMMA-1)*p+r*(u*u+v*v)/2)*Vn
    };
    return res;
}

constexpr std::vector<real> fEuler1D(const std::vector<real>& q ,std::array<real,3> norm)
{
    real r=q[0],u=q[1],p=q[2];
    //real H=GAMMA/(GAMMA-1)*p/r+(u*u)/2;
    std::vector<real> res={
        r*u,
        r*u*u+p,
        (GAMMA/(GAMMA-1)*p+r*(u*u)/2)*u
    };
    return res;
}
constexpr std::vector<real> fDefault(const std::vector<real>& q ,std::array<real,3> norm)
{
    return {q[0]};
}
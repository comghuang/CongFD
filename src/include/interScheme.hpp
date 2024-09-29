#pragma once
#include <array>
#include "macro.hpp"
constexpr real weno5_JSchen(std::array<real,5>);
constexpr real u1(real q1,real q2 ,real q3)
{
    return 3.0/8.0*q1-5.0/4.0*q2+15.0/8.0*q3;
}
constexpr real u2(real q1,real q2 ,real q3)
{
    return -1.0/8.0*q1+3.0/4.0*q2+3.0 /8.0*q3;
}
constexpr real u3(real q1,real q2 ,real q3)
{
    return 3.0/8.0*q1+3.0/4.0*q2-1.0 /8.0*q3;
}

constexpr real weno5_JSchen(std::array<real,5> q)
{
    real eps=1e-6;
    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);


    std::array<real,3> u;
    u[0]= 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
    u[1]=-1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
    u[2]= 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
    
    //std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};
    
    real sumbeta=0,result=0;
    for(int i=0;i<3;i++)
    {
        beta[i]=gamma[i]/pow(eps+beta[i],2);
        sumbeta+=beta[i];
    }
    for(int i=0;i<3;i++) result+=beta[i]*u[i];
    return result/sumbeta;
}

constexpr real weno5_Cong(std::array<real,5> q)
{
    real eps=1e-14;
    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta,u;
    // beta[0]= 0.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
    //          + 1.0/1.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    //  beta[1]= 0.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
    //          + 1.0/1.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    //  beta[2]= 0.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
    //          + 1.0/1.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);
     
    beta[0]=1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)/(q[2]*q[2])+10 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2)/(q[2]*q[2]);

    beta[1]=1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)/(q[2]*q[2])+10 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2)/(q[2]*q[2]);

    beta[2]=1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)/(q[2]*q[2])+10*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2)/(q[2]*q[2]);
    
    u[0]= 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
    u[1]=-1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
    u[2]= 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
    
    real sumbeta=0,result=0;
    for(int i=0;i<3;i++)
    {
        beta[i]=gamma[i]/pow(eps+beta[i],2.0);
        sumbeta+=beta[i];
    }
    for(int i=0;i<3;i++) result+=beta[i]*u[i];
    return result/sumbeta;
}
constexpr real weno5_Z(std::array<real,5> q)
{
    real eps=1e-40;
    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);

    std::array<real,3> u;
    u[0]= 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
    u[1]=-1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
    u[2]= 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
    
    //std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};
    
    real sumbeta=0,result=0;
    real C=1,qq=2,tau=abs(beta[2]-beta[0]);
    for(int i=0;i<3;i++)
    {
        beta[i]=gamma[i]*(C+pow(tau/(beta[i]+eps),qq));
        sumbeta+=beta[i];
    }
    for(int i=0;i<3;i++) result+=beta[i]*u[i];
    return result/sumbeta;
}

constexpr real Teno5_Z(std::array<real,5> q)
{
    real eps=1e-40;
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);

    real sumbeta=0;
    real C=1,qq=6,tau=abs(beta[2]-beta[0]);
    for(int i=0;i<3;i++)
    {
        real tempp=C+tau/(beta[i]+eps);
        tempp*=tempp;
        beta[i]=tempp*tempp*tempp;
        sumbeta+=beta[i];
    }
    real CT=1e-5*sumbeta;
    // if(isinf(sumbeta))
    // {
    //     std::cout<<"inf beta\n";
    // }
    // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);
    unsigned short flag=0;
    if(beta[0]<CT) flag+=1;
    if(beta[1]<CT) flag+=2;
    if(beta[2]<CT) flag+=4;
    switch (flag)
    {
    case 0:
        /* 1,1,1 */
        return 3.0/128.0*q[0]-5.0/32.0*q[1]+45.0/64.0*q[2]+15.0 /32.0*q[3]-5.0/128.0*q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0/16.0*q[1]+9.0/16.0*q[2]+9.0 /16.0*q[3]-1.0/16.0*q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0/16.0*q[0]-5.0/16.0*q[1]+15.0/16.0*q[2]+5.0/16.0*q[3];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}
constexpr real Teno5_Z2(std::array<real,5> q)
{
    real eps=1e-40;
    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);
    
    std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};
    
    real sumbeta=0,result=0;
    real C=1,qq=6,tau=abs(beta[2]-beta[0]);
    for(int i=0;i<3;i++)
    {
        real tempp=C+tau/(beta[i]+eps);
        tempp*=tempp;
        beta[i]=tempp*tempp*tempp;
        sumbeta+=beta[i];
    }


    real CT=1e-5*sumbeta,sumGamma=0;
    for(int i=0;i<3;i++) 
    {
        if(beta[i]>CT)
        {
            result+=gamma[i]*(*u[i])(q[i],q[i+1],q[i+2]);
            sumGamma+=gamma[i];
        }
    }
    result/=sumGamma;
    return result;
}


const static real CTi=pow(1.5*1e-5,1.0/6.0);
constexpr real Teno5_CongZ(std::array<real,5> q)
{
    real eps=1e-40;//1e-10;
    std::array<real,3> beta={
    1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2),

    1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2),

    1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2)};

    //int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1: 2):((beta[2]>beta[0])? 0 : 2);
    unsigned short minBeta=std::min_element(beta.begin(),beta.end())-beta.begin();
    constexpr real CT=0.15704178024750198;
    constexpr real CT_1=1-CT;
    real tau=abs(beta[2]-beta[0]);//,KK=0.15704178024750198*(beta[minBeta]+tau);
    real rr=CT*tau-CT_1*beta[minBeta];
    real ll=tau*beta[minBeta];
    //unsigned flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
    unsigned short flag=0;
    if(ll<rr*beta[0]) flag+=1;
    if(ll<rr*beta[1]) flag+=2;
    if(ll<rr*beta[2]) flag+=4;
    switch (flag)
    {
    case 0:
        /* 1,1,1 */
        return 3.0/128.0*q[0]-5.0/32.0*q[1]+45.0/64.0*q[2]+15.0 /32.0*q[3]-5.0/128.0*q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0/16.0*q[1]+9.0/16.0*q[2]+9.0 /16.0*q[3]-1.0/16.0*q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0/16.0*q[0]-5.0/16.0*q[1]+15.0/16.0*q[2]+5.0/16.0*q[3];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}
constexpr real Teno5_CongZ2(std::array<real,5> q)
{
    real eps=1e-40;//1e-10;
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);

    int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1: 2):((beta[2]>beta[0])? 0 : 2);
    real tau=abs(beta[2]-beta[0]),KK=CTi*(beta[minBeta]+tau);
    real sumGamma=0,result=0;
    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};
    for(int i=0;i<3;i++) 
    {
        if(minBeta==i||(beta[i]+tau)*beta[minBeta]>KK*beta[i])
        {
            sumGamma+=gamma[i];
            result+=gamma[i]*(*u[i])(q[i],q[i+1],q[i+2]);
        }
    }
    result/=sumGamma;
    return result;
}

constexpr real Teno5_Cong(std::array<real,5> q)
{
    real eps=1e-40;
    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);
    
    real sumbeta=0,result=0;

    real minBeta=*std::min_element(beta.begin(),beta.end());

    std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};
    real CT=3.5,sumGamma=0;
    for(int i=0;i<3;i++) 
    {
        if(beta[i]<CT*(minBeta+eps))
        {
            sumGamma+=gamma[i];
            result+=gamma[i]*(*u[i])(q[i],q[i+1],q[i+2]);
        }
    }
    result/=sumGamma;
    return result;
}

constexpr real Teno5_Cong2(std::array<real,5> q)
{
    real eps=1e-40;
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);
    
    real sumbeta=0,result=0;

    real minBeta=*std::min_element(beta.begin(),beta.end());
    if (minBeta<eps) return q[2];
    real CT=2,sumGamma=0;
    int flag=0,flags[3]={1,2,4};
    
    for(int i=0;i<3;i++) 
    {
        if(beta[i]<CT*(minBeta+eps))
        {
            flag+=flags[i];
        }
    }
    switch (flag)
    {
    case 1:
        /* 1,0,0 */
        result=3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
        break;
    case 2:
        /* 0,1,0 */
        result=-1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
        break;
    case 3:
        /* 1,1,0 */
        result=1.0/16.0*q[0]-5.0/16.0*q[1]+15.0/16.0*q[2]+5.0/16.0*q[3];
        break;
    case 4:
    case 5:
        /* 0,0,1 */
        /* 1,0,1 */
        result=3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
        break;
    case 6:
        /* 0,1,1 */
        result=-1.0/16.0*q[1]+9.0/16.0*q[2]+9.0 /16.0*q[3]-1.0/16.0*q[4];
        break;
    case 7:
        /* 1,1,1 */
        result=3.0/128.0*q[0]-5.0/32.0*q[1]+45.0/64.0*q[2]+15.0 /32.0*q[3]-5.0/128.0*q[4];
        break;
    default:
        /* 0,0,0 */
        result=q[2];
        break;
    }
    
    return result;
}


constexpr std::array<real,3> Teno5_BVDCong(std::array<real,5> q,bool &flag)
{
    real eps=1e-40;
    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);
    
    std::array<real,3> result={0,0}; 
    real sumbeta=0;
    real C=1,qq=6,tau=abs(beta[2]-beta[0]);
    for(int i=0;i<3;i++)
    {
        beta[i]=(C+pow(tau/(beta[i]+eps),qq));
        sumbeta+=beta[i];
    }
    for (int i = 0; i < 3; i++) beta[i]/=sumbeta;


    std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};
    real CT1=1e-7,CT2=1e-3,sumGamma=0;
    for(int i=0;i<3;i++) 
    {
        if(beta[i]>CT1)
        {
            sumGamma+=gamma[i];
            result[0]+=gamma[i]*(*u[i])(q[i],q[i+1],q[i+2]);
        }
        
    }
    result[0]/=sumGamma;
    
    real extremPoint=0.5*(q[1]-q[3])/(q[1]-2*q[2]+q[3]+eps);
    if(!(extremPoint>=1 || extremPoint<=0)) 
    {result[1]=q[2];
    result[2]=result[0];
    flag=false;}
    else
    {
        flag=true;
        sumGamma=0;
        for(int i=0;i<3;i++) 
        {
            if(beta[i]>CT2)
            {
                sumGamma+=gamma[i];
                result[1]+=gamma[i]*(*u[i])(q[i],q[i+1],q[i+2]);
            }
        }
        result[1]/=sumGamma;
        // auto minBetai=std::max_element(beta.begin(),beta.end())-beta.begin();
        // result[1]=(*u[minBetai])(q[minBetai],q[minBetai+1],q[minBetai+2]);
        // {result[2]=q[2];}
        
    }
    
    
    return result;
}

constexpr std::array<real,3> Teno5_BVDMR(std::array<real,5> q,bool &flag)
{
    real eps=1e-6;
    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);
    
    std::array<real,3> result={0,0}; 
    real sumbeta=0;

    auto minBetap=std::min_element(beta.begin(),beta.end());
    auto minBetai=minBetap-beta.begin();
    auto minBeta=*minBetap;
    if(minBeta<eps) return {q[2],q[2],q[2]};

    real CT1=20,CT2=1.5,sumGamma=0;
    int flag1=0,flag2=0,flags[]={1,2,4};
    for(int i=0;i<3;i++) 
    {
        if(beta[i]<CT1*(minBeta+eps))
        {
            flag1+=flags[i];
        }
        if(beta[i]<CT2*(minBeta+eps))//(extremPoint>=2-minBetai || extremPoint<=1-minBetai))
        {
            flag2+=flags[i];
        }
    }
    switch (flag1)
    {
    case 1:
        /* 1,0,0 */
        result[0]=3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
        break;
    case 2:
        /* 0,1,0 */
        result[0]=-1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
        break;
    case 3:
        /* 1,1,0 */
        result[0]=1.0/16.0*q[0]-5.0/16.0*q[1]+15.0/16.0*q[2]+5.0/16.0*q[3];
        break;
    case 4:
    case 5:
        /* 0,0,1 */
        /* 1,0,1 */
        result[0]=3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
        break;
    case 6:
        /* 0,1,1 */
        result[0]=-1.0/16.0*q[1]+9.0/16.0*q[2]+9.0 /16.0*q[3]-1.0/16.0*q[4];
        break;
    case 7:
        /* 1,1,1 */
        result[0]=3.0/128.0*q[0]-5.0/32.0*q[1]+45.0/64.0*q[2]+15.0 /32.0*q[3]-5.0/128.0*q[4];
        break;
    default:
        /* 0,0,0 */
        result[0]=q[2];
        break;
    }
    switch (flag2)
    {
    case 1:
        /* 1,0,0 */
        result[1]=-5.0/12.0*q[0]+1.0/3.0*q[1]+13.0/12.0*q[2];
        break;
    case 2:
        /* 0,1,0 */
        result[1]=1.0/12.0*q[1]+1.0/3.0*q[2]+7.0/12.0*q[3];
        break;
    case 3:
        /* 1,1,0 */
        result[1]=-9.0/80.0*q[0]-17.0/80.0*q[1]+33.0/80.0*q[2]+39.0/80.0*q[3];
        break;
    case 4:
    case 5:
        /* 0,0,1 */
        /* 1,0,1 */
        result[1]=7.0/12.0*q[2]+1.0/3.0*q[3]+1.0/12.0*q[4];
        break;
    case 6:
        /* 0,1,1 */
        result[1]=-1.0/16.0*q[1]+9.0/16.0*q[2]+9.0 /16.0*q[3]-1.0/16.0*q[4];
        break;
    case 7:
        /* 1,1,1 */
        result[1]=-3.0/160.0*q[0]+1.0/80.0*q[1]+9.0/20.0*q[2]+51.0 /80.0*q[3]-13.0/160.0*q[4];
        break;
    default:
        /* 0,0,0 */
        result[1]=q[2];
        break;
    }
    result[1]=result[0];
    
    return result;
}

constexpr real Teno5_CongSort(std::array<real,5> q)
{
    real eps=1e-12;
    std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};

    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);

    real CC=(beta[1]*2-beta[2]-beta[0])/2;
    // beta[0]+=2.0/3.0*CC;
    // beta[1]-=1.0/3.0*CC;
    // beta[2]+=2.0/3.0*CC;
    
    // 排序
    std::array<real,3> index={0,1,2};
    std::sort(index.begin(), index.end(),
        [&](const int& a, const int& b) {
            return (beta[a] < beta[b]);
        }
    );
    //if(beta[index[2]]<eps) return q[2];
    int ii=index[0];
    real result=(*u[ii])(q[ii],q[ii+1],q[ii+2])*gamma[index[0]];

    real extremPoint=0.5*(q[1]-q[3])/(q[1]-2*q[2]+q[3]+eps);
    real extremPoint2=0.5+(q[2]-q[3])/(q[2]-2*q[3]+q[4]+eps);

    bool flags[3]={true//(u1<*std::max_element(q.begin(),q.end())&&(u1>*std::min_element(q.begin(),q.end())))
            ,(extremPoint>=1 || extremPoint<=0)
            ,(extremPoint2>=1 || extremPoint2<=0)};//;|| (factor>6);
    // bool flag=flags[1]||flags[2];//;|| (factor>6);
    real CT2,CT=3,sumGamma=gamma[index[0]];
    // if(!flag) 
    // {
    //     CT2=3;
    //     CT=CT2;//-beta[index[1]]/(beta[index[0]]+eps)/1.5;
    // }
    // else 
    // {
    //     CT2=3;
    //     CT=CT2;//-beta[index[1]]/(beta[index[0]]+eps)/1.5;
    // }
    //return result/sumGamma;
    
    ii=index[1];
    if (beta[ii]<CT*(beta[index[0]]+eps))
    {
        sumGamma+=gamma[ii];
        result+=gamma[ii]*(*u[ii])(q[ii],q[ii+1],q[ii+2]);
    }
    else 
    {
        //real extremPoint=0.5*(q[1]-q[3])/(q[1]-2*q[2]+q[3]+eps);
        // int i0=index[0];
        // if(!flags[i0]) return q[2];
        return result/sumGamma;
    }
    // {if(flag) return result/sumGamma;
    //  else 
    //  return q[2];}

    ii=index[2];
    if (beta[ii]<(CT)*(beta[index[0]]+eps))
    {
        sumGamma+=gamma[ii];
        result+=gamma[ii]*(*u[ii])(q[ii],q[ii+1],q[ii+2]);
    }
    result/=sumGamma;
    
    return result;
}

constexpr real Teno5_CongSortPositive(std::array<real,5> q)
{
    real eps=1e-12;
    std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};

    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);
    
    // 排序
    std::array<real,3> index={0,1,2};
    std::sort(index.begin(), index.end(),
        [&](const int& a, const int& b) {
            return (beta[a] < beta[b]);
        }
    );
    //if(beta[index[2]]<eps) return q[2];
    int ii=index[0];
    real result=(*u[ii])(q[ii],q[ii+1],q[ii+2])*gamma[index[0]];
    if (result<0) return q[2];

    

    real extremPoint=0.5*(q[1]-q[3])/(q[1]-2*q[2]+q[3]+eps);
    real extremPoint2=0.5+(q[2]-q[3])/(q[2]-2*q[3]+q[4]+eps);
    real factor=abs((q[2]-q[1]+eps)/(q[2]-q[3]+eps));
    bool flag=(extremPoint>=1+eps && extremPoint>=0-eps)||(extremPoint2<=1+eps && extremPoint2>=0-eps);//;|| (factor>6);
    real CT2=factor,CT,sumGamma=gamma[index[0]];
    if(flag) 
    {
        CT2=4;
        CT=CT2-beta[index[1]]/(beta[index[0]]+eps);
    }
    else 
    {
        CT2=20;
        CT=CT2-beta[index[1]]/(beta[index[0]]+eps)/1.5;
    }
    
    //return result/sumGamma;

    
    ii=index[1];
    real u1=(*u[ii])(q[ii],q[ii+1],q[ii+2]);
    if (beta[ii]<CT*(beta[index[0]]) && u1>0)
    {
        sumGamma+=gamma[ii];
        result+=gamma[ii]*u1;
    }
    else
    {if(!flag) return result/sumGamma;
     else return q[2];}
    

    ii=index[2];
    real u2=(*u[ii])(q[ii],q[ii+1],q[ii+2]);
    if (beta[ii]<(CT)*(beta[index[0]]+eps)&& u1>0)
    {
        sumGamma+=gamma[ii];
        result+=gamma[ii]*(*u[ii])(q[ii],q[ii+1],q[ii+2]);
    }
    result/=sumGamma;
    
    return result;
}

constexpr real Teno5_CongSortabs(std::array<real,5> q)
{
    real eps=1e-12;
    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);
    
    // 排序
    std::array<real,3> index={0,1,2};
    std::sort(index.begin(), index.end(),
        [&](const int& a, const int& b) {
            return (beta[a] < beta[b]);
        }
    );

    std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};


    real factor=abs((q[2]-q[1]+eps)/(q[2]-q[3]+eps));
    real CT2=factor,CT,sumGamma=gamma[index[0]];
    if(factor>6) 
    {
        CT2=6;
        CT=CT2-beta[index[1]]/(beta[index[0]]+eps);
    }
    else 
    {
        CT2=20;
        CT=CT2-beta[index[1]]/(beta[index[0]]+eps)/1.5;
    }
    int ii=index[0];
    real result=(*u[ii])(q[ii],q[ii+1],q[ii+2])*gamma[index[0]];
    //return result/sumGamma;

    if(beta[index[2]]<eps) return q[2];
    //real critical=pow(q[2]-q[1],2);
    //CT2=std::max(20.0*pow(critical,2.0),4.0);

    
    //CT=4;
    ii=index[1];
    if (beta[ii]<CT*(beta[index[0]]))
    {
        sumGamma+=gamma[ii];
        result+=gamma[ii]*(*u[ii])(q[ii],q[ii+1],q[ii+2]);
    }
    else
    {if(factor<10) return result/sumGamma;
     else return q[2];}
    

    ii=index[2];
    if (beta[ii]<(CT)*(beta[index[0]]+eps))
    {
        sumGamma+=gamma[ii];
        result+=gamma[ii]*(*u[ii])(q[ii],q[ii+1],q[ii+2]);
    }
    result/=sumGamma;
    
    return result;
}

constexpr real Teno5_CongIncrease(std::array<real,5> q)
{
    real eps=1e-20;
    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);
    
    // 排序
    std::array<real,3> index={0,1,2};
    std::sort(index.begin(), index.end(),
        [&](const int& a, const int& b) {
            return (beta[a] < beta[b]);
        }
    );

    //过于光滑
     if(  beta[index[2]]<=1e-10 ) return q[2];

    std::array<real,3> u;
    u[0]= 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
    u[1]=-1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
    u[2]= 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
    real u0= 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
    real sumbeta=beta[index[0]],sumGamma=gamma[index[0]];
    real result=u[index[0]]*gamma[index[0]];

    int ii=index[1];
    real CT=10,a1=(beta[index[0]]+eps);
    if (beta[ii]<CT*a1) 
    {
        sumGamma+=gamma[ii];
        result+=gamma[ii]*u[ii];
        sumbeta+=beta[ii];
    }
    else 
    {
        return result/sumGamma;
    }

    real CT2=beta[ii];
    ii=index[2];
    real ddd=3;
    if (pow(beta[ii]*a1,ddd)-5*pow(CT2*CT2,ddd)<0) 
    {
        sumGamma+=gamma[ii];
        result+=gamma[ii]*u[ii];
        sumbeta+=beta[ii];
    }
    else
    {
        return result/sumGamma;
    }

    return result/sumGamma;
}

constexpr real musclInterpolation(real q1,real q2,real q3)
{
    
    real delta;
    real deltam,deltap;
    deltam=q2-q1;
    deltap=q3-q2;
    
    //minmod
    real beta=1.0;
    if (deltap>0)
    {
        delta=std::max(0.0,std::max(std::min(beta*deltam,deltap),std::min(deltam,beta*deltap)));
    }
    else
    {
        delta=std::min(0.0,std::min(std::max(beta*deltam,deltap),std::max(deltam,beta*deltap)));
    }
    return q2+delta*0.5;
}
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
    std::array<real,3> beta,u;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);
    
    u[0]= 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
    u[1]=-1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
    u[2]= 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
    
    real sumbeta=0,result=0;
    for(int i=0;i<3;i++)
    {
        beta[i]=gamma[i]/pow(eps+beta[i],2.0);
        sumbeta+=beta[i];
    }
    for(int i=0;i<3;i++) result+=beta[i]/sumbeta*u[i];
    return result;
}

constexpr real weno5_Cong(std::array<real,5> q)
{
    real eps=1e-40;
    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta,u;
//     beta[0]= 0.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
//              + 1.0/1.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

//      beta[1]= 0.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
//              + 1.0/1.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

//      beta[2]= 0.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
//              + 1.0/1.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);
     
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
    for(int i=0;i<3;i++) result+=beta[i]/sumbeta*u[i];
    return result;
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
    
    std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};
    
    real sumbeta=0,result=0;
    real C=1,qq=2,tau=abs(beta[2]-beta[0]);
    for(int i=0;i<3;i++)
    {
        beta[i]=gamma[i]*(C+pow(tau/(beta[i]+eps),qq));
        sumbeta+=beta[i];
    }
    for(int i=0;i<3;i++) result+=beta[i]/sumbeta*(*u[i])(q[i],q[i+1],q[i+2]);
    return result;
}

constexpr real Teno5_Z(std::array<real,5> q)
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
    
    //std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};
    std::array<real,3> u;
    u[0]= 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
    u[1]=-1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
    u[2]= 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
    
    real sumbeta=0,result=0;
    real C=1,qq=6,tau=abs(beta[2]-beta[0]);
    for(int i=0;i<3;i++)
    {
        beta[i]=(C+pow(tau/(beta[i]+eps),qq));
        sumbeta+=beta[i];
    }
    for (int i = 0; i < 3; i++) beta[i]/=sumbeta;
    
    real CT=1e-5,sumGamma=0;
    for(int i=0;i<3;i++) 
    {
        if(beta[i]>CT)
        {
            sumGamma+=gamma[i];
            //result+=gamma[i]*(*u[i])(q[i],q[i+1],q[i+2]);
            result+=gamma[i]*u[i];
        }
    }
    result/=sumGamma;
    return result;
}

constexpr real Teno5_Cong(std::array<real,5> q)
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
    
    real sumbeta=0,result=0;

    real minBeta=*std::min_element(beta.begin(),beta.end());

    std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};
    real CT=12.0,sumGamma=0;
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

constexpr real Teno5_CongSort(std::array<real,5> q)
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
    //for(int i=0;i<3;i++) beta[i]=beta[i]/(q[2]*q[2]+eps);
    //过于光滑
    //if(  beta[index[2]]<=eps ) return q[2];

    std::array<real,3> u;
    u[0]= 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
    u[1]=-1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
    u[2]= 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
    real sumbeta=beta[index[0]],CT=9,sumGamma=gamma[index[0]];
    real result=u[index[0]]*gamma[index[0]];
    bool corner=true;
    for(int i=1;i<3;i++) 
    {
        int ii=index[i];
        if(beta[ii]<CT*(beta[index[0]]+eps))
        {
            sumGamma+=gamma[ii];
            result+=gamma[ii]*u[ii];
            sumbeta+=beta[ii];
            CT=15.0-5.0/3.0*beta[ii]/beta[ii-1]; //CT=CT/2;
            if(i==3) corner=false;
        }
        else{break;}
    }
    result/=sumGamma;
    //int icentral=index[0]+1;
    //std::array<real,3> qs={q[icentral+1],q[icentral],q[icentral-1]};
    if(0) //|| pow((q[2]-qs[3])/(qs[2]-qs[1]),2)>2)
    {
        int icentral=index[0]+1;
        std::array<real,3> qs={q[icentral+1],q[icentral],q[icentral-1]};
        real critical=beta[index[0]]/pow(q[2]-q[3],2);//pow((qs[0]-qs[1])/(qs[1]-qs[2]),2);
        if(critical>10 ) return (q[2]+q[3])*0.5;//|| critical<1.0/20.0
    }
    
    return result;
}

constexpr real Teno5_CongSortPositive(std::array<real,5> q)
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
    //if(  beta[index[2]]<=eps ) return q[2];

    std::array<real,3> u;
    u[0]= 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
    u[1]=-1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
    u[2]= 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
    real sumbeta=beta[index[0]],CT=15,sumGamma=gamma[index[0]];
    real result=u[index[0]]*gamma[index[0]];
    for(int i=1;i<3;i++) 
    {
        int ii=index[i];
        if(beta[ii]<CT*(sumbeta+eps) && u[ii]>0)
        {
            sumGamma+=gamma[ii];
            result+=gamma[ii]*u[ii];
            sumbeta+=beta[ii];
            CT=5; //CT=CT/2;
        }
        else{break;}
    }
    if (result<0) return q[2];
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
    if(factor>1) factor=1.0/factor;
    real CT2=std::max(6.0*factor,2.0),CT,sumGamma=gamma[index[0]];
    int ii=index[0];
    real result=(*u[ii])(q[ii],q[ii+1],q[ii+2])*gamma[index[0]];
    //return result/sumGamma;

    if(beta[index[2]]<eps) return q[2];
    //real critical=pow(q[2]-q[1],2);
    //CT2=std::max(20.0*pow(critical,2.0),4.0);

    CT=CT2-beta[index[1]]/(beta[index[0]]+eps);
    //CT=4;
    ii=index[1];
    if (beta[ii]<CT*(beta[index[0]]))
    {
        sumGamma+=gamma[ii];
        result+=gamma[ii]*(*u[ii])(q[ii],q[ii+1],q[ii+2]);
    }
    else{return result/sumGamma;}
    

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
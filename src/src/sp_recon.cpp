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
    auto start = std::chrono::high_resolution_clock::now(); 
    real Q1=(*this->inter5Positive)(q1);
    real Q2=(*this->inter5)(q2);
    real Q3=(*this->inter5)(q3);
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count(); 
    timep+=duration;
    return {Q1+(Q2+Q3)/(cRef*cRef),(Q2-Q3)/cRef/rRef,Q2+Q3};
}
std::vector<real> SpaceDis::reconRChar1D(int i)
{
    enum{R,U,P};
    real rRef=at(i-1,R);
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
    auto start = std::chrono::high_resolution_clock::now(); 
    real Q1=(*this->inter5Positive)(q1);
    real Q2=(*this->inter5)(q2);
    real Q3=(*this->inter5)(q3);
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count(); 
    timep+=duration;
    return {Q1+(Q2+Q3)/(cRef*cRef),(Q2-Q3)/cRef/rRef,Q2+Q3};
}

static std::array<real,2> minDif(real u1L,real u2L,real u1R,real u2R)
{
    std::array<real,4> difs={abs(u1L-u1R),abs(u1L-u2R),abs(u2L-u1L),abs(u2L-u2R)};
    int index=std::min_element(difs.begin(),difs.end())-difs.begin();
    if(index==0) return{u1L,u1R};
    if(index==1) return{u1L,u2R};
    if(index==2) return{u2L,u1R};
    if(index==3) return{u2L,u2R};
    std::cout<<"spRecon Error: minDif index out\n";
    return{0,0};
}

static std::array<real,2> minDif(std::array<real,3> uL,std::array<real,3> uR)
{
    std::array<real,2> difs={abs(uL[0]-uR[0]),abs(uL[1]-uR[1])};
    int index=std::min_element(difs.begin(),difs.end())-difs.begin();
         if(index==0) return{uL[0],uR[0]};
    else if(index==1) 
    return{uL[1],uR[1]};
    std::cout<<"spRecon Error: minDif index out\n";
    return{0,0};
}
// static std::array<real,2> minDif(std::array<real,3> uL,std::array<real,3> uR)
// {
//     std::array<real,4> difs={abs(uL[0]-uR[0]),abs(uL[1]-uR[1]),abs(uL[0]-uR[1]),abs(uL[1]-uR[0])};
//     int index=std::min_element(difs.begin(),difs.end())-difs.begin();
//          if(index==0) return{uL[0],uR[0]};
//     else if(index==1) 
//     return{uL[1],uR[1]};
//     else if(index==2) 
//     return{uL[0],uR[1]};
//     else if(index==3) 
//     return{uL[1],uR[0]};
//     std::cout<<"spRecon Error: minDif index out\n";
//     return{0,0};
// }
static std::array<real,2> minDif2(std::array<real,3> uL,std::array<real,3> uR)
{
    std::array<real,3> difs={abs(uL[0]-uR[0]),abs(uL[0]-uR[1]),abs(uL[1]-uR[0])};
    int index=std::min_element(difs.begin(),difs.end())-difs.begin();
         if(index==0) return{uL[0],uR[0]};
    else if(index==1) 
    return{uL[0],uR[1]};
    else if(index==2) 
    return{uL[1],uR[0]};
    std::cout<<"spRecon Error: minDif index out\n";
    return{0,0};
}

std::vector<real> SpaceDis::recon1DBVD(int i)
{
    assert(info->dim==1);
    assert(info->eqType==EULER);
    enum{R,U,P};
    real rRef=(at(i-1,R)+at(i,R))/2;
    real pRef=(at(i-1,P)+at(i,P))/2;
    real cRef=sqrt(GAMMA*(pRef/rRef));
    if(isnan(cRef))
    std::cout<<"cRef error nan\n";
    std::array<real,5> q1L,q2L,q3L,q1R,q2R,q3R;
    for(int j=i-3;j<i+3;j++)
    {
        real q1tmp=at(j,R)-at(j,P)/(cRef*cRef);
        real q2tmp=(at(j,P)+at(j,U)*cRef*rRef)/2;
        real q3tmp=(at(j,P)-at(j,U)*cRef*rRef)/2;

        int iLocal=j-i+3;
        if(iLocal<5){
        q1L[iLocal]=q1tmp;
        q2L[iLocal]=q2tmp;
        q3L[iLocal]=q3tmp;}

        iLocal=i+2-j;
        if(iLocal<5){
        q1R[iLocal]=q1tmp;
        q2R[iLocal]=q2tmp;
        q3R[iLocal]=q3tmp;}
    }
    bool flagL1,flagR1,flagL2,flagR2,flagL3,flagR3;
    auto Q1LL=Teno5_BVDMR(q1L,flagL1);
    auto Q1RR=Teno5_BVDMR(q1R,flagR1);
    auto Q2LL=Teno5_BVDMR(q2L,flagL2);
    auto Q2RR=Teno5_BVDMR(q2R,flagR2);
    auto Q3LL=Teno5_BVDMR(q3L,flagL3);
    auto Q3RR=Teno5_BVDMR(q3R,flagR3);


    std::array<real,2> Q1,Q2,Q3;
    Q1=minDif(Q1LL,Q1RR);
    Q2=minDif(Q2LL,Q2RR);
    Q3=minDif(Q3LL,Q3RR);
    // Q1=(flagL1&&flagR1)? minDif(Q1LL,Q1RR):minDif2(Q1LL,Q1RR);
    // Q2=(flagL2&&flagR2)? minDif(Q2LL,Q2RR):minDif2(Q2LL,Q2RR);
    // Q3=(flagL3&&flagR3)? minDif(Q3LL,Q3RR):minDif2(Q3LL,Q3RR);
    return {Q1[0]+(Q2[0]+Q3[0])/(cRef*cRef),Q1[1]+(Q2[1]+Q3[1])/(cRef*cRef)
           ,(Q2[0]-Q3[0])/cRef/rRef        ,(Q2[1]-Q3[1])/cRef/rRef
           ,Q2[0]+Q3[0]                    ,Q2[1]+Q3[1]};
}
std::vector<real> SpaceDis::recon1DBVDPrim(int i)
{
    assert(info->dim==1);
    std::vector<real> res(nPrim*2);
    auto p=res.begin();
    for (int ivar = 0; ivar < nPrim; ivar++)
    {
        std::array q1L={at(i-3,ivar),at(i-2,ivar),at(i-1,ivar),at(i,ivar),at(i+1,ivar)};
        std::array q1R={at(i+2,ivar),at(i+1,ivar),at(i,ivar),at(i-1,ivar),at(i-2,ivar)};
        bool flagL1,flagR1;
        auto Q1LL=Teno5_BVDMR(q1L,flagL1);
        auto Q1RR=Teno5_BVDMR(q1R,flagR1);
        auto Q1=minDif(Q1LL,Q1RR);
        (*p++)=Q1[0];
        (*p++)=Q1[1];
    }
    return res;
}

std::vector<real> SpaceDis::recon1DBVD2(int i)
{
    assert(info->dim==1);
    assert(info->eqType==EULER);
    enum{R,U,P};
    real rRefL=at(i-1,R);
    real pRefL=at(i-1,P);
    real cRefL=sqrt(GAMMA*(pRefL/rRefL));

    real rRefR=at(i,R);
    real pRefR=at(i,P);
    real cRefR=sqrt(GAMMA*(pRefR/rRefR));
    if(isnan(cRefL)||isnan(cRefR))
    std::cout<<"cRef error nan\n";
    std::array<real,5> q1L,q2L,q3L,q1R,q2R,q3R;
    for(int j=i-3;j<i+3;j++)
    {
        int iLocal=j-i+3;
        if(iLocal<5){
        q1L[iLocal]=at(j,R)-at(j,P)/(cRefL*cRefL);
        q2L[iLocal]=(at(j,P)+at(j,U)*cRefL*rRefL)/2;
        q3L[iLocal]=(at(j,P)-at(j,U)*cRefL*rRefL)/2;}

        iLocal=i+2-j;
        if(iLocal<5){
        q1R[iLocal]=at(j,R)-at(j,P)/(cRefR*cRefR);
        q2R[iLocal]=(at(j,P)+at(j,U)*cRefR*rRefR)/2;
        q3R[iLocal]=(at(j,P)-at(j,U)*cRefR*rRefR)/2;}
    }
    bool flagL1,flagR1,flagL2,flagR2,flagL3,flagR3;
    auto Q1LL=Teno5_BVDMR(q1L,flagL1);
    auto Q1RR=Teno5_BVDMR(q1R,flagR1);
    auto Q2LL=Teno5_BVDMR(q2L,flagL2);
    auto Q2RR=Teno5_BVDMR(q2R,flagR2);
    auto Q3LL=Teno5_BVDMR(q3L,flagL3);
    auto Q3RR=Teno5_BVDMR(q3R,flagR3);


    std::array<real,3> Q1LPrim,Q1RPrim,Q2LPrim,Q2RPrim,Q3LPrim,Q3RPrim;
    for(int j=0;j<3;j++)
    {
        Q1LPrim[j]=Q1LL[j]+(Q2LL[j]+Q3LL[j])/(cRefL*cRefL);
        Q2LPrim[j]=(Q2LL[j]-Q3LL[j])/cRefL/rRefL;
        Q3LPrim[j]=Q2LL[j]+Q3LL[j];

        Q1RPrim[j]=Q1RR[j]+(Q2RR[j]+Q3RR[j])/(cRefR*cRefR);
        Q2RPrim[j]=(Q2RR[j]-Q3RR[j])/cRefR/rRefR;
        Q3RPrim[j]=Q2RR[j]+Q3RR[j];
    }

    std::array<real,2> Q1,Q2,Q3;
    Q1=minDif(Q1LPrim,Q1RPrim);
    Q2=minDif(Q2LPrim,Q2RPrim);
    Q3=minDif(Q3LPrim,Q3RPrim);
    // Q1=!(flagL1||flagR1)? minDif(Q1LPrim,Q1RPrim):minDif2(Q1LPrim,Q1RPrim);
    // Q2=!(flagL2||flagR2)? minDif(Q2LPrim,Q2RPrim):minDif2(Q2LPrim,Q2RPrim);
    // Q3=!(flagL3||flagR3)? minDif(Q3LPrim,Q3RPrim):minDif2(Q3LPrim,Q3RPrim);
    return {Q1[0],Q1[1],Q2[0],Q2[1],Q3[0],Q3[1]};
}


std::vector<real> SpaceDis::recon2DBVD(int i)
{
    assert(info->dim==2);
    assert(info->eqType==EULER);
    enum{R,U,V,P};
    real rRef=(at(i-1,R)+at(i,R))/2;
    real pRef=(at(i-1,P)+at(i,P))/2;
    real cRef=sqrt(GAMMA*(pRef/rRef));
    std::array<real,5> q1L,q2L,q3L,q4L,q1R,q2R,q3R,q4R;
    for(int j=i-3;j<i+3;j++)
    {
        
        real Vn=norm[0]*at(j,U)+norm[1]*at(j,V);
        real q1tmp=norm[0]*at(j,V)-norm[1]*at(j,U);
        real q2tmp=at(j,R)-at(j,P)/(cRef*cRef);
        real q3tmp=(at(j,P)+Vn*cRef*rRef)/2;
        real q4tmp=(at(j,P)-Vn*cRef*rRef)/2;

        int iLocal=j-i+3;
        if(iLocal<5){
        q1L[iLocal]=q1tmp;
        q2L[iLocal]=q2tmp;
        q3L[iLocal]=q3tmp;
        q4L[iLocal]=q4tmp;}

        iLocal=i+2-j;
        if(iLocal<5){
        q1R[iLocal]=q1tmp;
        q2R[iLocal]=q2tmp;
        q3R[iLocal]=q3tmp;
        q4R[iLocal]=q4tmp;}
    }
    bool flagL1,flagR1,flagL2,flagR2,flagL3,flagR3,flagL4,flagR4;
    auto Q1LL=Teno5_BVDMR(q1L,flagL1);
    auto Q1RR=Teno5_BVDMR(q1R,flagR1);
    auto Q2LL=Teno5_BVDMR(q2L,flagL2);
    auto Q2RR=Teno5_BVDMR(q2R,flagR2);
    auto Q3LL=Teno5_BVDMR(q3L,flagL3);
    auto Q3RR=Teno5_BVDMR(q3R,flagR3);
    auto Q4LL=Teno5_BVDMR(q4L,flagL4);
    auto Q4RR=Teno5_BVDMR(q4R,flagR4);

    
    std::array<real,2> Q1,Q2,Q3,Q4;
    // Q1=(flagL1&&flagR1)? minDif(Q1LL,Q1RR):minDif2(Q1LL,Q1RR);
    // Q2=(flagL2&&flagR2)? minDif(Q2LL,Q2RR):minDif2(Q2LL,Q2RR);
    // Q3=(flagL3&&flagR3)? minDif(Q3LL,Q3RR):minDif2(Q3LL,Q3RR);
    // Q4=(flagL4&&flagR4)? minDif(Q4LL,Q4RR):minDif2(Q4LL,Q4RR);

    Q1=minDif(Q1LL,Q1RR);
    Q2=minDif(Q2LL,Q2RR);
    Q3=minDif(Q3LL,Q3RR);
    Q4=minDif(Q4LL,Q4RR);
    std::vector<real> res=
    {Q2[0]+(Q3[0]+Q4[0])/(cRef*cRef),Q2[1]+(Q3[1]+Q4[1])/(cRef*cRef)
    ,-norm[1]*Q1[0]+norm[0]*(Q3[0]-Q4[0])/cRef/rRef,-norm[1]*Q1[1]+norm[0]*(Q3[1]-Q4[1])/cRef/rRef
     ,norm[0]*Q1[0]+norm[1]*(Q3[0]-Q4[0])/cRef/rRef, norm[0]*Q1[1]+norm[1]*(Q3[1]-Q4[1])/cRef/rRef
     ,Q3[0]+Q4[0],Q3[1]+Q4[1]};
    return res;
}

std::vector<real> SpaceDis::recon2DBVD2(int i)
{
    assert(info->dim==2);
    assert(info->eqType==EULER);
    enum{R,U,V,P};
    real rRefL=at(i-1,R);
    real pRefL=at(i-1,P);
    real cRefL=sqrt(GAMMA*(pRefL/rRefL));

    real rRefR=at(i,R);
    real pRefR=at(i,P);
    real cRefR=sqrt(GAMMA*(pRefR/rRefR));
    std::array<real,5> q1L,q2L,q3L,q4L,q1R,q2R,q3R,q4R;
    for(int j=i-3;j<i+3;j++)
    {
        
        real Vn=norm[0]*at(j,U)+norm[1]*at(j,V);

        int iLocal=j-i+3;
        if(iLocal<5){
        q1L[iLocal]=norm[0]*at(j,V)-norm[1]*at(j,U);
        q2L[iLocal]=at(j,R)-at(j,P)/(cRefL*cRefL);
        q3L[iLocal]=(at(j,P)+Vn*cRefL*rRefL)/2;
        q4L[iLocal]=(at(j,P)-Vn*cRefL*rRefL)/2;}

        iLocal=i+2-j;
        if(iLocal<5){
        q1R[iLocal]=norm[0]*at(j,V)-norm[1]*at(j,U);
        q2R[iLocal]=at(j,R)-at(j,P)/(cRefR*cRefR);
        q3R[iLocal]=(at(j,P)+Vn*cRefR*rRefR)/2;
        q4R[iLocal]=(at(j,P)-Vn*cRefR*rRefR)/2;}
    }
    bool flagL1,flagR1,flagL2,flagR2,flagL3,flagR3,flagL4,flagR4;
    auto Q1LL=Teno5_BVDMR(q1L,flagL1);
    auto Q1RR=Teno5_BVDMR(q1R,flagR1);
    auto Q2LL=Teno5_BVDMR(q2L,flagL2);
    auto Q2RR=Teno5_BVDMR(q2R,flagR2);
    auto Q3LL=Teno5_BVDMR(q3L,flagL3);
    auto Q3RR=Teno5_BVDMR(q3R,flagR3);
    auto Q4LL=Teno5_BVDMR(q4L,flagL4);
    auto Q4RR=Teno5_BVDMR(q4R,flagR4);

    std::array<real,3> Q1LPrim,Q1RPrim,Q2LPrim,Q2RPrim,Q3LPrim,Q3RPrim,Q4LPrim,Q4RPrim;
    for(int j=0;j<3;j++)
    {
        Q1LPrim[j]=Q2LL[j]+(Q3LL[j]+Q4LL[j])/(cRefL*cRefL);
        Q2LPrim[j]=-norm[1]*Q1LL[j]+norm[0]*(Q3LL[j]-Q4LL[j])/cRefL/rRefL;
        Q3LPrim[j]=norm[0]*Q1LL[j]+norm[1]*(Q3LL[j]-Q4LL[j])/cRefL/rRefL;
        Q4LPrim[j]=Q3LL[j]+Q4LL[j];

        Q1RPrim[j]=Q2RR[j]+(Q3RR[j]+Q4RR[j])/(cRefR*cRefR);
        Q2RPrim[j]=-norm[1]*Q1RR[j]+norm[0]*(Q3RR[j]-Q4RR[j])/cRefR/rRefR;
        Q3RPrim[j]=norm[0]*Q1RR[j]+norm[1]*(Q3RR[j]-Q4RR[j])/cRefR/rRefR;
        Q4RPrim[j]=Q3RR[j]+Q4RR[j];
    }
    std::array<real,2> Q1,Q2,Q3,Q4;

    Q1=minDif(Q1LPrim,Q1RPrim);
    Q2=minDif(Q2LPrim,Q2RPrim);
    Q3=minDif(Q3LPrim,Q3RPrim);
    Q4=minDif(Q4LPrim,Q4RPrim);
    return {Q1[0],Q1[1],Q2[0],Q2[1],Q3[0],Q3[1],Q4[0],Q4[1]};
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
        real Vn=(norm[1]>norm[0])?at(j,V):at(j,U);
        q1[iLocal]=(norm[1]>norm[0])?at(j,U):at(j,V);
        q2[iLocal]=at(j,R)-at(j,P)/(cRef*cRef);
        q3[iLocal]=(Vn*cRef*rRef+at(j,P))/2;
        q4[iLocal]=(-Vn*cRef*rRef+at(j,P))/2;
    }
    auto start = std::chrono::high_resolution_clock::now(); 
    real Q1=(*this->inter5)(q1);
    real Q2=(*this->inter5Positive)(q2);
    real Q3=(*this->inter5)(q3);
    real Q4=(*this->inter5)(q4);
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count(); 
    timep+=duration;
    std::vector<real> res={Q2+(Q3+Q4)/(cRef*cRef)
            ,(norm[1]>norm[0]? Q1 : (Q3-Q4)/cRef/rRef )
            ,(norm[1]>norm[0]? (Q3-Q4)/cRef/rRef  : Q1)
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
        real Vn=(norm[1]>norm[0])?at(j,V):at(j,U);
        q1[iLocal]=(norm[1]>norm[0])?at(j,U):at(j,V);
        q2[iLocal]=at(j,R)-at(j,P)/(cRef*cRef);
        q3[iLocal]=(Vn*cRef*rRef+at(j,P))/2;
        q4[iLocal]=(-Vn*cRef*rRef+at(j,P))/2;
    }
    auto start = std::chrono::high_resolution_clock::now(); 
    real Q1=(*this->inter5)(q1);
    real Q2=(*this->inter5Positive)(q2);
    real Q3=(*this->inter5)(q3);
    real Q4=(*this->inter5)(q4);
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count(); 
    timep+=duration;

    std::vector<real> res={Q2+(Q3+Q4)/(cRef*cRef)
            ,(norm[1]>norm[0]? Q1 : (Q3-Q4)/cRef/rRef )
            ,(norm[1]>norm[0]? (Q3-Q4)/cRef/rRef  : Q1)
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

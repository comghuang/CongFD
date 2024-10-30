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
        // q2[iLocal]=(at(j,P)+at(j,U)*cRef*rRef)/2;
        // q3[iLocal]=(at(j,P)-at(j,U)*cRef*rRef)/2;
        q2[iLocal]=at(j,P);
        q3[iLocal]=at(j,U);
    }
    auto start = std::chrono::high_resolution_clock::now(); 
    real Q1=(*this->inter5Positive)(q1);
    real Q2=(*this->inter5)(q2);
    real Q3=(*this->inter5)(q3);
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count(); 
    timep+=duration;
    // return {Q1+(Q2+Q3)/(cRef*cRef),(Q2-Q3)/cRef/rRef,Q2+Q3};
    return {Q1+(Q2)/(cRef*cRef),Q3,Q2};
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
        // q2[iLocal]=(at(j,P)+at(j,U)*cRef*rRef)/2;
        // q3[iLocal]=(at(j,P)-at(j,U)*cRef*rRef)/2;
        q2[iLocal]=at(j,P);
        q3[iLocal]=at(j,U);
    }
    auto start = std::chrono::high_resolution_clock::now(); 
    real Q1=(*this->inter5Positive)(q1);
    real Q2=(*this->inter5)(q2);
    real Q3=(*this->inter5)(q3);
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count(); 
    timep+=duration;
    // return {Q1+(Q2+Q3)/(cRef*cRef),(Q2-Q3)/cRef/rRef,Q2+Q3};
    return {Q1+(Q2)/(cRef*cRef),Q3,Q2};
}

static std::array<real,2> minDif(real u1L,real u2L,real u1R,real u2R)
{
    std::array<real,3> difs={abs(u1L-u1R),abs(u1L-u2R),abs(u2L-u1R)};
    int index=std::min_element(difs.begin(),difs.end())-difs.begin();
    if(index==0) return{u1L,u1R};
    if(index==1) return{u1L,u2R};
    if(index==2) return{u2L,u1R};
    std::cout<<"spRecon Error: minDif index out\n";
    return{0,0};
}

static std::array<real,2> minDif(std::array<real,3> uL,std::array<real,3> uR)
{
    std::array<real,2> difs={abs(uL[0]-uR[0]),abs(uL[1]-uR[1])};
    int index=std::min_element(difs.begin(),difs.end())-difs.begin();
         if(index==0) return{uL[0],uR[0]};
    else if(index==1) return{uL[1],uR[1]};
    std::cout<<"spRecon Error: minDif index out\n";
    return{0,0};
}

static std::array<real,2> minDif(std::array<real,2> uL,std::array<real,2> uR)
{
    // std::array<real,4> difs={abs(uL[0]-uR[0]),abs(uL[1]-uR[1]),abs(uL[0]-uR[1]),abs(uL[1]-uR[0])};
    std::array<real,2> difs={abs(uL[0]-uR[0]),abs(uL[1]-uR[1])};
    int index=std::min_element(difs.begin(),difs.end())-difs.begin();
         if(index==0) return{uL[0],uR[0]};
    else if(index==1) 
    return{uL[1],uR[1]};
    else if(index==2) return{uL[0],uR[1]};
    else if(index==3) return{uL[1],uR[0]};
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
std::vector<real> SpaceDis::recon1DFaceCenter(int i)
{
    std::array<real,3> primL,primR;
    memcpy(&primL[0],&at(i-1,0),nVar*sizeof(real));
    memcpy(&primR[0],&at(i,0),nVar*sizeof(real));
    eigensystemEuler1D eig=eigensystemEuler1D(primL,primR);
    std::array<real,5> q1L,q2L,q3L,q1R,q2R,q3R;
    for(int j=i-3;j<i+3;j++)
    {
        enum{R,U,P};
        auto charTemp=eig.primToChar({at(j,R),at(j,U),at(j,P)});

        int iLocal=j-i+3;
        if(iLocal<5){
        q1L[iLocal]=charTemp[0];
        q2L[iLocal]=charTemp[1];
        q3L[iLocal]=charTemp[2];}

        iLocal=i+2-j;
        if(iLocal<5){
        q1R[iLocal]=charTemp[0];
        q2R[iLocal]=charTemp[1];
        q3R[iLocal]=charTemp[2];}
    }
    // auto start = std::chrono::steady_clock::now(); 
    auto Q1LL=inter5(q1L);
    auto Q1RR=inter5(q1R);
    auto Q2LL=inter5(q2L);
    auto Q2RR=inter5(q2R);
    auto Q3LL=inter5(q3L);
    auto Q3RR=inter5(q3R);

    auto Q1LLTHINC=THINC1(q1L[1],q1L[2],q1L[3]);
    auto Q1RRTHINC=THINC1(q1R[1],q1R[2],q1R[3]);
    auto Q2LLTHINC=THINC1(q2L[1],q2L[2],q2L[3]);
    auto Q2RRTHINC=THINC1(q2R[1],q2R[2],q2R[3]);
    auto Q3LLTHINC=THINC1(q3L[1],q3L[2],q3L[3]);
    auto Q3RRTHINC=THINC1(q3R[1],q3R[2],q3R[3]);

    std::array<real,2> Q1,Q2,Q3;
    Q1=minDif((std::array<real,2>){Q1LL,Q1LLTHINC},(std::array<real,2>){Q1RR,Q1RRTHINC});
    Q2=minDif((std::array<real,2>){Q2LL,Q2LLTHINC},(std::array<real,2>){Q2RR,Q2RRTHINC});
    Q3=minDif((std::array<real,2>){Q3LL,Q3LLTHINC},(std::array<real,2>){Q3RR,Q3RRTHINC});
    Q1LL=Q1[0];Q1RR=Q1[1];
    Q2LL=Q2[0];Q2RR=Q2[1];
    Q3LL=Q3[0];Q3RR=Q3[1];
    // auto stop = std::chrono::steady_clock::now(); 
    // auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count(); 
    // timep+=duration;

    auto resTempL=eig.charToPrim({Q1LL,Q2LL,Q3LL});
    auto resTempR=eig.charToPrim({Q1RR,Q2RR,Q3RR});
    return {resTempL[0],resTempL[1],resTempL[2],resTempR[0],resTempR[1],resTempR[2]};
}


std::vector<real> SpaceDis::recon2DFaceCenter(int i)
{
    std::array<real,4> primL,primR;
    memcpy(&primL[0],&at(i-1,0),nVar*sizeof(real));
    memcpy(&primR[0],&at(i,0),nVar*sizeof(real));
    eigensystemEuler2D eig=eigensystemEuler2D(primL,primR,norm);
    std::array<real,5> q1L,q2L,q3L,q4L,q1R,q2R,q3R,q4R;
    for(int j=i-3;j<i+3;j++)
    {
        enum{R,U,V,P};
        auto charTemp=eig.primToChar({at(j,R),at(j,U),at(j,V),at(j,P)});

        int iLocal=j-i+3;
        if(iLocal<5){
        q1L[iLocal]=charTemp[0];
        q2L[iLocal]=charTemp[1];
        q3L[iLocal]=charTemp[2];
        q4L[iLocal]=charTemp[3];}

        iLocal=i+2-j;
        if(iLocal<5){
        q1R[iLocal]=charTemp[0];
        q2R[iLocal]=charTemp[1];
        q3R[iLocal]=charTemp[2];
        q4R[iLocal]=charTemp[3];}
    }
    bool flagL1,flagR1,flagL2,flagR2,flagL3,flagR3,flagL4,flagR4;
    auto start = std::chrono::steady_clock::now(); 
    auto Q1LL=inter5(q1L);
    auto Q1RR=inter5(q1R);
    auto Q2LL=inter5(q2L);
    auto Q2RR=inter5(q2R);
    auto Q3LL=inter5(q3L);
    auto Q3RR=inter5(q3R);
    auto Q4LL=inter5(q4L);
    auto Q4RR=inter5(q4R);
    auto stop = std::chrono::steady_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count(); 
    timep+=duration;

    auto resTempL=eig.charToPrim({Q1LL,Q2LL,Q3LL,Q4LL});
    auto resTempR=eig.charToPrim({Q1RR,Q2RR,Q3RR,Q4RR});
    return {resTempL[0],resTempL[1],resTempL[2],resTempL[3],resTempR[0],resTempR[1],resTempR[2],resTempR[3]};
}
std::vector<real> SpaceDis::reconLChar2D(int i)
{
    assert(info->dim==2);
    assert(info->eqType==EULER);
    enum{R,U,V,P};
    eigensystemEuler2D eig=eigensystemEuler2D({at(i-1,R),at(i-1,U),at(i-1,V),at(i-1,P)},norm);
    std::array<real,5> q1,q2,q3,q4;
    for(int j=i-3;j<i+2;j++)
    {
        int iLocal=j-i+3;
        auto charTemp=eig.primToChar({at(j,R),at(j,U),at(j,V),at(j,P)});
        q1[iLocal]=charTemp[0];
        q2[iLocal]=charTemp[1];
        q3[iLocal]=charTemp[2];
        q4[iLocal]=charTemp[3];
    }
    auto start = std::chrono::high_resolution_clock::now(); 
    real Q1=(*this->inter5)(q1);
    real Q2=(*this->inter5Positive)(q2);
    real Q3=(*this->inter5)(q3);
    real Q4=(*this->inter5)(q4);
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count(); 
    timep+=duration;
    auto resTemp=eig.charToPrim({Q1,Q2,Q3,Q4});
    return {resTemp[0],resTemp[1],resTemp[2],resTemp[3]};
}

std::vector<real> SpaceDis::reconRChar2D(int i)
{
    assert(info->dim==2);
    
    enum{R,U,V,P};
    eigensystemEuler2D eig=eigensystemEuler2D({at(i,R),at(i,U),at(i,V),at(i,P)},norm);
    std::array<real,5> q1,q2,q3,q4;
    for(int j=i+2;j>i-3;j--)
    {
        int iLocal=i+2-j;
        auto charTemp=eig.primToChar({at(j,R),at(j,U),at(j,V),at(j,P)});
        q1[iLocal]=charTemp[0];
        q2[iLocal]=charTemp[1];
        q3[iLocal]=charTemp[2];
        q4[iLocal]=charTemp[3];
    }
    auto start = std::chrono::high_resolution_clock::now(); 
    real Q1=(*this->inter5)(q1);
    real Q2=(*this->inter5Positive)(q2);
    real Q3=(*this->inter5)(q3);
    real Q4=(*this->inter5)(q4);
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count(); 
    timep+=duration;
    auto resTemp=eig.charToPrim({Q1,Q2,Q3,Q4});
    return {resTemp[0],resTemp[1],resTemp[2],resTemp[3]};
}

// std::vector<real> SpaceDis::reconLChar2D(int i)
// {
//     assert(info->dim==2);
//     assert(info->eqType==EULER);
//     enum{R,U,V,P};
//     double rhoAvg,uAvg,vAvg,HAvg,cAvg,VnAvg,q_2Avg,coef1,coef2;
//     std::array<real,2> H;
//     real rl=at(i-1,R),ul=at(i-1,U),vl=at(i-1,V),pl=at(i-1,P);
//     real rr=at(i,R),ur=at(i,U),vr=at(i,V),pr=at(i,P);
//     enum{LL,RR};
//     H[LL]=(ul*ul+vl*vl)/2+pl/rl*GAMMA/(GAMMA-1);
//     H[RR]=(ur*ur+vr*vr)/2+pr/rr*GAMMA/(GAMMA-1);
//     coef1=sqrt(rl);
//     coef2=sqrt(rr);
//     real divisor=1.0/(sqrt(rl)+sqrt(rr));
//     real rRef=sqrt(rl*rr);
//     uAvg=(coef1*ul+coef2*ur)*divisor;
//     vAvg=(coef1*vl+coef2*vr)*divisor;
//     HAvg=(coef1*H[LL]+coef2*H[RR])*divisor;
//     q_2Avg=(uAvg*uAvg+vAvg*vAvg)*0.5;
//     real cRef=sqrt((GAMMA-1)*(HAvg-q_2Avg));
//     std::array<real,5> q1,q2,q3,q4;
//     for(int j=i-3;j<i+2;j++)
//     {
//         int iLocal=j-i+3;
//         real Vn=(norm[1]>norm[0])?at(j,V):at(j,U);
//         q1[iLocal]=(norm[1]>norm[0])?at(j,U):at(j,V);
//         q2[iLocal]=at(j,R)-at(j,P)/(cRef*cRef);
//         q3[iLocal]=Vn+at(j,P)/(cRef*rRef);
//         q4[iLocal]=-Vn+at(j,P)/(cRef*rRef);
//     }
//     auto start = std::chrono::high_resolution_clock::now(); 
//     real Q1=(*this->inter5)(q1);
//     real Q2=(*this->inter5Positive)(q2);
//     real Q3=(*this->inter5)(q3);
//     real Q4=(*this->inter5)(q4);
//     auto stop = std::chrono::high_resolution_clock::now(); 
//     auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count(); 
//     timep+=duration;
//     std::vector<real> res={Q2+(Q3+Q4)*0.5*rRef/cRef
//             ,(norm[1]>norm[0]? Q1 : (Q3-Q4)*0.5 )
//             ,(norm[1]>norm[0]? (Q3-Q4)*0.5  : Q1)
//             ,(Q3+Q4)*(cRef*rRef)*0.5};
//     return res;
// }

// std::vector<real> SpaceDis::reconRChar2D(int i)
// {
//     assert(info->dim==2);
//     enum{R,U,V,P};

//     double rhoAvg,uAvg,vAvg,HAvg,cAvg,VnAvg,q_2Avg,coef1,coef2;
//     std::array<real,2> H;
//     real rl=at(i-1,R),ul=at(i-1,U),vl=at(i-1,V),pl=at(i-1,P);
//     real rr=at(i,R),ur=at(i,U),vr=at(i,V),pr=at(i,P);
//     enum{LL,RR};
//     H[LL]=(ul*ul+vl*vl)/2+pl/rl*GAMMA/(GAMMA-1);
//     H[RR]=(ur*ur+vr*vr)/2+pr/rr*GAMMA/(GAMMA-1);
//     coef1=sqrt(rl);
//     coef2=sqrt(rr);
//     real divisor=1.0/(sqrt(rl)+sqrt(rr));
//     real rRef=sqrt(rl*rr);
//     uAvg=(coef1*ul+coef2*ur)*divisor;
//     vAvg=(coef1*vl+coef2*vr)*divisor;
//     HAvg=(coef1*H[LL]+coef2*H[RR])*divisor;
//     q_2Avg=(uAvg*uAvg+vAvg*vAvg)*0.5;
//     real cRef=sqrt((GAMMA-1)*(HAvg-q_2Avg));
    

//     std::array<real,5> q1,q2,q3,q4;
//     for(int j=i+2;j>i-3;j--)
//     {
//         int iLocal=i+2-j;
//         real Vn=(norm[1]>norm[0])?at(j,V):at(j,U);
//         q1[iLocal]=(norm[1]>norm[0])?at(j,U):at(j,V);
//         q2[iLocal]=at(j,R)-at(j,P)/(cRef*cRef);
//         q3[iLocal]=Vn+at(j,P)/(cRef*rRef);
//         q4[iLocal]=-Vn+at(j,P)/(cRef*rRef);
//     }
//     auto start = std::chrono::high_resolution_clock::now(); 
//     real Q1=(*this->inter5)(q1);
//     real Q2=(*this->inter5Positive)(q2);
//     real Q3=(*this->inter5)(q3);
//     real Q4=(*this->inter5)(q4);
//     auto stop = std::chrono::high_resolution_clock::now(); 
//     auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count(); 
//     timep+=duration;
//     std::vector<real> res={Q2+(Q3+Q4)*0.5*rRef/cRef
//             ,(norm[1]>norm[0]? Q1 : (Q3-Q4)*0.5 )
//             ,(norm[1]>norm[0]? (Q3-Q4)*0.5  : Q1)
//             ,(Q3+Q4)*(cRef*rRef)*0.5};
//     return res;
// }

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

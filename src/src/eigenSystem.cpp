#include"eigenSystem.hpp"

eigensystemEuler2D::eigensystemEuler2D(const std::array<real,4> &prim,const std::array<real,3> & norm_)
{
    r=prim[0],u=prim[1],v=prim[2],p=prim[3];
    gamma=GAMMA,ek=(u*u+v*v)/2;
    h=p/r*gamma/(1-gamma);
    c=sqrt(gamma*p/r);
    norm=norm_;
    Vn=norm[0]*u+norm[1]*v;
}

eigensystemEuler2D::eigensystemEuler2D(const std::array<real,4> &priml,const std::array<real,4> &primr,const std::array<real,3> & norm_)
{
      norm=norm_;
      
      gamma=GAMMA;
      enum{LL,RR};
      real rl=priml[0],ul=priml[1],vl=priml[2],pl=priml[3];
      real rr=primr[0],ur=primr[1],vr=primr[2],pr=primr[3];

      std::array<real,2> H;
      H[LL]=(ul*ul+vl*vl)/2+pl/rl*gamma/(gamma-1);
      H[RR]=(ur*ur+vr*vr)/2+pr/rr*gamma/(gamma-1);
      real coef1=sqrt(rl);
      real coef2=sqrt(rr);
      real divisor=1.0/(sqrt(rl)+sqrt(rr));

      r=sqrt(rl*rr);
      u=(coef1*ul+coef2*ur)*divisor;
      v=(coef1*vl+coef2*vr)*divisor;
      Vn=norm[0]*u+norm[1]*v;
      real ht=(coef1*H[LL]+coef2*H[RR])*divisor;
      ek=(u*u+v*v)*0.5;
      h=ht-ek;
      c=sqrt((gamma-1)*h);
      p=r*h*((gamma-1)/gamma);

      leftEig={norm[Y]              ,(-norm[X])             ,0          ,(-norm[Y]*u+norm[X]*v),
               u/h                  , v/h                   ,-1.0/h     , (h-ek)/h,
               (-norm[X]/c-u/h)/2   ,(-norm[Y]/c-v/h)/2     ,1.0/(2.0*h),(Vn/c+ek/h)/2,
               (norm[X]/c-u/h)/2    ,(norm[Y]/c-v/h)/2      ,1.0/(2.0*h),(-Vn/c+ek/h)/2};
      
      rightEig={1 ,1          ,1          ,0,
                u ,u-norm[X]*c,u+norm[X]*c,norm[Y] ,
                v ,v-norm[Y]*c,v+norm[Y]*c,-norm[X],
                ek,h+ek-Vn*c  ,h+ek+Vn*c  ,norm[Y]*u-norm[X]*v};
}


std::array<real,4> eigensystemEuler2D::primToChar(const std::array<real,4> &prim)
{
    real rt=prim[0],ut=prim[1],vt=prim[2],pt=prim[3];
    real ekt=(ut*ut+vt*vt)/2;
    real rut=rt*ut,rvt=rt*vt,ret=pt/(gamma-1)+rt*ekt;

    std::array<real,4> res;

    res[0]=rut* leftEig[0]
          +rvt* leftEig[1]
          +ret* leftEig[2]
          +rt * leftEig[3];

    res[1]=rut* leftEig[4]
          +rvt* leftEig[5]
          +ret* leftEig[6]
          +rt * leftEig[7];

    res[2]=rut* leftEig[8]
          +rvt* leftEig[9]
          +ret* leftEig[10]
          +rt * leftEig[11];

    res[3]=rut* leftEig[12]
          +rvt* leftEig[13]
          +ret* leftEig[14]
          +rt * leftEig[15];
    
    return res;

}
std::array<real,4> eigensystemEuler2D::charToPrim(const std::array<real,4> & chars)
{
    real ch1=chars[0],ch2=chars[1],ch3=chars[2],ch4=chars[3],rt,rut,rvt,ret;
    rt    =ch2* rightEig[0]
          +ch3* rightEig[1]
          +ch4* rightEig[2]
          +ch1* rightEig[3];

    rut   =ch2* rightEig[4]
          +ch3* rightEig[5]
          +ch4* rightEig[6]
          +ch1* rightEig[7];

    rvt   =ch2* rightEig[8]
          +ch3* rightEig[9]
          +ch4* rightEig[10]
          +ch1* rightEig[11];

    ret   =ch2* rightEig[12]
          +ch3* rightEig[13]
          +ch4* rightEig[14]
          +ch1* rightEig[15];

    real ut=rut/rt;
    real vt=rvt/rt;
    real Et=ret/rt;
    real ekt=(ut*ut+vt*vt)/2;
    real pt=(gamma-1)*(ret-rt*ekt);
    return{rt,ut,vt,pt};
}

eigensystemEuler1D::eigensystemEuler1D(const std::array<real,3> &priml,const std::array<real,3> &primr)
{
      gamma=GAMMA;
      enum{LL,RR};
      real rl=priml[0],ul=priml[1],pl=priml[2];
      real rr=primr[0],ur=primr[1],pr=primr[2];

      std::array<real,2> H;
      H[LL]=(ul*ul)/2+pl/rl*gamma/(gamma-1);
      H[RR]=(ur*ur)/2+pr/rr*gamma/(gamma-1);
      real coef1=sqrt(rl);
      real coef2=sqrt(rr);
      real divisor=1.0/(sqrt(rl)+sqrt(rr));

      r=sqrt(rl*rr);
      u=(coef1*ul+coef2*ur)*divisor;
      real ht=(coef1*H[LL]+coef2*H[RR])*divisor;
      ek=(u*u)/2;
      h=ht-ek;
      c=sqrt((gamma-1)*h);
      p=r*h*((gamma-1)/gamma);

      real gamma_1=gamma-1;
      leftEig={ht+c*(u-c)/gamma_1      , -u-c/gamma_1   , 1.0,
               -2.0*ht+4.0*c*c/gamma_1 , 2.0*u          ,-2.0,
               ht-c*(u+c)/gamma_1      , -u+c/gamma_1   , 1.0};
      real factorEig=0.5*gamma_1/(c*c);


      for(int ii=0;ii<9;ii++) leftEig[ii]*=factorEig;

      rightEig={1.0  , 1.0     , 1.0,
               u-c   , u       , u+c,
               ht-u*c , 0.5*u*u , ht+u*c};
}


std::array<real,3> eigensystemEuler1D::primToChar(const std::array<real,3> &prim)
{
      real rt=prim[0],ut=prim[1],pt=prim[2];
      real ekt=(ut*ut)/2;
      real rut=rt*ut,ret=pt/(gamma-1)+rt*ekt;

      std::array<real,3> res;

      res[0]=rt * leftEig[0]
            +rut* leftEig[1]
            +ret* leftEig[2];

      res[1]=rt * leftEig[3]
            +rut* leftEig[4]
            +ret* leftEig[5];

      res[2]=rt * leftEig[6]
            +rut* leftEig[7]
            +ret* leftEig[8];

    return res;
}

std::array<real,3> eigensystemEuler1D::charToPrim(const std::array<real,3> & chars)
{
    real ch1=chars[0],ch2=chars[1],ch3=chars[2],rt,rut,ret;

    rt    =ch1* rightEig[0]
          +ch2* rightEig[1]
          +ch3* rightEig[2];

    rut   =ch1* rightEig[3]
          +ch2* rightEig[4]
          +ch3* rightEig[5];

    ret   =ch1* rightEig[6]
          +ch2* rightEig[7]
          +ch3* rightEig[8];

    real ut=rut/rt;
    real Et=ret/rt;
    real ekt=(ut*ut)/2;
    real pt=(gamma-1)*(ret-rt*ekt);
    return{rt,ut,pt};
}
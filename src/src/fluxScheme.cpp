
#include "macro.hpp"
#include "fluxScheme.hpp"
#include <algorithm>


enum{
    L,
    R
};

typedef std::array<real,2> arr2;
std::vector<real> roeFlux1D(arr2 r, arr2 u, arr2 p, arr2 H, arr2 RT)
{
    std::vector<real> res;
    res.resize(3);
    
    real gamma=GAMMA;

    real rBar=sqrt(r[L]*r[R]);
    real uBar=(u[L]*sqrt(r[L])+u[R]*sqrt(r[R]))/(sqrt(r[L])+sqrt(r[R]));
    real HBar=(H[L]*sqrt(r[L])+H[R]*sqrt(r[R]))/(sqrt(r[L])+sqrt(r[R]));
    real cBar=sqrt((gamma-1)*(HBar-uBar*uBar/2));
    real cBar2=(gamma-1)*(HBar-uBar*uBar/2);

    real dr=r[R]-r[L];
    real du=u[R]-u[L];
    real dp=p[R]-p[L];

    real K1[3]={1,uBar,uBar*uBar/2};
    real K2[3]={1,uBar-cBar,HBar-uBar*cBar};
    real K3[3]={1,uBar+cBar,HBar+uBar*cBar};

    real alpha[3]={
        dr-dp/cBar2,
        1.0/(2.0*cBar2)*(dp-rBar*cBar*du),
        1.0/(2.0*cBar2)*(dp+rBar*cBar*du)
    };
    real lambda[3]={
        abs(uBar),abs(uBar-cBar),abs(uBar+cBar)
    };
    /*
    real eps=0.1*(abs(uBar)+abs(uBar)+cBar);
    for(int i=0;i<3;i++)
    {
        if(lambda[i]<eps)
        {
            lambda[i]=(lambda[i]*lambda[i]+eps*eps)/(2.0*eps);
        }
    }*/

    real FL[3]={
        r[L]*u[L],
        r[L]*u[L]*u[L]+p[L],
        r[L]*H[L]*u[L]
    };
    real FR[3]={
        r[R]*u[R],
        r[R]*u[R]*u[R]+p[R],
        r[R]*H[R]*u[R]
    };

    for (int i = 0; i < 3; i++)
    {
        res[i]=(FL[i]+FR[i]
              -lambda[0]*alpha[0]*K1[i]
              -lambda[1]*alpha[1]*K2[i]
              -lambda[2]*alpha[2]*K3[i]
                )/2.0;
    }
    return res;


}

std::vector<real> roeFlux1D2(arr2 r, arr2 u, arr2 p, arr2 H, arr2 RT)
{
    //reference: https://blog.csdn.net/Tankrun1997/article/details/132743487
    std::vector<real> res;
    res.resize(3);
    
    real gamma=GAMMA;

    real rBar=sqrt(r[L]*r[R]);
    real uBar=(u[L]*sqrt(r[L])+u[R]*sqrt(r[R]))/(sqrt(r[L])+sqrt(r[R]));
    real HBar=(H[L]*sqrt(r[L])+H[R]*sqrt(r[R]))/(sqrt(r[L])+sqrt(r[R]));
    real cBar=sqrt((gamma-1)*(HBar-uBar*uBar/2));
    real cBar2=(gamma-1)*(HBar-uBar*uBar/2);

    real dr=r[R]-r[L];
    real du=u[R]-u[L];
    real dp=p[R]-p[L];

    real K1[3]={1,uBar,uBar*uBar/2};
    real K2[3]={1,uBar-cBar,HBar-uBar*cBar};
    real K3[3]={1,uBar+cBar,HBar+uBar*cBar};


    real lambda[3]={
        abs(uBar),abs(uBar+cBar),abs(uBar-cBar)
    };
    
    real eps=0.05*(abs(uBar)+cBar);
    for(int i=0;i<3;i++)
    {
        if(lambda[i]<eps)
        {
            lambda[i]=(lambda[i]*lambda[i]+eps*eps)/(2.0*eps);
        }
    }

    real alpha[5]={
        lambda[0]*(dr-dp/cBar2),
        lambda[1]/(2.0*cBar2)*(dp+rBar*cBar*du),
        lambda[2]/(2.0*cBar2)*(dp-rBar*cBar*du)
    };
    alpha[3]=alpha[0]+alpha[1]+alpha[2];
    alpha[4]=cBar*(alpha[1]-alpha[2]);

    real FL[3]={
        r[L]*u[L],
        r[L]*u[L]*u[L]+p[L],
        r[L]*H[L]*u[L]
    };
    real FR[3]={
        r[R]*u[R],
        r[R]*u[R]*u[R]+p[R],
        r[R]*H[R]*u[R]
    };

    for (int i = 0; i < 3; i++)
    {
        res[i]=(FL[i]+FR[i])/2.0;
    }
    res[0]-=0.5*alpha[3];
    res[1]-=0.5*(alpha[3]*uBar+alpha[4]);
    res[2]-=0.5*(HBar*alpha[3]+uBar*alpha[4]-cBar2*alpha[0]/(gamma-1));



    return res;


}





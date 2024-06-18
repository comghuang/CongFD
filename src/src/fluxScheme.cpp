
#include "macro.hpp"
#include "fluxScheme.hpp"
#include <algorithm>


typedef std::array<real,2> arr2;
std::vector<real> roeFlux1D(arr2 r, arr2 u, arr2 p, arr2 H, arr2 RT)
{
    std::vector<real> res;
    res.resize(3);
    
    real gamma=GAMMA;

    real rBar=sqrt(r[LEFTT]*r[RIGHT]);
    real uBar=(u[LEFTT]*sqrt(r[LEFTT])+u[RIGHT]*sqrt(r[RIGHT]))/(sqrt(r[LEFTT])+sqrt(r[RIGHT]));
    real HBar=(H[LEFTT]*sqrt(r[LEFTT])+H[RIGHT]*sqrt(r[RIGHT]))/(sqrt(r[LEFTT])+sqrt(r[RIGHT]));
    real cBar=sqrt((gamma-1)*(HBar-uBar*uBar/2));
    real cBar2=(gamma-1)*(HBar-uBar*uBar/2);

    real dr=r[RIGHT]-r[LEFTT];
    real du=u[RIGHT]-u[LEFTT];
    real dp=p[RIGHT]-p[LEFTT];

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
        r[LEFTT]*u[LEFTT],
        r[LEFTT]*u[LEFTT]*u[LEFTT]+p[LEFTT],
        H[LEFTT]*u[LEFTT]
    };
    real FR[3]={
        r[RIGHT]*u[RIGHT],
        r[RIGHT]*u[RIGHT]*u[RIGHT]+p[RIGHT],
        H[RIGHT]*u[RIGHT]
    };

    for (ind i = 0; i < 3; i++)
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

    real rBar=sqrt(r[LEFTT]*r[RIGHT]);
    real uBar=(u[LEFTT]*sqrt(r[LEFTT])+u[RIGHT]*sqrt(r[RIGHT]))/(sqrt(r[LEFTT])+sqrt(r[RIGHT]));
    real HBar=(H[LEFTT]*sqrt(r[LEFTT])+H[RIGHT]*sqrt(r[RIGHT]))/(sqrt(r[LEFTT])+sqrt(r[RIGHT]));
    real cBar=sqrt((gamma-1)*(HBar-uBar*uBar/2));
    real cBar2=(gamma-1)*(HBar-uBar*uBar/2);

    real dr=r[RIGHT]-r[LEFTT];
    real du=u[RIGHT]-u[LEFTT];
    real dp=p[RIGHT]-p[LEFTT];

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
        r[LEFTT]*u[LEFTT],
        r[LEFTT]*u[LEFTT]*u[LEFTT]+p[LEFTT],
        r[LEFTT]*H[LEFTT]*u[LEFTT]
    };
    real FR[3]={
        r[RIGHT]*u[RIGHT],
        r[RIGHT]*u[RIGHT]*u[RIGHT]+p[RIGHT],
        r[RIGHT]*H[RIGHT]*u[RIGHT]
    };

    for (ind i = 0; i < 3; i++)
    {
        res[i]=(FL[i]+FR[i])/2.0;
    }
    res[0]-=0.5*alpha[3];
    res[1]-=0.5*(alpha[3]*uBar+alpha[4]);
    res[2]-=0.5*(HBar*alpha[3]+uBar*alpha[4]-cBar2*alpha[0]/(gamma-1));



    return res;


}

std::vector<real> HLLCFlux1D(arr2 r, arr2 u, arr2 p, arr2 H, arr2 RT)
{
    //reference:https://zhuanlan.zhihu.com/p/583555029
    std::vector<real> res;
    res.resize(3);
    real gamma=GAMMA;
    real cl=sqrt(gamma*RT[LEFTT]);
    real cr=sqrt(gamma*RT[RIGHT]);

    real SL=std::min(u[LEFTT]-cl,u[RIGHT]-cr);
    real SR=std::max(u[LEFTT]+cl,u[RIGHT]+cr);

    real Sstar=(r[LEFTT]*u[LEFTT]*(SL-u[LEFTT])-p[LEFTT]
                +p[RIGHT]-r[RIGHT]*u[RIGHT]*(SR-u[RIGHT]))/
               (r[LEFTT]*(SL-u[LEFTT])-r[RIGHT]*(SR-u[RIGHT]));
    if (SL>=0)
    {
        res[0]=r[LEFTT]*u[LEFTT];
        res[1]=r[LEFTT]*u[LEFTT]*u[LEFTT]+p[LEFTT];
        res[2]=r[LEFTT]*H[LEFTT]*u[LEFTT];
    }
    else if (Sstar>=0)
    {
        real pStar=p[LEFTT]+r[LEFTT]*(SL-u[LEFTT])*(Sstar-u[LEFTT]);
        real U[3]={r[LEFTT],r[LEFTT]*u[LEFTT],p[LEFTT]/(gamma-1)+r[LEFTT]*u[LEFTT]*u[LEFTT]/2};
        real F[3]={r[LEFTT]*u[LEFTT],r[LEFTT]*u[LEFTT]*u[LEFTT]+p[LEFTT],r[LEFTT]*H[LEFTT]*u[LEFTT]};
        real D[3]={0,1,Sstar};
        for (ind i = 0; i < 3; i++)
        {
            res[i]=(Sstar*(SL*U[i]-F[i])+SL*pStar*D[i])/(SL-Sstar);
        }
        
    }
    else if (SR>=0)
    {
        real pStar=p[RIGHT]+r[RIGHT]*(SR-u[RIGHT])*(Sstar-u[RIGHT]);
        real U[3]={r[RIGHT],r[RIGHT]*u[RIGHT],p[RIGHT]/(gamma-1)+r[RIGHT]*u[RIGHT]*u[RIGHT]/2};
        real F[3]={r[RIGHT]*u[RIGHT],r[RIGHT]*u[RIGHT]*u[RIGHT]+p[RIGHT],r[RIGHT]*H[RIGHT]*u[RIGHT]};
        real D[3]={0,1,Sstar};
        for (ind i = 0; i < 3; i++)
        {
            res[i]=(Sstar*(SR*U[i]-F[i])+SR*pStar*D[i])/(SR-Sstar);
        }
    }
    else
    {
        res[0]=r[RIGHT]*u[RIGHT];
        res[1]=r[RIGHT]*u[RIGHT]*u[RIGHT]+p[RIGHT];
        res[2]=r[RIGHT]*H[RIGHT]*u[RIGHT];
    }
    
    
    


    return res;
}

#include <iostream>
#include "oneDDiscrete.hpp"
#include "macro.hpp"
#include "helloConfig.h"
#include "burgers.hpp"
int main(){
    
    std::vector<real> data;
    std::vector<real> x;
    std::fstream file("data.txt",std::ios::out);
    ind n=200;
    ind nVar=1;
    ind T=1;
    data.resize(n*nVar,0.0);
    x.resize(n,0.0);
    for(ind i=0;i<n;i++)
    {
        real h=2.0/n;
        real xi=h/2+i*h-1;
        file<<xi<<' ';
        x[i]=xi;
        data[i]=0.5*sin(M_PI*xi)+0.25;
    }
    file<<'\n';
    oneDDiscrete discrete(&data,n,nVar);
    
    
    real dt=0.001,t=0.0,tend=4-dt/4;
    while(t<tend)
    {
        std::vector<real> tempdata=data;

        //third order RK
        //stage 1
        std::vector<real>rhs=discrete.difference();
        for(ind i=0;i<n;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            data[i*nVar+ivar]=tempdata[i*nVar+ivar]-dt*rhs[i*nVar+ivar];
        }

        //stage 2
        rhs=discrete.difference();
        for(ind i=0;i<n;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            data[i*nVar+ivar]=0.75*tempdata[i*nVar+ivar]
                             -0.25*dt*rhs[i*nVar+ivar]
                             +0.25*data[i*nVar+ivar];
        }

        //stage 3
        rhs=discrete.difference();
        for(ind i=0;i<n;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            data[i*nVar+ivar]=1.0/3.0*tempdata[i*nVar+ivar]
                             -2.0/3.0*dt*rhs[i*nVar+ivar]
                             +2.0/3.0*data[i*nVar+ivar];
        }

        t=t+dt;
        real error=0.0;
        for(ind i=0;i<n;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            real exactSol=0;
            real err=abs(data[i*nVar+ivar]-exactSol);
            error+=err;

            file<<data[i*nVar+ivar]<<' ';
        }
        file<<'\n';
        std::cout<<"t=  "<<t<<"  error:  "<<error/n<<'\n';
        
    }


    //for (auto it = data.begin(); it != data.end() ; ++it) std::cout << *it << " ";


    return 0;
    
}
#include <iostream>
#include "oneDDiscrete.hpp"
#include "macro.hpp"
#include "helloConfig.h"
int main(){
    
    std::vector<real> data;
    std::vector<real> x;
    std::fstream file("test.txt",std::ios::out);
    ind n=40;
    ind nVar=1;
    ind T=1;
    data.resize(n*nVar,0.0);
    x.resize(n,0.0);
    for(ind i=0;i<n;i++)
    {
        real h=1.0/(n+1);
        real xi=h/2+i*h;
        file<<xi<<' ';
        x[i]=xi;
        data[i]=sin(2*M_PI/T*xi)+2;
    }
    file<<'\n';
    oneDDiscrete discrete(&data,n,nVar);
    
    
    real dt=0.01,t=0.0,tend=1.0+dt/2;
    while(t<tend)
    {
        std::vector<real>rhs=discrete.difference();
        real error=0.0;
        for(ind i=0;i<n;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            data[i*nVar+ivar]-=dt*rhs[i*nVar+ivar];
            t=t+dt;

            real exactSol=sin(2*M_PI/T*(x[i]-t))+2;
            error+=abs(data[i*nVar+ivar]-exactSol);

            file<<data[i*nVar+ivar]<<' ';
        }
        file<<'\n';
        std::cout<<"t=  "<<t<<"  error:  "<<error/n<<'\n';
        
    }


    for (auto it = data.begin(); it != data.end() ; ++it) std::cout << *it << " ";


    return 0;
    
}
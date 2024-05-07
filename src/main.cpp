#include <iostream>
#include "SpaceDis.hpp"
#include "macro.hpp"
#include "helloConfig.h"
#include "burgers.hpp"
#include "data.hpp"
int main(){
    
    //std::fstream file("data.txt",std::ios::out);
    real dt=0.001,t=0.0,tend=4-dt/4;
    ind nstep=floor((tend-t)/dt);
    ind n=200;
    ind nVar=1;

    Data dat;
    dat.solInit(n,nVar);

    Data coor;
    coor.init(n,1);
    coor.uniMesh();

    std::string name="CoordinateX";

    coor.cgnsoutputInit();
    dat.oneDsolOutput(0);


    /*file<<n<<' '<<nstep+1<<'\n';
    coor.output(&file,0);
    dat.output(&file,0);*/

    SpaceDis discrete(&dat,&coor,n,nVar);
    discrete.setMethod(WCNSJS5,BURGERS);
    ind step=0;
    while(t<tend)
    {
        Data tempdata=dat;

        //third order RK
        //stage 1
        std::vector<real>rhs=discrete.difference();
        for(ind i=0;i<n;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            dat[i*nVar+ivar]=tempdata[i*nVar+ivar]-dt*rhs[i*nVar+ivar];
        }

        //stage 2
        rhs=discrete.difference();
        for(ind i=0;i<n;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            dat[i*nVar+ivar]=0.75*tempdata[i*nVar+ivar]
                             -0.25*dt*rhs[i*nVar+ivar]
                             +0.25*dat[i*nVar+ivar];
        }

        //stage 3
        rhs=discrete.difference();
        for(ind i=0;i<n;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            dat[i*nVar+ivar]=1.0/3.0*tempdata[i*nVar+ivar]
                             -2.0/3.0*dt*rhs[i*nVar+ivar]
                             +2.0/3.0*dat[i*nVar+ivar];
        }

        t=t+dt;
        real error=0.0;
        for(ind i=0;i<n;i++)
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            real exactSol=0;
            real err=abs(dat[i*nVar+ivar]-exactSol);
            error+=err;

            
        }
        std::cout<<"t=  "<<t<<"  error:  "<<error/n<<'\n';
        step++;
        if(step%100==0)for(ind ivar=0;ivar<nVar;ivar++) dat.oneDsolOutput(t);

        
    }


    //for (auto it = data.begin(); it != data.end() ; ++it) std::cout << *it << " ";


    return 0;
    
}
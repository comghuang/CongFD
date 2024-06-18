
#include "zone.hpp"
#include "initializer.hpp"


int main()
{
    /*
    real dt=0.001,t=0.0,tend=0.5-dt/4;
    ind nstep=floor((tend-t)/dt);
    ind n=200;
    ind nVar=3;
    FluxType ftype=EULER1D;
    Zone zone;
    zone.init(n,1,1,1,ftype);
    ind step=0;
    while(t<tend)
    {
        zone.RK3(dt);

        t=t+dt;
        std::cout<<"t=  "<<t<<'\n';
        step++;
        if(step%10==0) zone.oneDsolOutput(t);
        
    }
    */
   std::array<int,3> iMax{100,2,2};
   std::array<double,6> calZone{-1,1,0,2,0,2};
   std::shared_ptr<Block> grid(new Block);
   grid->initUniform(iMax,calZone);
   grid->outputCgns();
   Info info;
   std::shared_ptr<Zone> zone(new Zone);
   zone->init(info,grid);
   Initializer init;
   init.solInit(grid,zone->getCons(),info);

   std::cout<<"Finish\n";
}
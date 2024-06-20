
#include "blockSolver.hpp"


int main()
{
    /*
    real dt=0.001,t=0.0,tend=0.5-dt/4;
    ind nstep=floor((tend-t)/dt);
    ind n=200;
    ind nVar=3;
    EquationType ftype=EULER1D;
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
//    std::shared_ptr<Info> info(new Info);
//    std::shared_ptr<Block> block(new Block);
//    std::shared_ptr<Initializer> initer(new Initializer(info));
//    std::shared_ptr<Equation> eqn(new Equation);
//    std::shared_ptr<Bnds> bnds(new Bnds);
//    std::shared_ptr<SpDistributor> spDis(new SpDistributor);

//    initer->initUniformBlock(block);
//    initer->initEqution(eqn,block);
//    initer->initBnds(bnds,eqn);
//    initer->initSpDistributor(spDis,eqn,block,bnds);

//     eqn->consToPrim();
//     bnds->update();
//     spDis->rhsSolve();
    BlockSolver bSolver;
    bSolver.solve();


   std::cout<<"Finish\n";
}
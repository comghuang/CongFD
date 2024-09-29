
#include "blockSolver.hpp"
#include "eigenSystem.hpp"


int main()
{
    // std::array<real,4> prim0={1.0,0.75,-0.5,1.0};
    // std::array<real,4> prim={2.0,-0.75,0.5,1.0};
    // eigensystemEuler2D eig=eigensystemEuler2D(prim0,{1,0,0});
    // auto eigValues=eig.primToChar(prim);
    // auto prim2=eig.charToPrim(eigValues);
    // std::cout<<"finish\n";

    omp_set_num_threads(10);

    Info* info=new Info;

    info->eqType=EULER;
    info->spMethod=WCNS5;


    info->diffMethod=MND6;
    info->interMethod=TCNS5;
    //info->interMethod=WCNSZ5Char;
    //info->BVD=true;
    //info->interMethod=WCNS5Char;
    //info->interMethod=WCNS5CONG;
    //info->interMethod=WCNSCONGPOLY;
    info->interMethod=WCNS5CONGZ;

    //Shu-Osher
    info->endStep=1;
    info->CFL=0.5;
    info->outputDt=1.8;
    info->nCase=1;
    info->calZone={0,10.0,0,0,0,0};
    info->iMax={201,2,2};
    info->dim=1;

    //sod tube
    info->CFL=0.5;
    info->endStep=20;
    info->outputDt=0.01;
    info->nCase=0;
    info->calZone={-0.5,0.5,0,0,0,0};
    info->iMax={101,2,2};
    info->dim=1;


    //lax sod tube
    // info->endStep=14;
    // info->outputDt=0.01;
    // info->CFL=0.1;
    // info->nCase=2;
    // info->calZone={-0.5,0.5,0,0,0,0};
    // info->iMax={2001,2,2};
    // info->dim=1;

    //lax sod tube speed test
    // info->endStep=14;
    // info->outputDt=0.01;
    // info->CFL=0.1;
    // info->nCase=2;
    // info->calZone={-0.5,0.5,0,0,0,0};
    // info->iMax={2001,2,2};
    // info->dim=1;

    //sedov
    info->endStep=1;
    info->outputDt=0.001;
    info->CFL=0.5;
    info->nCase=3;
    info->calZone={-2,2,0,0,0,0};
    info->iMax={400,2,2};
    info->dim=1;

    //Woodward-Colella
    // info->endStep=38;
    // info->outputDt=0.001;
    // info->CFL=0.5;
    // info->nCase=4;
    // info->calZone={0,1,0,0,0,0};
    // info->iMax={401,2,2};
    // info->dim=1;

    //双稀疏波
    // info->endStep=100;
    // info->outputDt=0.01;
    // info->CFL=0.5;
    // info->nCase=5;
    // info->calZone={-5,5,0,0,0,0};
    // info->iMax={401,2,2};
    // info->dim=1;

    //implosion
    // info->dt=0.0001;
    // info->endStep=25;
    // info->outputDt=0.1;
    // info->CFL=0.5;
    // info->nCase=2;
    // info->calZone={0,0.3,0,0.3,0,0};
    // info->iMax={201,201,2};
    // info->dim=2;

    //Riemann 1
    // real endt=0.3;
    // int outputsteps=1;
    // info->endStep=outputsteps;
    // info->outputDt=endt/outputsteps;
    // info->CFL=0.5;
    // info->nCase=0;
    // info->calZone={-0.5,0.5,-0.5,0.5,0,0};
    // info->iMax={401,401,2};
    // info->dim=2;

    //Riemann 2 vortex
    info->endStep=1;
    info->outputDt=0.3;
    info->CFL=0.5;
    info->nCase=1;
    info->calZone={-0.5,0.5,-0.5,0.5,0,0};
    info->iMax={401,401,2};
    info->dim=2;

    //RT instability
    //记得改GAMMA
    // info->endStep=1;
    // info->outputDt=1.95;
    // info->CFL=0.5;
    // info->nCase=3;
    // info->calZone={0,0.25,0,1,0,0};
    // info->iMax={101,401,2};
    // info->iMax={65,257,2};
    // info->dim=2;
    // info->sourceType=GRAVITY;


    //info->diffMethod=HDS6;
    //Double Mach
    // info->endStep=20;
    // info->outputDt=0.01;
    // info->CFL=0.5;
    // info->nCase=4;
    // info->calZone={0,4,0,1,0,0};
    // info->iMax={801,201,2};
    // info->dim=2;

    BlockSolver bSolver(info);
    auto start = std::chrono::high_resolution_clock::now(); 
    if(info->eqType!=EULER)
    bSolver.stepsLoop();
    else
    bSolver.stepsLoopCFL();
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count(); 
    //bSolver.stepsLoopDTS();
    //bSolver.solve();
    //bSolver.outputPrim();
    //bSolver.Test();
    
    
    

   std::cout<<"totaltime= "<<duration<<"   Finish\n";
   std::cout<<"time= "<<timepp/1e6<<"   Finish\n";
   std::cout<<"timesteps= "<<bSolver.timesteps<<"   Finish\n";
   std::cout<<"solvertime= "<<timesss<<'\n';
}
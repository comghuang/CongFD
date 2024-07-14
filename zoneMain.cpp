
#include "blockSolver.hpp"


int main()
{
    omp_set_num_threads(18);

    Info* info=new Info;

    info->eqType=EULER;
    info->spMethod=WCNS5;


    info->diffMethod=MND6;
    //info->interMethod=TCNS5;
    // info->interMethod=WCNS5Char;
    info->interMethod=WCNS5CONGABS;


    //Shu-Osher
    // info->endStep=18;
    // info->CFL=0.1;
    // info->outputDt=0.1;
    // info->nCase=1;
    // info->calZone={0,10.0,0,0,0,0};
    // info->iMax={201,2,2};
    // info->dim=1;

    //sod tube
    // info->dt=0.0001;
    // info->CFL=0.5;
    // info->endStep=20;
    // info->outputInterval=1000;
    // info->outputDt=0.01;
    // info->nCase=0;
    // info->calZone={-0.5,0.5,0,0,0,0};
    // info->iMax={101,2,2};
    // info->dim=1;

    //lax sod tube
    // info->dt=0.0001;
    // info->endStep=14;
    // info->outputDt=0.01;
    // info->CFL=0.1;
    // info->nCase=2;
    // info->calZone={-0.5,0.5,0,0,0,0};
    // info->iMax={101,2,2};
    // info->dim=1;

    //implosion
    // info->dt=0.0001;
    // info->endStep=250;
    // info->outputDt=0.01;
    // info->CFL=0.5;
    // info->nCase=2;
    // info->calZone={0,0.3,0,0.3,0,0};
    // info->iMax={201,201,2};
    // info->dim=2;

    //Riemann 1
    // info->endStep=60;
    // info->outputDt=0.01;
    // info->CFL=0.5;
    // info->nCase=0;
    // info->calZone={-0.5,0.5,-0.5,0.5,0,0};
    // info->iMax={401,401,2};
    // info->dim=2;

    //RT instability
    //记得改GAMMA
    // info->endStep=200;
    // info->outputDt=0.01;
    // info->CFL=0.5;
    // info->nCase=3;
    // info->calZone={0,0.25,0,1,0,0};
    // info->iMax={101,401,2};
    // info->dim=2;
    // info->sourceType=GRAVITY;


    //info->diffMethod=HDS6;
    //Double Mach
    //记得改GAMMA
    
    info->endStep=200;
    info->outputDt=0.001;
    info->CFL=0.5;
    info->nCase=4;
    info->calZone={0,4,0,1,0,0};
    info->iMax={801,201,2};
    info->dim=2;

    BlockSolver bSolver(info);
    if(info->eqType!=EULER)
    bSolver.stepsLoop();
    else
    bSolver.stepsLoopCFL();
    //bSolver.stepsLoopDTS();
    //bSolver.solve();
    //bSolver.outputPrim();
    //bSolver.Test();
    
    
    


   std::cout<<"Finish\n";
}
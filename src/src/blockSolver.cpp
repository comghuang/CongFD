#include "blockSolver.hpp"

BlockSolver::BlockSolver()
{


   info=new Info();
   block=new Block();
   initer=new Initializer(info);
   eqn=new Equation();
   bnds=new Bnds();
   spDis=new SpDistributor();

   initer->initUniformBlock(block);
   initer->initEqution(eqn,block);
   initer->initBnds(bnds,eqn,block->getICMax());
   initer->initSpDistributor(spDis,eqn,block,bnds);

   cons=eqn->getCons();
   rhs=eqn->getRhs();
}

BlockSolver::BlockSolver(Info info_)
{
   (*info)=info_;
   block=new Block();
   initer=new Initializer(info);
   eqn=new Equation();
   bnds=new Bnds();
   spDis=new SpDistributor();

   initer->initUniformBlock(block);
   cgnsIO.BlockCgnsOutput(block,info);
   initer->initEqution(eqn,block);
   initer->initBnds(bnds,eqn,block->getICMax());
   initer->initSpDistributor(spDis,eqn,block,bnds);
}

void BlockSolver::RK3()
{
    
    Data tempdata(*cons);
    int n=cons->size();
    double dt=info->dt;

    //third order RK
    //stage 1
    rhs->setZeros();
    eqn->consToPrim();
    bnds->update();
    spDis->rhsSolve();
    omp_set_num_threads(20);
    #pragma omp parallel for
    for(int i=0;i<n;i++)
    {
        (*cons)[i]=tempdata[i]-dt*(*rhs)[i];
    }

    //stage 2
    rhs->setZeros();
    eqn->consToPrim();
    bnds->update();
    spDis->rhsSolve();
    #pragma omp parallel for
    for(int i=0;i<n;i++)
    {
        (*cons)[i]=0.75*tempdata[i]-0.25*dt*(*rhs)[i]+0.25*(*cons)[i];
    }

    //stage 3
    rhs->setZeros();
    eqn->consToPrim();
    bnds->update();
    spDis->rhsSolve();
    #pragma omp parallel for
    for(int i=0;i<n;i++)
    {
        (*cons)[i]=1.0/3.0*tempdata[i]-2.0/3.0*dt*(*rhs)[i]+2.0/3.0*(*cons)[i];
    }
}

void BlockSolver::solve()
{
    RK3();
    info->t+=info->dt;
}

BlockSolver::~BlockSolver()
{
    delete info;
    delete block;
    delete initer;
    delete eqn;
    delete bnds;
    delete spDis;
}

void BlockSolver::outputGrid()
{
    cgnsIO.BlockCgnsOutput(block,info);
}

void BlockSolver::outputCons()
{
    cgnsIO.solCgnsOutput(eqn->getCons(),info);
}
void BlockSolver::outputPrim()
{
    eqn->consToPrim();
    cgnsIO.solCgnsOutput(eqn->getPrim(),info);
}

void BlockSolver::stepsLoop()
{
    
    for(;info->step<info->endStep;info->step++)
    {
        if(info->step%info->outputInterval==0)
        {
            outputGrid();
            outputPrim();
        }
        solve();
        std::cout<<std::format("time= {:.4f} \n",info->t);
    }
    outputGrid();
    outputPrim();
}
void BlockSolver::Test()
{
    outputGrid();
    RK3();

    //outputCons();
    //outputPrim();
    cgnsIO.solCgnsOutput(rhs,info);

}
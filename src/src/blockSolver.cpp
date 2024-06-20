#include "blockSolver.hpp"

BlockSolver::BlockSolver()
{


   info=std::make_shared<Info>();
   block=std::make_shared<Block>();
   initer=std::make_shared<Initializer>(info);
   eqn=std::make_shared<Equation>();
   bnds=std::make_shared<Bnds>();
   spDis=std::make_shared<SpDistributor>();

   initer->initUniformBlock(block);
   initer->initEqution(eqn,block);
   initer->initBnds(bnds,eqn,block->getICMax());
   initer->initSpDistributor(spDis,eqn,block,bnds);

   cons=eqn->getCons();
   rhs=eqn->getRhs();
}

BlockSolver::BlockSolver(std::shared_ptr<Info> info_)
{
   info=info_;
   block=std::make_shared<Block>();
   initer=std::make_shared<Initializer>(info);
   eqn=std::make_shared<Equation>();
   bnds=std::make_shared<Bnds>();
   spDis=std::make_shared<SpDistributor>();

   initer->initUniformBlock(block);
   initer->initEqution(eqn,block);
   initer->initBnds(bnds,eqn,block->getICMax());
   initer->initSpDistributor(spDis,eqn,block,bnds);
}

void BlockSolver::RK3()
{
    
    Data tempdata=*(cons);
    int n=cons->size();

    //third order RK
    //stage 1
    rhs->setZeros();
    eqn->consToPrim();
    bnds->update();
    spDis->rhsSolve();
    for(ind i=0;i<n;i++)
    {
        (*cons)[i]=tempdata[i]-dt*(*rhs)[i];
    }

    //stage 2
    rhs->setZeros();
    eqn->consToPrim();
    bnds->update();
    spDis->rhsSolve();
    for(ind i=0;i<n;i++)
    {
        (*cons)[i]=0.75*tempdata[i]-0.25*dt*(*rhs)[i]+0.25*(*cons)[i];
    }

    //stage 3
    rhs->setZeros();
    eqn->consToPrim();
    bnds->update();
    spDis->rhsSolve();
    for(ind i=0;i<n;i++)
    {
        (*cons)[i]=1.0/3.0*tempdata[i]-2.0/3.0*dt*(*rhs)[i]+2.0/3.0*(*cons)[i];
    }
}

void BlockSolver::solve()
{
    RK3();
}


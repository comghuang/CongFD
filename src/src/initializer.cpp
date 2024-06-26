#include"initializer.hpp"

void Initializer::solInit(Block* grid,Data* sol)
{
    std::vector<real> tempsol;
    switch (info->eqType)
    {
        /*case begin*/
    case LINEARCONV1D:
        if (grid->dim!=1)
        { std::cout<<"initialize: dim error \n";return;}
        switch (info->nCase)
        {
        case 0:
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                tempsol.push_back(1+sin(x));
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        
        default:
            break;
        }
        break;
    /*case end*/
    
    /*case begin*/
    case BURGERS1D:
        if (grid->dim!=1)
        { std::cout<<"initialize: dim error \n";return;}
        switch (info->nCase)
        {
        /*case 0 begin*/
        case 0:
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                tempsol.push_back(-sin(M_PI*x));
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        /*case 0 end*/
        
        default:
            break;
        }
        break;
    /*case end*/

    case EULER:
        /*case begin*/
        if (grid->dim==1)
        switch (info->nCase)
        {
        /*case 0 begin*/
        case 0:
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                real gamma=1.4;
                if (x<0)
                {
                    tempsol.push_back(1);
                    tempsol.push_back(0);
                    tempsol.push_back(1.0/(gamma-1)*1);
                }
                else
                {
                    tempsol.push_back(0.125);
                    tempsol.push_back(0);
                    tempsol.push_back(1.0/(gamma-1)*0.1);
                }
                
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        /*case 0 end*/
        default:
            break;
        }
        else if (grid->dim==2)
        switch (info->nCase)
        {
        case 0:
            //2D Riemann Problem;
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                if (x>0.25)
                {
                    if (y>0.25)
                    {
                        tempsol.push_back(1.5);
                        tempsol.push_back(0);
                        tempsol.push_back(0);
                        tempsol.push_back(1.0/(gamma-1)*1.5);
                    }
                    else
                    {
                        tempsol.push_back(0.5323);
                        tempsol.push_back(0);
                        tempsol.push_back(0.5323*1.206);
                        tempsol.push_back(1.0/(gamma-1)*0.3+0.5323*1.206*1.206/2);
                    }
                }
                else
                {
                    if (y>0.25)
                    {
                        tempsol.push_back(0.5323);
                        tempsol.push_back(0.5323*1.206);
                        tempsol.push_back(0);
                        tempsol.push_back(1.0/(gamma-1)*0.3+0.5323*1.206*1.206/2);
                    }
                    else
                    {
                        tempsol.push_back(0.138);
                        tempsol.push_back(0.138*1.206);
                        tempsol.push_back(0.138*1.206);
                        tempsol.push_back(1.0/(gamma-1)*0.029+0.138*1.206*1.206);
                    }
                }
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;
        
        default:
            break;
        }
        
    /*case end*/
    default:
        break;
    }
}


void Initializer::initUniformBlock(Block* block)
{
    if(!block)
    {
        std::cout<<"Initializer error: empty shared_ptr Block\n";
        return;
    }
    auto iMax=info->iMax;
    int dim=info->dim();
    auto icMax=info->icMax();
    int nVer=1,nCel=1;
    for (int i = 0; i < 3; i++)
    {
        nVer*=iMax[i];
        nCel*=icMax[i];
    }

    block->dim=dim;
    block->nVer=nVer;
    block->nCel=nCel;
    block->icMax=icMax;
    block->iMax=iMax;
    block->coorVer.init(nVer,dim);
    block->coorCel.init(nCel,dim);
    block->inited=true;

    //for vertex
    for (int idim = 0; idim < dim; idim++)
    {
        double cmin=info->calZone[idim*2];
        double cmax=info->calZone[idim*2+1];
        double interval=(cmax-cmin)/(iMax[idim]-1);

        int l,m,n;
        int* onedIndex=((idim == 0 ) ? &l : ( idim == 1 ? &m : &n));
        
        for (l = 0; l < iMax[0]; l++)
        for (m = 0; m < iMax[1]; m++)
        for (n = 0; n < iMax[2]; n++)
        {
            int globalIndex=l+m*iMax[0]+n*iMax[0]*iMax[1];
            block->coorVer(globalIndex,idim)=cmin+(*onedIndex)*interval;
        }
    }
    
    //for cellcenter
    for (int idim = 0; idim < dim; idim++)
    {
        int l,m,n,iLen=(dim==1?2:dim==2? 4:8);
        std::vector<int> index;
        index.resize(iLen);
        for (l = 0; l < icMax[0]; l++)
        for (m = 0; m < icMax[1]; m++)
        for (n = 0; n < icMax[2]; n++)
        {
            int iVerGlobal=l+m*iMax[0]+n*iMax[0]*iMax[1];
            int iCelGlobal=l+m*icMax[0]+n*icMax[0]*icMax[1];
            double temp=0;
            index[0]=iVerGlobal;
            index[1]=iVerGlobal+1;
            if (dim>=2)
            {
                index[2]=iVerGlobal+iMax[0];
                index[3]=iVerGlobal+iMax[0]+1;
            }
            if (dim>=3)
            {
                index[4]=iVerGlobal+iMax[0]*iMax[1];
                index[5]=iVerGlobal+1+iMax[0]*iMax[1];
                index[6]=iVerGlobal+iMax[0]+iMax[0]*iMax[1];
                index[7]=iVerGlobal+iMax[0]+1+iMax[0]*iMax[1];
            }

            for(auto iver:index) 
            temp+=block->coorVer(iver,idim);
            
            
            block->coorCel(iCelGlobal,idim)=temp/iLen;
            
        }
    }
}


Initializer::Initializer()
{

}

Initializer::Initializer(Info* info_)
{
    info=info_;
}


void Initializer::initEqution(Equation* eq,Block* block)
{
    if(!eq||!block)
    {
        std::cout<<"Initializer error: empty shared_ptr Equation or Block\n";
        return;
    }
    eq->n=block->icMax[0]*block->icMax[1]*block->icMax[2];
    eq->nPrim=info->nPrim();
    eq->nCons=info->nCons();
    eq->type=info->eqType;
    eq->dim=block->dim;

    
    eq->rhs=new Data(eq->n,eq->nCons);
    eq->cons=new Data(eq->n,eq->nCons);
    eq->prim=new Data(eq->n,eq->nPrim);
    eq->inited=true;
    eq->cons->setvarName(info->getVarNameListCons());
    eq->prim->setvarName(info->getVarNameListPrim());
    eq->rhs->setvarName(info->getVarNameListRhs());

    solInit(block,eq->cons);
}


void Initializer::initBnds(Bnds* bnds,Equation* eqn,std::array<int,3> iMax)
{
    bnds->iMax=iMax;
    bnds->dim=info->dim();

    int nGhost=info->nGhostCell();
    int nCons=info->nCons();
    int nPrim=info->nPrim();

    std::array<int,3> nBnds  {bnds->iMax[1]*bnds->iMax[2],
                            bnds->iMax[0]*bnds->iMax[2],
                            bnds->iMax[0]*bnds->iMax[1]};
    int nBnd=0;
    for (int i = 0; i < bnds->dim; i++)
    {
        nBnd+=nBnds[i]*2;
    }
    
    

    bnds->oneDBnds.resize(nBnd);
    std::array<int,2> offsets;
    switch (info->eqType)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        {
            for (int i = 0; i < 2; i++)
            {
                bnds->oneDBnds.at(i)=std::make_shared<OneDBnd>(nGhost,nPrim,PERIODIC1D);
            }
            offsets=calOffsetInverse(1,0,0,bnds->iMax);
            bnds->oneDBnds.at(0)->setUpdate(eqn->prim,offsets[0],offsets[1]);

            offsets=calOffset(1,0,0,bnds->iMax);
            bnds->oneDBnds.at(1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
        }
        break;
    case EULER:
        if(eqn->dim==1)
        {
            bnds->oneDBnds.at(0)=std::make_shared<OneDBnd>(nGhost,nPrim,SUPERSONICOUTLET);
            offsets=calOffset(1,0,0,bnds->iMax);
            bnds->oneDBnds.at(0)->setUpdate(eqn->prim,offsets[0],offsets[1]);


            bnds->oneDBnds.at(1)=std::make_shared<OneDBnd>(nGhost,nPrim,SUPERSONICOUTLET);
            offsets=calOffsetInverse(1,0,0,bnds->iMax);
            bnds->oneDBnds.at(1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
        }
        else if (eqn->dim==2)
        {
            for (int i = 0; i < iMax[1]; i++)
            {
                bnds->oneDBnds.at(2*i)=std::make_shared<OneDBnd>(nGhost,nPrim,SUPERSONICOUTLET);
                offsets=calOffset(1,i,0,bnds->iMax);
                bnds->oneDBnds.at(2*i)->setUpdate(eqn->prim,offsets[0],offsets[1]);

                bnds->oneDBnds.at(2*i+1)=std::make_shared<OneDBnd>(nGhost,nPrim,SUPERSONICOUTLET);
                offsets=calOffsetInverse(1,i,0,bnds->iMax);
                bnds->oneDBnds.at(2*i+1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
            }

            for (int i = 0; i < iMax[0]; i++)
            {
                bnds->oneDBnds.at(2*i+iMax[1]*2)=std::make_shared<OneDBnd>(nGhost,nPrim,SUPERSONICOUTLET);
                offsets=calOffset(2,i,0,bnds->iMax);
                bnds->oneDBnds.at(2*i+iMax[1]*2)->setUpdate(eqn->prim,offsets[0],offsets[1]);

                bnds->oneDBnds.at(2*i+1+iMax[1]*2)=std::make_shared<OneDBnd>(nGhost,nPrim,SUPERSONICOUTLET);
                offsets=calOffsetInverse(2,i,0,bnds->iMax);
                bnds->oneDBnds.at(2*i+1+iMax[1]*2)->setUpdate(eqn->prim,offsets[0],offsets[1]);
            }
            
        }
        

        break;
    
    default:
        break;
    }
    
}

void Initializer::initSpDistributor(SpDistributor* spDis,Equation* eqn
                                    ,Block* block,Bnds* bnds)
{
    
    spDis->nCons=info->nCons();
    spDis->nPrim=info->nPrim();
    spDis->iMax=block->icMax;
    spDis->prim=eqn->prim;
    spDis->cons=eqn->cons;
    spDis->dim=info->dim();
    spDis->bnds=bnds;
    spDis->rhs=eqn->rhs;
    spDis->info=info;
    

}
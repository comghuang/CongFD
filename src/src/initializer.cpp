#include"initializer.hpp"

void Initializer::solInit(std::shared_ptr<Block> grid,std::shared_ptr<Data> sol)
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
        { std::cout<<"dim error \n";return;}
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

    case EULER1D:
        if (grid->dim!=1)
        { std::cout<<"dim error \n";return;}
        /*case begin*/
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
    /*case end*/
    default:
        break;
    }
}


void Initializer::initUniformBlock(std::shared_ptr<Block> block)
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

Initializer::Initializer(std::shared_ptr<Info> info_)
{
    info=info_;
}


void Initializer::initEqution(std::shared_ptr<Equation> eq,std::shared_ptr<Block> block)
{
    if(!eq||!block)
    {
        std::cout<<"Initializer error: empty shared_ptr Equation or Block\n";
        return;
    }
    eq->n=block->icMax[0]*block->icMax[1]*block->icMax[2];
    eq->nPrim=info->nPrim();
    eq->nCons=info->nCons();
    
    eq->rhs=std::make_shared<Data>(eq->n,eq->nCons);
    eq->cons=std::make_shared<Data>(eq->n,eq->nCons);
    eq->prim=std::make_shared<Data>(eq->n,eq->nPrim);
    eq->inited=true;

    solInit(block,eq->cons);
}


void Initializer::initBnds(std::shared_ptr<Bnds> bnds)
{
    bnds->iMax=info->icMax();
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

    switch (info->eqType)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        for (int i = 0; i < 2; i++)
        {
            bnds->oneDBnds.at(i)=std::make_shared<OneDBnd>(nGhost,nPrim,PERIODIC1D);
        }
        break;
    case EULER1D:
        bnds->oneDBnds.at(0)=std::make_shared<OneDBnd>(nGhost,nPrim,DIRICLET_SODL);
        bnds->oneDBnds.at(1)=std::make_shared<OneDBnd>(nGhost,nPrim,DIRICLET_SODR);
        break;
    
    default:
        break;
    }
    
}

void Initializer::initSpDistributor(std::shared_ptr<SpDistributor> spDis,std::shared_ptr<Equation> eqn
                                    ,std::shared_ptr<Block> block,std::shared_ptr<Bnds> bnds)
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
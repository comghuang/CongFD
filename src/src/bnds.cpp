#include"bnds.hpp"

void Bnds::initFromCode(std::shared_ptr<Block> block,Info info)
{
    iMax=block->icMax;
    dim=block->dim;

    int nGhost=info.nGhostCell();
    int nCons=info.nCons();
    int nPrim=info.nPrim();



    if(info.eqType==LINEARCONV1D)
    {
        switch (info.nCase)
        {
        case 0:
            {
                bnds.reserve(2);
                for (int i = 0; i < 2; i++)
                {
                    bnds.at(i)=std::make_shared<OneDBnd>();
                    bnds.at(i)->init(nGhost,nPrim,PERIODIC1D);
                }
                
            }
            break;
        
        default:
            break;
        }
    }
    else if(info.eqType==BURGERS1D)
    {
        switch (info.nCase)
        {
        case 0:
            /* code */
            break;
        
        default:
            break;
        }
    }
    else if(info.eqType==EULER1D)
    {
        switch (info.nCase)
        {
        case 0:
            /* code */
            break;
        
        default:
            break;
        }
    }
}
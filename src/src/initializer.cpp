#include"initializer.hpp"

void Initializer::solInit(std::shared_ptr<Block> grid,std::shared_ptr<Data> sol,Info info)
{
    std::vector<real> tempsol;
    switch (info.fType)
    {
        /*case begin*/
    case LINEARCONV1D:
        if (grid->dim!=1)
        { std::cout<<"initialize: dim error \n";return;}
        switch (info.nCase)
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
        switch (info.nCase)
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
        switch (info.nCase)
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
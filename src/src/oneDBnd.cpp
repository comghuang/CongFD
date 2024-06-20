#include "oneDBnd.hpp"

OneDBnd::OneDBnd(ind n_,ind nVar_,BndType bType_)
{
    n=n_;
    nVar=nVar_;
    data.resize(n*nVar,0.0);
    type=bType_;
}
real& OneDBnd::operator()(ind i,ind ivar)
{
    return data[i*nVar+ivar];
}

ind OneDBnd::getN()
{
    return n;
}


void OneDBnd::setValue(std::vector<real> value)
{
    if(value.size()!=data.size())
    {
        std::cout<<"OneDBnd error: setValue() incorrect size\n";
    }
    std::copy(value.begin(),value.end(),data.begin());
}
BndType OneDBnd::getType()
{
    return type;
}

void OneDBnd::setUpdate(std::shared_ptr<Data> prim_,int i0_,int offset_)
{
    prim=prim_;
    i0=i0_;
    offset=offset_;
}

void OneDBnd::update()
{
    switch (type)
    {
    case PERIODIC1D:
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nVar; j++)
            {
                data.at(i*nVar+j)=(*prim)(i0+i*offset,j);
            }
            
        }
        break;
    case SUPERSONICOUTLET:
        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < nVar; j++)
            {
                data.at(i*nVar+j)=(*prim)(i0,j);
            }
        }
        break;
    
    default:
        break;
    }
    
    
}
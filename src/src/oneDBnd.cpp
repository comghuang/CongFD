#include "oneDBnd.hpp"


real& OneDBnd::operator()(ind i,ind ivar)
{
    return data[i*nVar+ivar];
}

ind OneDBnd::getN()
{
    return n;
}


void OneDBnd::init(ind n_,ind nVar_,BndType bType_)
{
    n=n_;
    nVar=nVar_;
    data.resize(n*nVar,0.0);
    type=bType_;
}

void OneDBnd::setGhostValue(real* value)
{
    for(ind i=0;i<n;i++)
    {
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            (*this)(i,ivar)=value[i*nVar+ivar];
        }
    }
}
BndType OneDBnd::getType()
{
    return type;
}
#include "boundary.hpp"


real& OneDGhost::operator()(ind i,ind ivar)
{
    return data[i*nVar+ivar];
}

ind OneDGhost::getN()
{
    return n;
}


void OneDGhost::init(ind n_,ind nVar)
{
    n=n_;
    nVar=nVar;
    data.resize(n*nVar,0.0);
}

void OneDGhost::setGhostValue(std::vector<real> value)
{
    for(ind i=0;i<n;i++)
    {
        for(ind ivar=0;ivar<nVar;ivar++)
        {
            (*this)(i,ivar)=value[i*nVar+ivar];
        }
    }
}
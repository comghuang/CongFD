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
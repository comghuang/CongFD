#include "data.hpp"



void Data::init(ind n_,ind nvar_)
{
    n=n_;
    nVar=nvar_;
    data.resize(nVar*n,0.0);

    for(ind i=0;i<n;i++)
    {
        real h=2.0/n;
        real xi=h/2.0+i*h-1.0;
        data[i]=-sin(M_PI*xi);
    }
}

real& Data::operator() (ind i,ind j)
{
    if (i<0) return(this->data)[(n+i)*this->nVar+j];
    if (i>=n) return(this->data)[(i-n)*this->nVar+j];
    return (this->data)[i*this->nVar+j];
}

real& Data::operator[] (ind i)
{
    return (this->data)[i];
}

void Data::uniMesh()
{
    real h=2.0/n;
    for(ind i=0;i<n;i++)
    {
        real xi=h/2.0+i*h-1.0;
        (*this)(i,0)=xi;
    }
}

void Data::output(std::fstream* f,ind iVar)
{
    for(ind i=0;i<n;i++)
    {
        *f<<(*this)(i,iVar)<<' ';
    }
    *f<<'\n';
}


real Data::maxElement(ind ivar)
{
    if(data.empty()) return 0;
    real res=data[0];
    for(int i=1;i<n;i++)
    {
        res=(res<data[i]) ? data[i] : res;
    }
    return res;
}
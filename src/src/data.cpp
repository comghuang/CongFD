#include "data.hpp"

Data::Data(int n_,int nVar_)
{
    n=n_;
    nVar=nVar_;
    data.resize(nVar*n,0.0);
}

void Data::solInit(ind n_,ind nvar_)
{
    n=n_;
    nVar=nvar_;
    data.resize(nVar*n,0.0);

    for(ind i=0;i<n;i++)
    {
        real h=2.0/n;
        real xi=h/2.0+i*h-1.0;
        //data[i]=-sin(M_PI*xi);//for burgers equation

        //for sod tube 1D
        
        real gamma=GAMMA;
        if(xi<0)
        {
            (*this)(i,0)=1.0;
            (*this)(i,1)=0;
            (*this)(i,2)=1.0/(gamma-1)*1;
        }
        else
        {
            (*this)(i,0)=0.125;
            (*this)(i,1)=0;
            (*this)(i,2)=1.0/(gamma-1)*0.1;
        }
    }
}
void Data::init(ind n_,ind nvar_)
{
    n=n_;
    nVar=nvar_;
    data.resize(nVar*n,0.0);
}

real& Data::operator() (ind i,ind ivar)
{
    return data[i*nVar+ivar];
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



void Data::setValue(std::vector<real> value)
{
    if(value.size()!=n*nVar)
    {
        std::cout<<"vector length incorrect \n";
        return;
    }
    std::copy(value.begin(),value.end(),data.begin());
}

void Data::operator= (Data& dat)
{
    if(this->n==dat.n&&this->nVar==dat.nVar)
    {
        for(ind i = 0; i < n*nVar; i++)
        {
            (*this)[i]=dat[i];
        }
        
    }
    else
    {
        std::cout<<"incorrect size at class Data=Data\n";
    }
}

void Data::operator+= (std::vector<real> arr)
{
    if (arr.size()==n*nVar)
    {
        for (ind i = 0; i < n*nVar; i++)
        {
            (*this)[i]+=arr[i];
        }
        
    }
    else
    {
        std::cout<<"incorrect size at class Data += vector<real>\n";
    }
}

void Data::setZeros()
{
    std::fill(data.begin(),data.end(),0.0);
}



std::vector<real>::iterator Data::begin()
{
    return data.begin();
}
std::vector<real>::iterator Data::end()
{
    return data.end();
}

real Data::size()
{
    return data.size();
}
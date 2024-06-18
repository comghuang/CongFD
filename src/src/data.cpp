#include "data.hpp"



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
    if (i<0) 
    {
        return (*ghVertex[LEFTT])(-(i+1),ivar);
    }
    if (i>=n) 
    {
        return (*ghVertex[RIGHT])(i-n,ivar);
    }
    return data[(i0+offset*i)*nVar+ivar];
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

void Data::setGhostVertex(OneDBnd* bndl,OneDBnd* bndr)
{
    ghVertex[LEFTT]=bndl;
    ghVertex[RIGHT]=bndr;
}

void Data::updateGhostVertex()
{
    for (ind ilr = 0; ilr < 2; ilr++)
    {
        std::vector<real> res;
        switch (ghVertex[ilr]->getType())
        {
        case PERIODIC:
            res.reserve(ghVertex[ilr]->getN()*nVar);
            for (ind i=0;i<ghVertex[ilr]->getN();i++)
            {
                for(ind ivar=0;ivar<nVar;ivar++)
                {
                    if (ilr==LEFTT) res.push_back((*this)(n-1-i,ivar));
                    else res.push_back((*this)(i,ivar));
                }
            }
            ghVertex[ilr]->setGhostValue(&(res[0]));
            break;
        case DIRICLET_SODL:
            res.reserve(ghVertex[ilr]->getN()*nVar);
            for (ind i=0;i<ghVertex[ilr]->getN();i++)
            {
                res.push_back(1.0);
                res.push_back(0.0);
                res.push_back(1.0);
                res.push_back(GAMMA/(GAMMA-1));
                res.push_back(1.0);
            }
            ghVertex[ilr]->setGhostValue(&(res[0]));
            break;
        case DIRICLET_SODR:
            res.reserve(ghVertex[ilr]->getN()*nVar);
            for (ind i=0;i<ghVertex[ilr]->getN();i++)
            {
                res.push_back(0.125);
                res.push_back(0.0);
                res.push_back(0.1);
                res.push_back(GAMMA/(GAMMA-1)*0.8);
                res.push_back(0.8);
            }
            ghVertex[ilr]->setGhostValue(&(res[0]));
            break;
        default:
            break;
        }
    }

}


std::array<ind,2> Data::getNGhost()
{
    std::array<ind,2> res;
    res[LEFTT]=ghVertex[LEFTT]->getN();
    res[RIGHT]=ghVertex[RIGHT]->getN();
    return res;
}

void Data::setValue(real* value,ind len)
{
    if(len!=n*nVar)
    {
        std::cout<<"vector length incorrect \n";
        return;
    }
    for(int i=0;i<len;i++)
    {
        data[i]=value[i];
    }
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
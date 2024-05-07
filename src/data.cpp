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
        data[i]=-sin(M_PI*xi);//for burgers equation

        //for sod tube 1D
        /*
        real gamma=GAMMA;
        if(xi<0)
        {
            (*this)(i,0)=1;
            (*this)(i,1)=0;
            (*this)(i,2)=1.0/(gamma-1)*1;
        }
        else
        {
            (*this)(i,0)=0.125;
            (*this)(i,1)=0;
            (*this)(i,2)=1.0/(gamma-1)*0.1;
        }*/
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
    if (i<0) return ghVertex[LEFTT](1-i,ivar);
    if (i>=n) return ghVertex[RIGHT](n-i,ivar);
    return (this->data)[i*this->nVar+ivar];
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

void Data::setGhostVertex(ind nL,ind igh)
{
    ghVertex[igh].init(nL,nVar);
}

void Data::updateGhostVertex()
{
    std::vector<real> res;
    for (ind i=0;i<ghVertex[LEFTT].getN();i++)
    {
        for(ind ivar;ivar<nVar;ivar++)
        {
            res.push_back((*this)(n-1-i,ivar));
        }
    }
    ghVertex[LEFTT].setGhostValue(res);
    res.clear();
    for (ind i=0;i<ghVertex[RIGHT].getN();i++)
    {
        for(ind ivar;ivar<nVar;ivar++)
        {
            res.push_back((*this)(i,ivar));
        }
    }
    ghVertex[RIGHT].setGhostValue(res);

}


std::array<ind,2> Data::getNGhost()
{
    std::array<ind,2> res;
    res[LEFTT]=ghVertex[LEFTT].getN();
    res[RIGHT]=ghVertex[RIGHT].getN();
}
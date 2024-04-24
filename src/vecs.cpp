#include"vecs.hpp"

void vecs::accept(std::vector<real>* data_,ind nVar_)
{
    this->data=data_;
    this->nVar=nVar_;
    this->n=(data->size())/nVar_;
    
}
void vecs::allocate(ind n_,ind nVar_)
{
    if (this->data!=NULL)delete this->data;
    this->data=new std::vector<real>;
    this->data->resize(n_*nVar_,0.0);
    n=n_;
    nVar=nVar_;
}
void vecs::free()
{
    delete this->data;
}

real& vecs::operator() (ind i,ind j)
{
    if (i<0) return(*this->data)[(n+i)*this->nVar+j];
    if (i>=n) return(*this->data)[(i-n)*this->nVar+j];
    return (*this->data)[i*this->nVar+j];
}

ind vecs::getSize()
{
    return n*nVar;
}
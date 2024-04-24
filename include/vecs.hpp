#pragma once 

#include"macro.hpp"
#include<vector>
class vecs
{
    public:
    
    void accept(std::vector<real>*,ind);
    void allocate(ind,ind);
    void free();
    real& operator() (ind,ind);
    ind getSize();

    private:
    std::vector<real>* data=NULL;
    ind nVar;
    ind n;
};
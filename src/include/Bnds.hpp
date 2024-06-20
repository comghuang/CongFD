#pragma once
#include "info.hpp"
#include "oneDBnd.hpp"





class Bnds
{
    public:
    std::array<std::shared_ptr<OneDBnd>,2> getOneDBnd(int,int,int);

    private:
    friend class Initializer;
    ind dim;
    std::array<ind,3> iMax;
    //idim*
    std::vector<std::shared_ptr<OneDBnd>> oneDBnds;//JK*2+IK*2+IJ*2
};
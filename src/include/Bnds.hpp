#pragma once
#include "block.hpp"
#include "info.hpp"





class Bnds
{
    public:
    void initFromCode(std::shared_ptr<Block>,Info);

    private:
    ind dim;
    std::array<ind,3> iMax;
    //idim*
    std::vector<std::shared_ptr<OneDBnd>> bnds;
};
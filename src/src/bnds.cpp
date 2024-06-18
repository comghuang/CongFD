#include"bnds.hpp"

void Bnds::initFromCode(std::shared_ptr<Block> block)
{
    iMax=block->icMax;
    dim=block->dim;
}
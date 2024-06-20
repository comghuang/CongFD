#pragma once
#include"block.hpp"
#include"info.hpp"
#include"equation.hpp"
#include"Bnds.hpp"
#include"sp_distributor.hpp"

class Initializer
{
    public:
    Initializer();
    Initializer(std::shared_ptr<Info>);
    void solInit(std::shared_ptr<Block>,std::shared_ptr<Data>);
    void initUniformBlock(std::shared_ptr<Block>);
    void initEqution(std::shared_ptr<Equation>,std::shared_ptr<Block>);
    void initBnds(std::shared_ptr<Bnds> bnds,std::shared_ptr<Equation>,std::array<int,3> iMax);
    void initSpDistributor(std::shared_ptr<SpDistributor>,std::shared_ptr<Equation>,std::shared_ptr<Block>,std::shared_ptr<Bnds>);

    private:
    std::shared_ptr<Info> info;
    
    
};
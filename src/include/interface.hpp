#pragma once
#include "macro.hpp"
#include "Data.hpp"
#include "oneDBnd.hpp"

class InterfaceDonor
{
    /*Cache for Periodical Boundary*/
    public:
    InterfaceDonor(int n);
    void getValue();
    void sendValue();
    
    private:
    std::vector<int> indexes;
    std::vector<real> data;
    std::vector<std::shared_ptr<OneDBnd>> bnds; 


};
class InterfaceAccept
{
    /*Cache for Periodical Boundary*/
    public:
    InterfaceAccept(int n,int nVar);
    InterfaceAccept();
    void getValue(std::vector<real>);
    
    private:
    int n,nVar;
    std::vector<std::shared_ptr<OneDBnd>> bnds; 


};
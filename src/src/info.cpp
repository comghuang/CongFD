#include "info.hpp"

int Info::nGhostCell()
{
    
    if(WCNSJS5==spMethod && TRAD6 == diffMethod) return 5;
    else if(WCNSJS5==spMethod && HDS6 == diffMethod) return 3;
    else if(MUSCL==spMethod && TRAD2 == diffMethod) return 1;
    else if(FIRSTORDER==spMethod && TRAD2 == diffMethod) return 1;

    else
    {
        std::cout<<"Info error: undifined spMethod and diffMethod combination\n";
        
    }
    return 0;
}
int Info::nFluxPoint()
{
    switch (eqType)
    {
    case TRAD2:
    case HDS6:
        return 0;
        break;
    case TRAD6:
        return 2;
        break;
    
    default:
        break;
    }
}

int Info::nPrim()
{
    switch (eqType)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        return 1;
        break;
    case EULER1D:
        return 5;
        break;
    
    default:
        std::cout<<"Info error: undifined eqType in nPrim()\n";
        return 0;
        
        break;
    }
}
int Info::nCons()
{
    switch (eqType)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        return 1;
        break;
    case EULER1D:
        return 3;
        break;
    
    default:
        std::cout<<"Info error: undifined eqType in nCons()\n";
        return 0;
        
        break;
    }
}
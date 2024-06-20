#include "info.hpp"

int Info::nGhostCell()
{
    
    if(WCNSJS5==spMethod && TRAD6 == diffMethod) return 5;
    else if(WCNSJS5==spMethod && HDS6 == diffMethod) return 3;
    else if(MUSCL==spMethod && TRAD2 == diffMethod) return 2;
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
        return 0;
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


int Info::dim()
{
    int dim=0;
    for (int i = 0; i < 3; i++)
    {
        if (iMax[i]>2) dim++;
    }
    return dim;
}

std::array<int,3> Info::icMax()
{
    std::array<int,3> icMax;
    for (int i = 0; i < 3; i++)
    {
        if (iMax[i]<2) icMax[i]=1;
        else
        icMax[i]=iMax[i]-1;
    }
    return icMax;
}


BndType Info::defaultBndType()
{
    switch (eqType)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        return PERIODIC1D;
        break;
    
    default:
        return TYPENULL;
        break;
    }
}
#include "solverType.hpp"
#include "differ.hpp"
#include "fluxPointFlux.hpp"
#include "fluxSchemes.hpp"
#include "interScheme.hpp"
#include "reconBVD.hpp"
#include "reconTenoThinc.hpp"
#include "reconstructor.hpp"
#include "reconstructorCellCenter.hpp"
#include "reconstructorFaceCenter.hpp"
#include "solvePointFlux.hpp"
SolverType::SolverType(const SolverType& origin)
{
    info = origin.info;
    reconer = origin.reconer;
    solvePointSolver = origin.solvePointSolver;
    fluxPointSolver = origin.fluxPointSolver;
    differ = origin.differ;
}

SolverType::SolverType(Info* info_)
{
    info = info_;
    initReconer();
    initFluxPointSolver();
    initSolvePointSolver();
    initDiffer();
}

void SolverType::initReconer()
{

    if (info->interMethod == TRIAL) {
        reconer = pro::make_proxy<ProxyDataManipulator>(
            // Recon5OrderCellCenter<teno5_CongCell>());
            // Recon5OrderFaceCenter<weno5_JSchen>());
            ReconerTenoThinc<3, eigensystemEuler1D>());
        return;
    }
    switch (info->eqType) {
    case LINEARCONV1D:
    case BURGERS1D:
        switch (info->interMethod) {
        case WCNS5:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<weno5_JSchen>());
            break;
        case WCNSZ5:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<weno5_Z>());
            break;
        case TCNS5:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<Teno5_Z>());
            break;
        case TCNS5CT4:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<Teno5_ZCT4>());
            break;
        case TCNS5CT7:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<Teno5_ZCT7>());
            break;
        case WCNS5CONGZ:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<Teno5_CongZ>());
            break;
        case WCNS5CONGZCT4:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<Teno5_CongZCT4>());
            break;
        case WCNS5CONGZCT7:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<Teno5_CongZCT7>());
            break;
        case BVD5:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                ReconerBVD());
            break;

        default:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<weno5_JSchen>());
            break;
        }
        break;
    case EULER:
        if (info->dim == 1) {
            switch (info->interMethod) {
            case WCNS5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<weno5_JSchen>());
                break;
            case WCNSZ5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<weno5_Z>());
                break;
            case TCNS5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<Teno5_Z>());
                break;
            case TCNS5CT4:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<Teno5_ZCT4>());
                break;
            case TCNS5CT7:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<Teno5_ZCT7>());
                break;
            case WCNS5CONGZ:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    // Recon5OrderFaceCenter<weno5_Z>());
                    // Recon5OrderFaceCenter<Teno5_CongZ>());
                    Recon5Order1DEulerEig<Teno5_CongZ>());
                break;
            case WCNS5CONGZCT4:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<Teno5_CongZCT4>());
                break;
            case WCNS5CONGZCT7:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<Teno5_CongZCT7>());
                break;
            case BVD5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    ReconerBVD());
                break;

            default:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<weno5_JSchen>());
                break;
            }
        } else if (info->dim == 2) {
            switch (info->interMethod) {
            case WCNS5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<weno5_JSchen>());
                break;
            case WCNSZ5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<weno5_Z>());
                break;
            case TCNS5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<Teno5_Z>());
                break;
            case TCNS5CT4:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<Teno5_ZCT4>());
                break;
            case TCNS5CT7:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<Teno5_ZCT7>());
                break;
            case WCNS5CONGZ:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<Teno5_CongZ>());
                // Recon5OrderFaceCenter<weno5_Z>());
                break;
            case WCNS5CONGZCT4:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<Teno5_CongZCT4>());
                break;
            case WCNS5CONGZCT7:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<Teno5_CongZCT7>());
                break;
            case BVD5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    ReconerBVD());
                break;
            default:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<weno5_JSchen>());
                break;
            }
        }
        break;
    default:
        break;
    }
}

void SolverType::initSolvePointSolver()
{
    switch (info->eqType) {
    case LINEARCONV1D:
        solvePointSolver = pro::make_proxy<ProxySolvePointSolver>(
            SolvePointSolver<1, fluxSolveLinearConv>());
        break;
    case BURGERS1D:
        solvePointSolver = pro::make_proxy<ProxySolvePointSolver>(
            SolvePointSolver<1, fluxSolveBurgers>());
        break;
    case EULER:
        if (info->dim == 1) {
            solvePointSolver = pro::make_proxy<ProxySolvePointSolver>(
                SolvePointSolver<3, fluxSolveEuler1D>());
        } else if (info->dim == 2) {
            solvePointSolver = pro::make_proxy<ProxySolvePointSolver>(
                SolvePointSolver<4, fluxSolveEuler2D>());
        }
        break;
    default:
        break;
    }
}

void SolverType::initFluxPointSolver()
{
    switch (info->eqType) {
    case EULER:
        switch (info->fluxMethod) {
        case ROE:
            if (info->dim == 1) {
                fluxPointSolver = pro::make_proxy<ProxyFluxPointSolver>(
                    FluxPointSolver<3, roeFlux1D2>());
            } else if (info->dim == 2) {
                fluxPointSolver = pro::make_proxy<ProxyFluxPointSolver>(
                    FluxPointSolver<4, roeFlux2DSym>());
            }
            break;
        case HLLC:
            if (info->dim == 1) {
                fluxPointSolver = pro::make_proxy<ProxyFluxPointSolver>(
                    FluxPointSolver<3, HLLCFlux1D>());
            } else if (info->dim == 2) {
                fluxPointSolver = pro::make_proxy<ProxyFluxPointSolver>(
                    FluxPointSolver<4, HLLCFlux2D2>());
            }
            break;

        default:
            break;
        }
        break;
    default:
        break;
    }
}

void SolverType::initDiffer()
{
    switch (info->eqType) {
    case EULER:
        switch (info->diffMethod) {
        case TRAD2:
            differ = pro::make_proxy<ProxyDiffer>(Differ());
            solvePointSolver->setNlNr(0, 0);
            break;
        case TRAD6:
            differ = pro::make_proxy<ProxyDiffer>(ExplicitDif6());
            solvePointSolver->setNlNr(0, 0);
            break;
        case MND6:
            differ = pro::make_proxy<ProxyDiffer>(MidNodeAndNodeDif6());
            solvePointSolver->setNlNr(1, 1);
            break;
        default:
            break;
        }
        break;
    default:
        break;
    }
}
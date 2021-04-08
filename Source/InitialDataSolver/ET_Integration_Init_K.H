#ifndef ET_INTEGRATION_INIT_K_H
#define ET_INTEGRATION_INIT_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_init(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab, 
        amrex::Real time, const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo(); 
        
        amrex::Real x_U_0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x_U_1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x_U_2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        amrex::Real bareMp = 0.48559999999999998;

        amrex::Real bareMm = 0.48559999999999998;

        amrex::Real cp_U_0 = 4.891;
        amrex::Real cp_U_1 = 0;
        amrex::Real cp_U_2 = 0;

        amrex::Real cm_U_0 = -4.891;
        amrex::Real cm_U_1 = 0;
        amrex::Real cm_U_2 = 0;

        amrex::Real Pp_U_0 = 0;
        amrex::Real Pp_U_1 = 0.0969;
        amrex::Real Pp_U_2 = 0;

        amrex::Real Pm_U_0 = 0;
        amrex::Real Pm_U_1 = -0.0969;
        amrex::Real Pm_U_2 = 0;

        amrex::Real Pp = std::sqrt(std::pow(Pp_U_0, 2) + std::pow(Pp_U_1, 2) + std::pow(Pp_U_2, 2));

        amrex::Real Pm = std::sqrt(std::pow(Pm_U_0, 2) + std::pow(Pm_U_1, 2) + std::pow(Pm_U_2, 2));

        amrex::Real xpc_U_0 = -cp_U_0 + x_U_0;
        amrex::Real xpc_U_1 = -cp_U_1 + x_U_1;
        amrex::Real xpc_U_2 = -cp_U_2 + x_U_2;

        amrex::Real xmc_U_0 = -cm_U_0 + x_U_0;
        amrex::Real xmc_U_1 = -cm_U_1 + x_U_1;
        amrex::Real xmc_U_2 = -cm_U_2 + x_U_2;

        amrex::Real rpc = std::max(0.00100000000000000, std::sqrt(std::pow(xpc_U_0, 2) + std::pow(xpc_U_1, 2) + std::pow(xpc_U_2, 2)));

        amrex::Real rmc = std::max(0.00100000000000000, std::sqrt(std::pow(xmc_U_0, 2) + std::pow(xmc_U_1, 2) + std::pow(xmc_U_2, 2)));

        amrex::Real npc_U_0 = xpc_U_0/rpc;
        amrex::Real npc_U_1 = xpc_U_1/rpc;
        amrex::Real npc_U_2 = xpc_U_2/rpc;

        amrex::Real nmc_U_0 = xmc_U_0/rmc;
        amrex::Real nmc_U_1 = xmc_U_1/rmc;
        amrex::Real nmc_U_2 = xmc_U_2/rmc;

        amrex::Real Abarp_UU_00 = -3.0/2.0*Pp_U_0*npc_U_0*(1 - std::pow(npc_U_0, 2))/std::pow(rpc, 2) + 3*Pp_U_0*npc_U_0/std::pow(rpc, 2) - 3.0/2.0*Pp_U_1*npc_U_1*(1 - std::pow(npc_U_0, 2))/std::pow(rpc, 2) - 3.0/2.0*Pp_U_2*npc_U_2*(1 - std::pow(npc_U_0, 2))/std::pow(rpc, 2);
        amrex::Real Abarp_UU_01 = (3.0/2.0)*Pp_U_0*std::pow(npc_U_0, 2)*npc_U_1/std::pow(rpc, 2) + (3.0/2.0)*Pp_U_1*npc_U_0*std::pow(npc_U_1, 2)/std::pow(rpc, 2) + (3.0/2.0)*Pp_U_2*npc_U_0*npc_U_1*npc_U_2/std::pow(rpc, 2) + (3.0/2.0)*(Pp_U_0*npc_U_1 + Pp_U_1*npc_U_0)/std::pow(rpc, 2);
        amrex::Real Abarp_UU_02 = (3.0/2.0)*Pp_U_0*std::pow(npc_U_0, 2)*npc_U_2/std::pow(rpc, 2) + (3.0/2.0)*Pp_U_1*npc_U_0*npc_U_1*npc_U_2/std::pow(rpc, 2) + (3.0/2.0)*Pp_U_2*npc_U_0*std::pow(npc_U_2, 2)/std::pow(rpc, 2) + (3.0/2.0)*(Pp_U_0*npc_U_2 + Pp_U_2*npc_U_0)/std::pow(rpc, 2);
        amrex::Real Abarp_UU_11 = -3.0/2.0*Pp_U_0*npc_U_0*(1 - std::pow(npc_U_1, 2))/std::pow(rpc, 2) - 3.0/2.0*Pp_U_1*npc_U_1*(1 - std::pow(npc_U_1, 2))/std::pow(rpc, 2) + 3*Pp_U_1*npc_U_1/std::pow(rpc, 2) - 3.0/2.0*Pp_U_2*npc_U_2*(1 - std::pow(npc_U_1, 2))/std::pow(rpc, 2);
        amrex::Real Abarp_UU_12 = (3.0/2.0)*Pp_U_0*npc_U_0*npc_U_1*npc_U_2/std::pow(rpc, 2) + (3.0/2.0)*Pp_U_1*std::pow(npc_U_1, 2)*npc_U_2/std::pow(rpc, 2) + (3.0/2.0)*Pp_U_2*npc_U_1*std::pow(npc_U_2, 2)/std::pow(rpc, 2) + (3.0/2.0)*(Pp_U_1*npc_U_2 + Pp_U_2*npc_U_1)/std::pow(rpc, 2);
        amrex::Real Abarp_UU_22 = -3.0/2.0*Pp_U_0*npc_U_0*(1 - std::pow(npc_U_2, 2))/std::pow(rpc, 2) - 3.0/2.0*Pp_U_1*npc_U_1*(1 - std::pow(npc_U_2, 2))/std::pow(rpc, 2) - 3.0/2.0*Pp_U_2*npc_U_2*(1 - std::pow(npc_U_2, 2))/std::pow(rpc, 2) + 3*Pp_U_2*npc_U_2/std::pow(rpc, 2);

        amrex::Real Abarm_UU_00 = -3.0/2.0*Pm_U_0*nmc_U_0*(1 - std::pow(nmc_U_0, 2))/std::pow(rmc, 2) + 3*Pm_U_0*nmc_U_0/std::pow(rmc, 2) - 3.0/2.0*Pm_U_1*nmc_U_1*(1 - std::pow(nmc_U_0, 2))/std::pow(rmc, 2) - 3.0/2.0*Pm_U_2*nmc_U_2*(1 - std::pow(nmc_U_0, 2))/std::pow(rmc, 2);
        amrex::Real Abarm_UU_01 = (3.0/2.0)*Pm_U_0*std::pow(nmc_U_0, 2)*nmc_U_1/std::pow(rmc, 2) + (3.0/2.0)*Pm_U_1*nmc_U_0*std::pow(nmc_U_1, 2)/std::pow(rmc, 2) + (3.0/2.0)*Pm_U_2*nmc_U_0*nmc_U_1*nmc_U_2/std::pow(rmc, 2) + (3.0/2.0)*(Pm_U_0*nmc_U_1 + Pm_U_1*nmc_U_0)/std::pow(rmc, 2);
        amrex::Real Abarm_UU_02 = (3.0/2.0)*Pm_U_0*std::pow(nmc_U_0, 2)*nmc_U_2/std::pow(rmc, 2) + (3.0/2.0)*Pm_U_1*nmc_U_0*nmc_U_1*nmc_U_2/std::pow(rmc, 2) + (3.0/2.0)*Pm_U_2*nmc_U_0*std::pow(nmc_U_2, 2)/std::pow(rmc, 2) + (3.0/2.0)*(Pm_U_0*nmc_U_2 + Pm_U_2*nmc_U_0)/std::pow(rmc, 2);
        amrex::Real Abarm_UU_11 = -3.0/2.0*Pm_U_0*nmc_U_0*(1 - std::pow(nmc_U_1, 2))/std::pow(rmc, 2) - 3.0/2.0*Pm_U_1*nmc_U_1*(1 - std::pow(nmc_U_1, 2))/std::pow(rmc, 2) + 3*Pm_U_1*nmc_U_1/std::pow(rmc, 2) - 3.0/2.0*Pm_U_2*nmc_U_2*(1 - std::pow(nmc_U_1, 2))/std::pow(rmc, 2);
        amrex::Real Abarm_UU_12 = (3.0/2.0)*Pm_U_0*nmc_U_0*nmc_U_1*nmc_U_2/std::pow(rmc, 2) + (3.0/2.0)*Pm_U_1*std::pow(nmc_U_1, 2)*nmc_U_2/std::pow(rmc, 2) + (3.0/2.0)*Pm_U_2*nmc_U_1*std::pow(nmc_U_2, 2)/std::pow(rmc, 2) + (3.0/2.0)*(Pm_U_1*nmc_U_2 + Pm_U_2*nmc_U_1)/std::pow(rmc, 2);
        amrex::Real Abarm_UU_22 = -3.0/2.0*Pm_U_0*nmc_U_0*(1 - std::pow(nmc_U_2, 2))/std::pow(rmc, 2) - 3.0/2.0*Pm_U_1*nmc_U_1*(1 - std::pow(nmc_U_2, 2))/std::pow(rmc, 2) - 3.0/2.0*Pm_U_2*nmc_U_2*(1 - std::pow(nmc_U_2, 2))/std::pow(rmc, 2) + 3*Pm_U_2*nmc_U_2/std::pow(rmc, 2);

        amrex::Real Abar_UU_00 = Abarm_UU_00 + Abarp_UU_00;
        amrex::Real Abar_UU_01 = Abarm_UU_01 + Abarp_UU_01;
        amrex::Real Abar_UU_02 = Abarm_UU_02 + Abarp_UU_02;
        amrex::Real Abar_UU_11 = Abarm_UU_11 + Abarp_UU_11;
        amrex::Real Abar_UU_12 = Abarm_UU_12 + Abarp_UU_12;
        amrex::Real Abar_UU_22 = Abarm_UU_22 + Abarp_UU_22;

        state_fab(i, j, k, Idx::alphafunc) = 1.0/((1.0/2.0)*bareMm/rmc + (1.0/2.0)*bareMp/rpc);

        amrex::Real alphafunc = 1.0/((1.0/2.0)*bareMm/rmc + (1.0/2.0)*bareMp/rpc);

        state_fab(i, j, k, Idx::betafunc) = 0.125*std::pow(Abar_UU_00, 2)*std::pow(alphafunc, 7) + 0.25*std::pow(Abar_UU_01, 2)*std::pow(alphafunc, 7) + 0.25*std::pow(Abar_UU_02, 2)*std::pow(alphafunc, 7) + 0.125*std::pow(Abar_UU_11, 2)*std::pow(alphafunc, 7) + 0.25*std::pow(Abar_UU_12, 2)*std::pow(alphafunc, 7) + 0.125*std::pow(Abar_UU_22, 2)*std::pow(alphafunc, 7);

        state_fab(i, j, k, Idx::psi) = 1 + 1.0/alphafunc;

        amrex::Real psi = 1 + 1.0/alphafunc;

        state_fab(i, j, k, Idx::chi) = std::pow(psi, -4);

        amrex::Real Abar_LL_00 = Abar_UU_00;
        amrex::Real Abar_LL_01 = Abar_UU_01;
        amrex::Real Abar_LL_02 = Abar_UU_02;
        amrex::Real Abar_LL_11 = Abar_UU_11;
        amrex::Real Abar_LL_12 = Abar_UU_12;
        amrex::Real Abar_LL_22 = Abar_UU_22;

        state_fab(i, j, k, Idx::Atilde_LL_00) = Abar_LL_00/std::pow(psi, 6);
        state_fab(i, j, k, Idx::Atilde_LL_01) = Abar_LL_01/std::pow(psi, 6);
        state_fab(i, j, k, Idx::Atilde_LL_02) = Abar_LL_02/std::pow(psi, 6);
        state_fab(i, j, k, Idx::Atilde_LL_11) = Abar_LL_11/std::pow(psi, 6);
        state_fab(i, j, k, Idx::Atilde_LL_12) = Abar_LL_12/std::pow(psi, 6);
        state_fab(i, j, k, Idx::Atilde_LL_22) = Abar_LL_22/std::pow(psi, 6);

        amrex::Real up0 = std::pow(bareMp, 4) + 10*std::pow(bareMp, 3)*rpc + 40*std::pow(bareMp, 2)*std::pow(rpc, 2) + 80*bareMp*std::pow(rpc, 3) + 80*std::pow(rpc, 4);

        amrex::Real um0 = std::pow(bareMm, 4) + 10*std::pow(bareMm, 3)*rmc + 40*std::pow(bareMm, 2)*std::pow(rmc, 2) + 80*bareMm*std::pow(rmc, 3) + 80*std::pow(rmc, 4);

        amrex::Real up2 = (1.0/5.0)*bareMp*(42*std::pow(bareMp, 5)*rpc + 378*std::pow(bareMp, 4)*std::pow(rpc, 2) + 1316*std::pow(bareMp, 3)*std::pow(rpc, 3) + 2156*std::pow(bareMp, 2)*std::pow(rpc, 4) + 1536*bareMp*std::pow(rpc, 5) + 21*bareMp*std::pow(bareMp + 2*rpc, 5)*std::log(bareMp/(bareMp + 2*rpc)) + 240*std::pow(rpc, 6))/std::pow(rpc, 3);

        amrex::Real um2 = (1.0/5.0)*bareMm*(42*std::pow(bareMm, 5)*rmc + 378*std::pow(bareMm, 4)*std::pow(rmc, 2) + 1316*std::pow(bareMm, 3)*std::pow(rmc, 3) + 2156*std::pow(bareMm, 2)*std::pow(rmc, 4) + 1536*bareMm*std::pow(rmc, 5) + 21*bareMm*std::pow(bareMm + 2*rmc, 5)*std::log(bareMm/(bareMm + 2*rmc)) + 240*std::pow(rmc, 6))/std::pow(rmc, 3);

        amrex::Real cosp = xpc_U_1/rpc;

        amrex::Real cosm = xmc_U_1/rmc;

        amrex::Real Lp2 = 1.5*std::pow(cosp, 2) - 0.5;

        amrex::Real Lm2 = 1.5*std::pow(cosm, 2) - 0.5;

        amrex::Real uPp = (1.0/8.0)*std::pow(Pp, 2)*(Lp2*up2 + up0)/(bareMp*std::pow(bareMp + 2*rpc, 5));

        amrex::Real uPm = (1.0/8.0)*std::pow(Pm, 2)*(Lm2*um2 + um0)/(bareMm*std::pow(bareMm + 2*rmc, 5));

        state_fab(i, j, k, Idx::u) = uPm + uPp;
        
}
#endif
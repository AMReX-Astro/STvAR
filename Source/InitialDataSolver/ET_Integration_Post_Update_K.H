#ifndef ET_INTEGRATION_POST_UPDATE_K_H
#define ET_INTEGRATION_POST_UPDATE_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_post_update(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x_U_0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x_U_1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x_U_2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2];
        
        amrex::Real u = state_fab(i, j, k, Idx::u);
        /*
        amrex::Real bareMp = 0.5;

        amrex::Real bareMm = 0.5;

        amrex::Real cp_U_0 = 0; 
        amrex::Real cp_U_1 = 0;
        amrex::Real cp_U_2 = 1.1515;

        amrex::Real cm_U_0 = 0;
        amrex::Real cm_U_1 = 0;
        amrex::Real cm_U_2 = -1.1515;

        amrex::Real Pp_U_0 = 0;
        amrex::Real Pp_U_1 = 0;
        amrex::Real Pp_U_2 = 0;

        amrex::Real Pm_U_0 = 0;
        amrex::Real Pm_U_1 = 0;
        amrex::Real Pm_U_2 = 0;
        */
        
        amrex::Real bareMp = 0.48559999999999998;

        amrex::Real bareMm = 0.48559999999999998;

        amrex::Real cp_U_0 = 0; 
        amrex::Real cp_U_1 = 4.891;
        amrex::Real cp_U_2 = 0;

        amrex::Real cm_U_0 = 0;
        amrex::Real cm_U_1 = -4.891;
        amrex::Real cm_U_2 = 0;

        amrex::Real Pp_U_0 = -0.0969;
        amrex::Real Pp_U_1 = 0;
        amrex::Real Pp_U_2 = 0;

        amrex::Real Pm_U_0 = 0.0969;
        amrex::Real Pm_U_1 = 0;
        amrex::Real Pm_U_2 = 0;
        

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

        state_fab(i, j, k, Idx::psi) = u + 1 + 1.0/alphafunc;

        amrex::Real psi = u + 1 + 1.0/alphafunc;

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
        
}
#endif
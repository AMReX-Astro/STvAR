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
        
        amrex::Real chi = state_fab(i, j, k, Idx::chi);

        amrex::Real alpha = state_fab(i, j, k, Idx::alpha);

        amrex::Real gamtilde_LL_00 = state_fab(i, j, k, Idx::gamtilde_LL_00);
        amrex::Real gamtilde_LL_01 = state_fab(i, j, k, Idx::gamtilde_LL_01);
        amrex::Real gamtilde_LL_02 = state_fab(i, j, k, Idx::gamtilde_LL_02);
        amrex::Real gamtilde_LL_11 = state_fab(i, j, k, Idx::gamtilde_LL_11);
        amrex::Real gamtilde_LL_12 = state_fab(i, j, k, Idx::gamtilde_LL_12);
        amrex::Real gamtilde_LL_22 = state_fab(i, j, k, Idx::gamtilde_LL_22);

        amrex::Real detgamtilde = gamtilde_LL_00*gamtilde_LL_11*gamtilde_LL_22 - gamtilde_LL_00*std::pow(gamtilde_LL_12, 2) - std::pow(gamtilde_LL_01, 2)*gamtilde_LL_22 + 2*gamtilde_LL_01*gamtilde_LL_02*gamtilde_LL_12 - std::pow(gamtilde_LL_02, 2)*gamtilde_LL_11;

        amrex::Real scalefactor = std::pow(detgamtilde, -0.33333333333333331);

        state_fab(i, j, k, Idx::gamtilde_LL_00) = gamtilde_LL_00*scalefactor;
        state_fab(i, j, k, Idx::gamtilde_LL_01) = gamtilde_LL_01*scalefactor;
        state_fab(i, j, k, Idx::gamtilde_LL_02) = gamtilde_LL_02*scalefactor;
        state_fab(i, j, k, Idx::gamtilde_LL_11) = gamtilde_LL_11*scalefactor;
        state_fab(i, j, k, Idx::gamtilde_LL_12) = gamtilde_LL_12*scalefactor;
        state_fab(i, j, k, Idx::gamtilde_LL_22) = gamtilde_LL_22*scalefactor;

        amrex::Real Atilde_LL_00 = state_fab(i, j, k, Idx::Atilde_LL_00);
        amrex::Real Atilde_LL_01 = state_fab(i, j, k, Idx::Atilde_LL_01);
        amrex::Real Atilde_LL_02 = state_fab(i, j, k, Idx::Atilde_LL_02);
        amrex::Real Atilde_LL_11 = state_fab(i, j, k, Idx::Atilde_LL_11);
        amrex::Real Atilde_LL_12 = state_fab(i, j, k, Idx::Atilde_LL_12);
        amrex::Real Atilde_LL_22 = state_fab(i, j, k, Idx::Atilde_LL_22);

        amrex::Real invgamtilde_UU_00 = (gamtilde_LL_11*gamtilde_LL_22 - std::pow(gamtilde_LL_12, 2))/(gamtilde_LL_00*gamtilde_LL_11*gamtilde_LL_22 - gamtilde_LL_00*std::pow(gamtilde_LL_12, 2) - std::pow(gamtilde_LL_01, 2)*gamtilde_LL_22 + 2*gamtilde_LL_01*gamtilde_LL_02*gamtilde_LL_12 - std::pow(gamtilde_LL_02, 2)*gamtilde_LL_11);
        amrex::Real invgamtilde_UU_01 = (-gamtilde_LL_01*gamtilde_LL_22 + gamtilde_LL_02*gamtilde_LL_12)/(gamtilde_LL_00*gamtilde_LL_11*gamtilde_LL_22 - gamtilde_LL_00*std::pow(gamtilde_LL_12, 2) - std::pow(gamtilde_LL_01, 2)*gamtilde_LL_22 + 2*gamtilde_LL_01*gamtilde_LL_02*gamtilde_LL_12 - std::pow(gamtilde_LL_02, 2)*gamtilde_LL_11);
        amrex::Real invgamtilde_UU_02 = (gamtilde_LL_01*gamtilde_LL_12 - gamtilde_LL_02*gamtilde_LL_11)/(gamtilde_LL_00*gamtilde_LL_11*gamtilde_LL_22 - gamtilde_LL_00*std::pow(gamtilde_LL_12, 2) - std::pow(gamtilde_LL_01, 2)*gamtilde_LL_22 + 2*gamtilde_LL_01*gamtilde_LL_02*gamtilde_LL_12 - std::pow(gamtilde_LL_02, 2)*gamtilde_LL_11);
        amrex::Real invgamtilde_UU_11 = (gamtilde_LL_00*gamtilde_LL_22 - std::pow(gamtilde_LL_02, 2))/(gamtilde_LL_00*gamtilde_LL_11*gamtilde_LL_22 - gamtilde_LL_00*std::pow(gamtilde_LL_12, 2) - std::pow(gamtilde_LL_01, 2)*gamtilde_LL_22 + 2*gamtilde_LL_01*gamtilde_LL_02*gamtilde_LL_12 - std::pow(gamtilde_LL_02, 2)*gamtilde_LL_11);
        amrex::Real invgamtilde_UU_12 = (-gamtilde_LL_00*gamtilde_LL_12 + gamtilde_LL_01*gamtilde_LL_02)/(gamtilde_LL_00*gamtilde_LL_11*gamtilde_LL_22 - gamtilde_LL_00*std::pow(gamtilde_LL_12, 2) - std::pow(gamtilde_LL_01, 2)*gamtilde_LL_22 + 2*gamtilde_LL_01*gamtilde_LL_02*gamtilde_LL_12 - std::pow(gamtilde_LL_02, 2)*gamtilde_LL_11);
        amrex::Real invgamtilde_UU_22 = (gamtilde_LL_00*gamtilde_LL_11 - std::pow(gamtilde_LL_01, 2))/(gamtilde_LL_00*gamtilde_LL_11*gamtilde_LL_22 - gamtilde_LL_00*std::pow(gamtilde_LL_12, 2) - std::pow(gamtilde_LL_01, 2)*gamtilde_LL_22 + 2*gamtilde_LL_01*gamtilde_LL_02*gamtilde_LL_12 - std::pow(gamtilde_LL_02, 2)*gamtilde_LL_11);

        amrex::Real TrAtilde = Atilde_LL_00*invgamtilde_UU_00 + 2*Atilde_LL_01*invgamtilde_UU_01 + 2*Atilde_LL_02*invgamtilde_UU_02 + Atilde_LL_11*invgamtilde_UU_11 + 2*Atilde_LL_12*invgamtilde_UU_12 + Atilde_LL_22*invgamtilde_UU_22;

        state_fab(i, j, k, Idx::Atilde_LL_00) = Atilde_LL_00 - 0.33333333333333331*TrAtilde*gamtilde_LL_00;
        state_fab(i, j, k, Idx::Atilde_LL_01) = Atilde_LL_01 - 0.33333333333333331*TrAtilde*gamtilde_LL_01;
        state_fab(i, j, k, Idx::Atilde_LL_02) = Atilde_LL_02 - 0.33333333333333331*TrAtilde*gamtilde_LL_02;
        state_fab(i, j, k, Idx::Atilde_LL_11) = Atilde_LL_11 - 0.33333333333333331*TrAtilde*gamtilde_LL_11;
        state_fab(i, j, k, Idx::Atilde_LL_12) = Atilde_LL_12 - 0.33333333333333331*TrAtilde*gamtilde_LL_12;
        state_fab(i, j, k, Idx::Atilde_LL_22) = Atilde_LL_22 - 0.33333333333333331*TrAtilde*gamtilde_LL_22;

        state_fab(i, j, k, Idx::chi) = std::max(0.00100000000000000, chi);

        state_fab(i, j, k, Idx::alpha) = std::max(0.00100000000000000, alpha);
        
}
#endif
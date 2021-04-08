#ifndef ET_INTEGRATION_LOADINITDATA_K_H
#define ET_INTEGRATION_LOADINITDATA_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
initialize_from_data(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab, 
        amrex::Array4<amrex::Real const> const& initial_data,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x_U_0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x_U_1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x_U_2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        state_fab(i, j, k, Idx::chi) = initial_data(i, j, k, InitIdx::chi);

        amrex::Real chi = initial_data(i, j, k, InitIdx::chi);

        state_fab(i, j, k, Idx::gamtilde_LL_00) = 1;
        state_fab(i, j, k, Idx::gamtilde_LL_01) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_02) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_11) = 1;
        state_fab(i, j, k, Idx::gamtilde_LL_12) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_22) = 1;

        state_fab(i, j, k, Idx::Khat) = 0;

        state_fab(i, j, k, Idx::Atilde_LL_00) = initial_data(i, j, k, InitIdx::Atilde_LL_00);
        state_fab(i, j, k, Idx::Atilde_LL_01) = initial_data(i, j, k, InitIdx::Atilde_LL_01);
        state_fab(i, j, k, Idx::Atilde_LL_02) = initial_data(i, j, k, InitIdx::Atilde_LL_02);
        state_fab(i, j, k, Idx::Atilde_LL_11) = initial_data(i, j, k, InitIdx::Atilde_LL_11);
        state_fab(i, j, k, Idx::Atilde_LL_12) = initial_data(i, j, k, InitIdx::Atilde_LL_12);
        state_fab(i, j, k, Idx::Atilde_LL_22) = initial_data(i, j, k, InitIdx::Atilde_LL_22);

        state_fab(i, j, k, Idx::theta) = 0;

        state_fab(i, j, k, Idx::Gamtilde_U_0) = 0;
        state_fab(i, j, k, Idx::Gamtilde_U_1) = 0;
        state_fab(i, j, k, Idx::Gamtilde_U_2) = 0;

        state_fab(i, j, k, Idx::alpha) = std::sqrt(chi);

        state_fab(i, j, k, Idx::beta_U_0) = 0;
        state_fab(i, j, k, Idx::beta_U_1) = 0;
        state_fab(i, j, k, Idx::beta_U_2) = 0;

}
#endif
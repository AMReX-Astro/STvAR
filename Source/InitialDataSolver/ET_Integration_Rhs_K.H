#ifndef ET_INTEGRATION_RHS_K_H
#define ET_INTEGRATION_RHS_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_rhs(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& rhs_fab, 
        amrex::Array4<amrex::Real const> const& state_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        amrex::Real u = state_fab(i, j, k, Idx::u);

        amrex::Real dDu_L_0 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::u) - 1.0/12.0*state_fab(i + 2, j, k, Idx::u) - 2.0/3.0*state_fab(i - 1, j, k, Idx::u) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::u))/dx[0];
        amrex::Real dDu_L_1 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::u) - 1.0/12.0*state_fab(i, j + 2, k, Idx::u) - 2.0/3.0*state_fab(i, j - 1, k, Idx::u) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::u))/dx[1];
        amrex::Real dDu_L_2 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::u) - 1.0/12.0*state_fab(i, j, k + 2, Idx::u) - 2.0/3.0*state_fab(i, j, k - 1, Idx::u) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::u))/dx[2];

        amrex::Real dDDu_LL_00 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::u) - 1.0/12.0*state_fab(i + 2, j, k, Idx::u) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::u) - 1.0/12.0*state_fab(i - 2, j, k, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[0], 2);
        amrex::Real dDDu_LL_01 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::u) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::u) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::u) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::u) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::u) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::u) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::u) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::u) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::u) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::u) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::u) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::u) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::u) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::u) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::u) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::u))/(dx[0]*dx[1]);
        amrex::Real dDDu_LL_02 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::u) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::u) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::u) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::u) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::u) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::u) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::u) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::u) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::u) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::u) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::u) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::u) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::u) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::u) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::u) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::u))/(dx[0]*dx[2]);
        amrex::Real dDDu_LL_11 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::u) - 1.0/12.0*state_fab(i, j + 2, k, Idx::u) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::u) - 1.0/12.0*state_fab(i, j - 2, k, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[1], 2);
        amrex::Real dDDu_LL_12 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::u) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::u) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::u) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::u) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::u) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::u) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::u) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::u) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::u) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::u) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::u) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::u) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::u) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::u) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::u) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::u))/(dx[1]*dx[2]);
        amrex::Real dDDu_LL_22 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::u) - 1.0/12.0*state_fab(i, j, k + 2, Idx::u) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::u) - 1.0/12.0*state_fab(i, j, k - 2, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[2], 2);

        amrex::Real alphafunc = state_fab(i, j, k, Idx::alphafunc);

        amrex::Real betafunc = state_fab(i, j, k, Idx::betafunc);

        rhs_fab(i, j, k, Idx::u) = betafunc/std::pow(alphafunc*(u + 1) + 1, 7) + dDDu_LL_00 + dDDu_LL_11 + dDDu_LL_22;
        

}
#endif
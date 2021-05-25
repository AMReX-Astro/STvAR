#ifndef ET_INTEGRATION_DIAG_K_H
#define ET_INTEGRATION_DIAG_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_diagnostics(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& diag, 
        amrex::Array4<amrex::Real const> const& state_fab,
        const amrex::Real time_lev,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        amrex::Real v = 5.0;
        
        amrex::Real PhiReal = state_fab(i, j, k, Idx::PhiReal);
        amrex::Real PhiImaginary = state_fab(i, j, k, Idx::PhiImaginary);
        
        amrex::Real Chi = state_fab(i, j, k, Idx::Chi);

        amrex::Real Chi0 = state_fab(i, j, k, Idx::Chi0);

        amrex::Real Chi1 = state_fab(i, j, k, Idx::Chi1);
        
        amrex::Real PiAT = state_fab(i, j, k, Idx::PiAT);
        
        amrex::Real dDA_LL_00 = ((1.0/2.0)*state_fab(i + 1, j, k, Idx::A_L_0) - 1.0/2.0*state_fab(i - 1, j, k, Idx::A_L_0))/dx[0];
        amrex::Real dDA_LL_01 = ((1.0/2.0)*state_fab(i, j + 1, k, Idx::A_L_0) - 1.0/2.0*state_fab(i, j - 1, k, Idx::A_L_0))/dx[1];
        amrex::Real dDA_LL_02 = ((1.0/2.0)*state_fab(i, j, k + 1, Idx::A_L_0) - 1.0/2.0*state_fab(i, j, k - 1, Idx::A_L_0))/dx[2];
        amrex::Real dDA_LL_10 = ((1.0/2.0)*state_fab(i + 1, j, k, Idx::A_L_1) - 1.0/2.0*state_fab(i - 1, j, k, Idx::A_L_1))/dx[0];
        amrex::Real dDA_LL_11 = ((1.0/2.0)*state_fab(i, j + 1, k, Idx::A_L_1) - 1.0/2.0*state_fab(i, j - 1, k, Idx::A_L_1))/dx[1];
        amrex::Real dDA_LL_12 = ((1.0/2.0)*state_fab(i, j, k + 1, Idx::A_L_1) - 1.0/2.0*state_fab(i, j, k - 1, Idx::A_L_1))/dx[2];
        amrex::Real dDA_LL_20 = ((1.0/2.0)*state_fab(i + 1, j, k, Idx::A_L_2) - 1.0/2.0*state_fab(i - 1, j, k, Idx::A_L_2))/dx[0];
        amrex::Real dDA_LL_21 = ((1.0/2.0)*state_fab(i, j + 1, k, Idx::A_L_2) - 1.0/2.0*state_fab(i, j - 1, k, Idx::A_L_2))/dx[1];
        amrex::Real dDA_LL_22 = ((1.0/2.0)*state_fab(i, j, k + 1, Idx::A_L_2) - 1.0/2.0*state_fab(i, j, k - 1, Idx::A_L_2))/dx[2];
        
        diag(i, j, k, Diag::ModPhi) = v-std::sqrt(std::pow(PhiReal,2) + std::pow(PhiImaginary,2));
        
        diag(i, j, k, Diag::divA) = -PiAT + dDA_LL_00 + dDA_LL_11 + dDA_LL_22;
        
        diag(i, j, k, Diag::ModChi) = std::sqrt(std::pow(Chi,2)+std::pow(Chi0,2)+std::pow(Chi1,2));

}
#endif
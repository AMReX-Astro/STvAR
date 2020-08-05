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
        amrex::Real Psi = state_fab(i, j, k, Idx::Psi);

        amrex::Real Pi = state_fab(i, j, k, Idx::Pi);

        amrex::Real dPsi0 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::Psi) - 1.0/12.0*state_fab(i + 2, j, k, Idx::Psi) - 2.0/3.0*state_fab(i - 1, j, k, Idx::Psi) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::Psi))/dx[0];
        amrex::Real dPsi1 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::Psi) - 1.0/12.0*state_fab(i, j + 2, k, Idx::Psi) - 2.0/3.0*state_fab(i, j - 1, k, Idx::Psi) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::Psi))/dx[1];
        amrex::Real dPsi2 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::Psi) - 1.0/12.0*state_fab(i, j, k + 2, Idx::Psi) - 2.0/3.0*state_fab(i, j, k - 1, Idx::Psi) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::Psi))/dx[2];
        amrex::Real m = 1;

        diag(i, j, k, Diag::EnDens) = 0.5*std::pow(Pi, 2) + 0.5*std::pow(Psi, 2)*std::pow(m, 2) + 0.5*std::pow(dPsi0, 2) + 0.5*std::pow(dPsi1, 2);

}
#endif
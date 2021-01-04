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
        amrex::Array4<amrex::Real> const& stress_energy_arr,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real Psi = state_fab(i, j, k, Idx::Psi);

        amrex::Real Pi = state_fab(i, j, k, Idx::Pi);

        amrex::Real dPsi0 = ((1.0/2.0)*state_fab(i + 1, j, k, Idx::Psi) - 1.0/2.0*state_fab(i - 1, j, k, Idx::Psi))/dx[0];
        amrex::Real dPsi1 = ((1.0/2.0)*state_fab(i, j + 1, k, Idx::Psi) - 1.0/2.0*state_fab(i, j - 1, k, Idx::Psi))/dx[1];
        amrex::Real dPsi2 = ((1.0/2.0)*state_fab(i, j, k + 1, Idx::Psi) - 1.0/2.0*state_fab(i, j, k - 1, Idx::Psi))/dx[2];
        amrex::Real ddPsi00 = (state_fab(i + 1, j, k, Idx::Psi) + state_fab(i - 1, j, k, Idx::Psi) - 2*state_fab(i, j, k, Idx::Psi))/std::pow(dx[0], 2);
        amrex::Real ddPsi01 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::Psi) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::Psi) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::Psi) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::Psi))/(dx[0]*dx[1]);
        amrex::Real ddPsi02 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::Psi) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::Psi) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::Psi) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::Psi))/(dx[0]*dx[2]);
        amrex::Real ddPsi10 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::Psi) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::Psi) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::Psi) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::Psi))/(dx[0]*dx[1]);
        amrex::Real ddPsi11 = (state_fab(i, j + 1, k, Idx::Psi) + state_fab(i, j - 1, k, Idx::Psi) - 2*state_fab(i, j, k, Idx::Psi))/std::pow(dx[1], 2);
        amrex::Real ddPsi12 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::Psi) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::Psi) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::Psi) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::Psi))/(dx[1]*dx[2]);
        amrex::Real ddPsi20 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::Psi) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::Psi) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::Psi) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::Psi))/(dx[0]*dx[2]);
        amrex::Real ddPsi21 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::Psi) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::Psi) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::Psi) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::Psi))/(dx[1]*dx[2]);
        amrex::Real ddPsi22 = (state_fab(i, j, k + 1, Idx::Psi) + state_fab(i, j, k - 1, Idx::Psi) - 2*state_fab(i, j, k, Idx::Psi))/std::pow(dx[2], 2);
        amrex::Real m = 1;

        rhs_fab(i, j, k, Idx::Psi) = Pi;

        rhs_fab(i, j, k, Idx::Pi) = -Psi*std::pow(m, 2) + ddPsi00 + ddPsi11;


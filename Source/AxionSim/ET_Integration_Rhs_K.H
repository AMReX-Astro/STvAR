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
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        amrex::Real Psi1 = state_fab(i, j, k, Idx::Psi1);

        amrex::Real Pi1 = state_fab(i, j, k, Idx::Pi1);

        amrex::Real Psi2 = state_fab(i, j, k, Idx::Psi2);

        amrex::Real Pi2 = state_fab(i, j, k, Idx::Pi2);

        amrex::Real ddPsi100 = (state_fab(i + 1, j, k, Idx::Psi1) + state_fab(i - 1, j, k, Idx::Psi1) - 2*state_fab(i, j, k, Idx::Psi1))/std::pow(dx[0], 2);
        amrex::Real ddPsi101 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::Psi1) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::Psi1) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::Psi1) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::Psi1))/(dx[0]*dx[1]);
        amrex::Real ddPsi102 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::Psi1) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::Psi1) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::Psi1) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::Psi1))/(dx[0]*dx[2]);
        amrex::Real ddPsi110 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::Psi1) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::Psi1) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::Psi1) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::Psi1))/(dx[0]*dx[1]);
        amrex::Real ddPsi111 = (state_fab(i, j + 1, k, Idx::Psi1) + state_fab(i, j - 1, k, Idx::Psi1) - 2*state_fab(i, j, k, Idx::Psi1))/std::pow(dx[1], 2);
        amrex::Real ddPsi112 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::Psi1) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::Psi1) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::Psi1) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::Psi1))/(dx[1]*dx[2]);
        amrex::Real ddPsi120 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::Psi1) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::Psi1) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::Psi1) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::Psi1))/(dx[0]*dx[2]);
        amrex::Real ddPsi121 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::Psi1) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::Psi1) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::Psi1) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::Psi1))/(dx[1]*dx[2]);
        amrex::Real ddPsi122 = (state_fab(i, j, k + 1, Idx::Psi1) + state_fab(i, j, k - 1, Idx::Psi1) - 2*state_fab(i, j, k, Idx::Psi1))/std::pow(dx[2], 2);
        amrex::Real ddPsi200 = (state_fab(i + 1, j, k, Idx::Psi2) + state_fab(i - 1, j, k, Idx::Psi2) - 2*state_fab(i, j, k, Idx::Psi2))/std::pow(dx[0], 2);
        amrex::Real ddPsi201 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::Psi2) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::Psi2) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::Psi2) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::Psi2))/(dx[0]*dx[1]);
        amrex::Real ddPsi202 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::Psi2) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::Psi2) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::Psi2) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::Psi2))/(dx[0]*dx[2]);
        amrex::Real ddPsi210 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::Psi2) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::Psi2) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::Psi2) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::Psi2))/(dx[0]*dx[1]);
        amrex::Real ddPsi211 = (state_fab(i, j + 1, k, Idx::Psi2) + state_fab(i, j - 1, k, Idx::Psi2) - 2*state_fab(i, j, k, Idx::Psi2))/std::pow(dx[1], 2);
        amrex::Real ddPsi212 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::Psi2) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::Psi2) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::Psi2) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::Psi2))/(dx[1]*dx[2]);
        amrex::Real ddPsi220 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::Psi2) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::Psi2) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::Psi2) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::Psi2))/(dx[0]*dx[2]);
        amrex::Real ddPsi221 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::Psi2) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::Psi2) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::Psi2) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::Psi2))/(dx[1]*dx[2]);
        amrex::Real ddPsi222 = (state_fab(i, j, k + 1, Idx::Psi2) + state_fab(i, j, k - 1, Idx::Psi2) - 2*state_fab(i, j, k, Idx::Psi2))/std::pow(dx[2], 2);
        amrex::Real lambda = 1.0;

        amrex::Real eta = 1.0;

        rhs_fab(i, j, k, Idx::Psi1) = Pi1;

        rhs_fab(i, j, k, Idx::Pi1) = -2*Pi1/eta + ddPsi100 + ddPsi111 + ddPsi122;

        rhs_fab(i, j, k, Idx::Psi2) = Pi2;

        rhs_fab(i, j, k, Idx::Pi2) = -2*Pi2/eta + ddPsi200 + ddPsi211 + ddPsi222;

}
#endif
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
        
        
        amrex::Real u = state_fab(i, j, k, Idx::u);

        amrex::Real psi = state_fab(i, j, k, Idx::psi);

        amrex::Real chi = state_fab(i, j, k, Idx::chi);

        amrex::Real alphafunc = state_fab(i, j, k, Idx::alphafunc);

        amrex::Real betafunc = state_fab(i, j, k, Idx::betafunc);

        amrex::Real dDu_L_0 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::u) - 1.0/12.0*state_fab(i + 2, j, k, Idx::u) - 2.0/3.0*state_fab(i - 1, j, k, Idx::u) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::u))/dx[0];
        amrex::Real dDu_L_1 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::u) - 1.0/12.0*state_fab(i, j + 2, k, Idx::u) - 2.0/3.0*state_fab(i, j - 1, k, Idx::u) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::u))/dx[1];
        amrex::Real dDu_L_2 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::u) - 1.0/12.0*state_fab(i, j, k + 2, Idx::u) - 2.0/3.0*state_fab(i, j, k - 1, Idx::u) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::u))/dx[2];

        amrex::Real dDDu_LL_00 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::u) - 1.0/12.0*state_fab(i + 2, j, k, Idx::u) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::u) - 1.0/12.0*state_fab(i - 2, j, k, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[0], 2);
        amrex::Real dDDu_LL_01 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::u) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::u) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::u) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::u) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::u) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::u) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::u) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::u) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::u) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::u) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::u) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::u) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::u) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::u) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::u) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::u))/(dx[0]*dx[1]);
        amrex::Real dDDu_LL_02 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::u) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::u) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::u) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::u) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::u) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::u) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::u) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::u) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::u) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::u) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::u) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::u) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::u) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::u) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::u) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::u))/(dx[0]*dx[2]);
        amrex::Real dDDu_LL_11 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::u) - 1.0/12.0*state_fab(i, j + 2, k, Idx::u) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::u) - 1.0/12.0*state_fab(i, j - 2, k, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[1], 2);
        amrex::Real dDDu_LL_12 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::u) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::u) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::u) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::u) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::u) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::u) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::u) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::u) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::u) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::u) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::u) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::u) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::u) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::u) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::u) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::u))/(dx[1]*dx[2]);
        amrex::Real dDDu_LL_22 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::u) - 1.0/12.0*state_fab(i, j, k + 2, Idx::u) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::u) - 1.0/12.0*state_fab(i, j, k - 2, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[2], 2);

        amrex::Real dDpsi_L_0 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::psi) - 1.0/12.0*state_fab(i + 2, j, k, Idx::psi) - 2.0/3.0*state_fab(i - 1, j, k, Idx::psi) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::psi))/dx[0];
        amrex::Real dDpsi_L_1 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::psi) - 1.0/12.0*state_fab(i, j + 2, k, Idx::psi) - 2.0/3.0*state_fab(i, j - 1, k, Idx::psi) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::psi))/dx[1];
        amrex::Real dDpsi_L_2 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::psi) - 1.0/12.0*state_fab(i, j, k + 2, Idx::psi) - 2.0/3.0*state_fab(i, j, k - 1, Idx::psi) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::psi))/dx[2];

        amrex::Real dDDpsi_LL_00 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::psi) - 1.0/12.0*state_fab(i + 2, j, k, Idx::psi) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::psi) - 1.0/12.0*state_fab(i - 2, j, k, Idx::psi) - 5.0/2.0*state_fab(i, j, k, Idx::psi))/std::pow(dx[0], 2);
        amrex::Real dDDpsi_LL_01 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::psi) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::psi) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::psi) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::psi) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::psi) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::psi) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::psi) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::psi) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::psi) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::psi) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::psi) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::psi) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::psi) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::psi) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::psi) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::psi))/(dx[0]*dx[1]);
        amrex::Real dDDpsi_LL_02 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::psi) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::psi) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::psi) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::psi) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::psi) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::psi) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::psi) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::psi) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::psi) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::psi) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::psi) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::psi) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::psi) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::psi) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::psi) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::psi))/(dx[0]*dx[2]);
        amrex::Real dDDpsi_LL_11 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::psi) - 1.0/12.0*state_fab(i, j + 2, k, Idx::psi) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::psi) - 1.0/12.0*state_fab(i, j - 2, k, Idx::psi) - 5.0/2.0*state_fab(i, j, k, Idx::psi))/std::pow(dx[1], 2);
        amrex::Real dDDpsi_LL_12 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::psi) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::psi) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::psi) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::psi) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::psi) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::psi) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::psi) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::psi) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::psi) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::psi) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::psi) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::psi) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::psi) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::psi) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::psi) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::psi))/(dx[1]*dx[2]);
        amrex::Real dDDpsi_LL_22 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::psi) - 1.0/12.0*state_fab(i, j, k + 2, Idx::psi) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::psi) - 1.0/12.0*state_fab(i, j, k - 2, Idx::psi) - 5.0/2.0*state_fab(i, j, k, Idx::psi))/std::pow(dx[2], 2);

        amrex::Real dDchi_L_0 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::chi) - 1.0/12.0*state_fab(i + 2, j, k, Idx::chi) - 2.0/3.0*state_fab(i - 1, j, k, Idx::chi) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::chi))/dx[0];
        amrex::Real dDchi_L_1 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::chi) - 1.0/12.0*state_fab(i, j + 2, k, Idx::chi) - 2.0/3.0*state_fab(i, j - 1, k, Idx::chi) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::chi))/dx[1];
        amrex::Real dDchi_L_2 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::chi) - 1.0/12.0*state_fab(i, j, k + 2, Idx::chi) - 2.0/3.0*state_fab(i, j, k - 1, Idx::chi) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::chi))/dx[2];

        amrex::Real dDDchi_LL_00 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::chi) - 1.0/12.0*state_fab(i + 2, j, k, Idx::chi) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::chi) - 1.0/12.0*state_fab(i - 2, j, k, Idx::chi) - 5.0/2.0*state_fab(i, j, k, Idx::chi))/std::pow(dx[0], 2);
        amrex::Real dDDchi_LL_01 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::chi) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::chi) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::chi) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::chi) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::chi) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::chi) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::chi) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::chi) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::chi) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::chi) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::chi) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::chi) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::chi) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::chi) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::chi) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::chi))/(dx[0]*dx[1]);
        amrex::Real dDDchi_LL_02 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::chi) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::chi) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::chi) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::chi) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::chi) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::chi) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::chi) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::chi) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::chi) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::chi) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::chi) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::chi) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::chi) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::chi) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::chi) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::chi))/(dx[0]*dx[2]);
        amrex::Real dDDchi_LL_11 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::chi) - 1.0/12.0*state_fab(i, j + 2, k, Idx::chi) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::chi) - 1.0/12.0*state_fab(i, j - 2, k, Idx::chi) - 5.0/2.0*state_fab(i, j, k, Idx::chi))/std::pow(dx[1], 2);
        amrex::Real dDDchi_LL_12 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::chi) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::chi) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::chi) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::chi) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::chi) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::chi) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::chi) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::chi) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::chi) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::chi) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::chi) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::chi) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::chi) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::chi) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::chi) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::chi))/(dx[1]*dx[2]);
        amrex::Real dDDchi_LL_22 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::chi) - 1.0/12.0*state_fab(i, j, k + 2, Idx::chi) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::chi) - 1.0/12.0*state_fab(i, j, k - 2, Idx::chi) - 5.0/2.0*state_fab(i, j, k, Idx::chi))/std::pow(dx[2], 2);

        amrex::Real RsclrPsi = -8*dDDpsi_LL_00/std::pow(psi, 5) - 8*dDDpsi_LL_11/std::pow(psi, 5) - 8*dDDpsi_LL_22/std::pow(psi, 5);

        amrex::Real RsclrChi = 2*dDDchi_LL_00 + 2*dDDchi_LL_11 + 2*dDDchi_LL_22 - 5.0/2.0*std::pow(dDchi_L_0, 2)/chi - 5.0/2.0*std::pow(dDchi_L_1, 2)/chi - 5.0/2.0*std::pow(dDchi_L_2, 2)/chi;

        amrex::Real Atilde_LL_00 = state_fab(i, j, k, Idx::Atilde_LL_00);
        amrex::Real Atilde_LL_01 = state_fab(i, j, k, Idx::Atilde_LL_01);
        amrex::Real Atilde_LL_02 = state_fab(i, j, k, Idx::Atilde_LL_02);
        amrex::Real Atilde_LL_11 = state_fab(i, j, k, Idx::Atilde_LL_11);
        amrex::Real Atilde_LL_12 = state_fab(i, j, k, Idx::Atilde_LL_12);
        amrex::Real Atilde_LL_22 = state_fab(i, j, k, Idx::Atilde_LL_22);

        amrex::Real Atilde_UU_00 = Atilde_LL_00;
        amrex::Real Atilde_UU_01 = Atilde_LL_01;
        amrex::Real Atilde_UU_02 = Atilde_LL_02;
        amrex::Real Atilde_UU_11 = Atilde_LL_11;
        amrex::Real Atilde_UU_12 = Atilde_LL_12;
        amrex::Real Atilde_UU_22 = Atilde_LL_22;

        diag(i, j, k, Diag::Hpsi) = Atilde_LL_00*Atilde_UU_00 + 2*Atilde_LL_01*Atilde_UU_01 + 2*Atilde_LL_02*Atilde_UU_02 + Atilde_LL_11*Atilde_UU_11 + 2*Atilde_LL_12*Atilde_UU_12 + Atilde_LL_22*Atilde_UU_22 + RsclrPsi;

        diag(i, j, k, Diag::Hchi) = Atilde_LL_00*Atilde_UU_00 + 2*Atilde_LL_01*Atilde_UU_01 + 2*Atilde_LL_02*Atilde_UU_02 + Atilde_LL_11*Atilde_UU_11 + 2*Atilde_LL_12*Atilde_UU_12 + Atilde_LL_22*Atilde_UU_22 + RsclrChi;

        diag(i, j, k, Diag::Hu) = -8*(betafunc/std::pow(alphafunc*(u + 1) + 1, 7) + dDDu_LL_00 + dDDu_LL_11 + dDDu_LL_22)/std::pow(psi, 5);
        
}

AMREX_GPU_DEVICE 
inline 
amrex::Real
sum_C_constraint(const int i, const int j, const int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom) noexcept 
{
    
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x_U_0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x_U_1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x_U_2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        return 0;
        
}
#endif
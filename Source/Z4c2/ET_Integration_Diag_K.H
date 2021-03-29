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
        
        amrex::Real chi = state_fab(i, j, k, Idx::chi);

        amrex::Real gamtilde_LL_00 = state_fab(i, j, k, Idx::gamtilde_LL_00);
        amrex::Real gamtilde_LL_01 = state_fab(i, j, k, Idx::gamtilde_LL_01);
        amrex::Real gamtilde_LL_02 = state_fab(i, j, k, Idx::gamtilde_LL_02);
        amrex::Real gamtilde_LL_11 = state_fab(i, j, k, Idx::gamtilde_LL_11);
        amrex::Real gamtilde_LL_12 = state_fab(i, j, k, Idx::gamtilde_LL_12);
        amrex::Real gamtilde_LL_22 = state_fab(i, j, k, Idx::gamtilde_LL_22);

        amrex::Real Khat = state_fab(i, j, k, Idx::Khat);

        amrex::Real Atilde_LL_00 = state_fab(i, j, k, Idx::Atilde_LL_00);
        amrex::Real Atilde_LL_01 = state_fab(i, j, k, Idx::Atilde_LL_01);
        amrex::Real Atilde_LL_02 = state_fab(i, j, k, Idx::Atilde_LL_02);
        amrex::Real Atilde_LL_11 = state_fab(i, j, k, Idx::Atilde_LL_11);
        amrex::Real Atilde_LL_12 = state_fab(i, j, k, Idx::Atilde_LL_12);
        amrex::Real Atilde_LL_22 = state_fab(i, j, k, Idx::Atilde_LL_22);

        amrex::Real theta = state_fab(i, j, k, Idx::theta);

        amrex::Real Gamtilde_U_0 = state_fab(i, j, k, Idx::Gamtilde_U_0);
        amrex::Real Gamtilde_U_1 = state_fab(i, j, k, Idx::Gamtilde_U_1);
        amrex::Real Gamtilde_U_2 = state_fab(i, j, k, Idx::Gamtilde_U_2);

        amrex::Real alpha = state_fab(i, j, k, Idx::alpha);

        amrex::Real beta_U_0 = state_fab(i, j, k, Idx::beta_U_0);
        amrex::Real beta_U_1 = state_fab(i, j, k, Idx::beta_U_1);
        amrex::Real beta_U_2 = state_fab(i, j, k, Idx::beta_U_2);

        amrex::Real dDchi_L_0 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::chi) - 1.0/12.0*state_fab(i + 2, j, k, Idx::chi) - 2.0/3.0*state_fab(i - 1, j, k, Idx::chi) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::chi))/dx[0];
        amrex::Real dDchi_L_1 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::chi) - 1.0/12.0*state_fab(i, j + 2, k, Idx::chi) - 2.0/3.0*state_fab(i, j - 1, k, Idx::chi) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::chi))/dx[1];
        amrex::Real dDchi_L_2 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::chi) - 1.0/12.0*state_fab(i, j, k + 2, Idx::chi) - 2.0/3.0*state_fab(i, j, k - 1, Idx::chi) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::chi))/dx[2];

        amrex::Real dDDchi_LL_00 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::chi) - 1.0/12.0*state_fab(i + 2, j, k, Idx::chi) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::chi) - 1.0/12.0*state_fab(i - 2, j, k, Idx::chi) - 5.0/2.0*state_fab(i, j, k, Idx::chi))/std::pow(dx[0], 2);
        amrex::Real dDDchi_LL_01 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::chi) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::chi) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::chi) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::chi) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::chi) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::chi) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::chi) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::chi) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::chi) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::chi) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::chi) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::chi) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::chi) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::chi) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::chi) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::chi))/(dx[0]*dx[1]);
        amrex::Real dDDchi_LL_02 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::chi) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::chi) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::chi) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::chi) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::chi) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::chi) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::chi) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::chi) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::chi) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::chi) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::chi) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::chi) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::chi) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::chi) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::chi) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::chi))/(dx[0]*dx[2]);
        amrex::Real dDDchi_LL_11 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::chi) - 1.0/12.0*state_fab(i, j + 2, k, Idx::chi) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::chi) - 1.0/12.0*state_fab(i, j - 2, k, Idx::chi) - 5.0/2.0*state_fab(i, j, k, Idx::chi))/std::pow(dx[1], 2);
        amrex::Real dDDchi_LL_12 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::chi) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::chi) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::chi) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::chi) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::chi) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::chi) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::chi) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::chi) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::chi) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::chi) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::chi) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::chi) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::chi) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::chi) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::chi) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::chi))/(dx[1]*dx[2]);
        amrex::Real dDDchi_LL_22 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::chi) - 1.0/12.0*state_fab(i, j, k + 2, Idx::chi) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::chi) - 1.0/12.0*state_fab(i, j, k - 2, Idx::chi) - 5.0/2.0*state_fab(i, j, k, Idx::chi))/std::pow(dx[2], 2);

        amrex::Real dDgamtilde_LLL_000 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::gamtilde_LL_00) - 1.0/12.0*state_fab(i + 2, j, k, Idx::gamtilde_LL_00) - 2.0/3.0*state_fab(i - 1, j, k, Idx::gamtilde_LL_00) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::gamtilde_LL_00))/dx[0];
        amrex::Real dDgamtilde_LLL_001 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::gamtilde_LL_00) - 1.0/12.0*state_fab(i, j + 2, k, Idx::gamtilde_LL_00) - 2.0/3.0*state_fab(i, j - 1, k, Idx::gamtilde_LL_00) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::gamtilde_LL_00))/dx[1];
        amrex::Real dDgamtilde_LLL_002 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::gamtilde_LL_00) - 1.0/12.0*state_fab(i, j, k + 2, Idx::gamtilde_LL_00) - 2.0/3.0*state_fab(i, j, k - 1, Idx::gamtilde_LL_00) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::gamtilde_LL_00))/dx[2];
        amrex::Real dDgamtilde_LLL_010 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::gamtilde_LL_01) - 1.0/12.0*state_fab(i + 2, j, k, Idx::gamtilde_LL_01) - 2.0/3.0*state_fab(i - 1, j, k, Idx::gamtilde_LL_01) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::gamtilde_LL_01))/dx[0];
        amrex::Real dDgamtilde_LLL_011 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::gamtilde_LL_01) - 1.0/12.0*state_fab(i, j + 2, k, Idx::gamtilde_LL_01) - 2.0/3.0*state_fab(i, j - 1, k, Idx::gamtilde_LL_01) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::gamtilde_LL_01))/dx[1];
        amrex::Real dDgamtilde_LLL_012 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::gamtilde_LL_01) - 1.0/12.0*state_fab(i, j, k + 2, Idx::gamtilde_LL_01) - 2.0/3.0*state_fab(i, j, k - 1, Idx::gamtilde_LL_01) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::gamtilde_LL_01))/dx[2];
        amrex::Real dDgamtilde_LLL_020 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::gamtilde_LL_02) - 1.0/12.0*state_fab(i + 2, j, k, Idx::gamtilde_LL_02) - 2.0/3.0*state_fab(i - 1, j, k, Idx::gamtilde_LL_02) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::gamtilde_LL_02))/dx[0];
        amrex::Real dDgamtilde_LLL_021 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::gamtilde_LL_02) - 1.0/12.0*state_fab(i, j + 2, k, Idx::gamtilde_LL_02) - 2.0/3.0*state_fab(i, j - 1, k, Idx::gamtilde_LL_02) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::gamtilde_LL_02))/dx[1];
        amrex::Real dDgamtilde_LLL_022 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::gamtilde_LL_02) - 1.0/12.0*state_fab(i, j, k + 2, Idx::gamtilde_LL_02) - 2.0/3.0*state_fab(i, j, k - 1, Idx::gamtilde_LL_02) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::gamtilde_LL_02))/dx[2];
        amrex::Real dDgamtilde_LLL_110 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::gamtilde_LL_11) - 1.0/12.0*state_fab(i + 2, j, k, Idx::gamtilde_LL_11) - 2.0/3.0*state_fab(i - 1, j, k, Idx::gamtilde_LL_11) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::gamtilde_LL_11))/dx[0];
        amrex::Real dDgamtilde_LLL_111 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::gamtilde_LL_11) - 1.0/12.0*state_fab(i, j + 2, k, Idx::gamtilde_LL_11) - 2.0/3.0*state_fab(i, j - 1, k, Idx::gamtilde_LL_11) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::gamtilde_LL_11))/dx[1];
        amrex::Real dDgamtilde_LLL_112 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::gamtilde_LL_11) - 1.0/12.0*state_fab(i, j, k + 2, Idx::gamtilde_LL_11) - 2.0/3.0*state_fab(i, j, k - 1, Idx::gamtilde_LL_11) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::gamtilde_LL_11))/dx[2];
        amrex::Real dDgamtilde_LLL_120 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::gamtilde_LL_12) - 1.0/12.0*state_fab(i + 2, j, k, Idx::gamtilde_LL_12) - 2.0/3.0*state_fab(i - 1, j, k, Idx::gamtilde_LL_12) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::gamtilde_LL_12))/dx[0];
        amrex::Real dDgamtilde_LLL_121 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::gamtilde_LL_12) - 1.0/12.0*state_fab(i, j + 2, k, Idx::gamtilde_LL_12) - 2.0/3.0*state_fab(i, j - 1, k, Idx::gamtilde_LL_12) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::gamtilde_LL_12))/dx[1];
        amrex::Real dDgamtilde_LLL_122 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::gamtilde_LL_12) - 1.0/12.0*state_fab(i, j, k + 2, Idx::gamtilde_LL_12) - 2.0/3.0*state_fab(i, j, k - 1, Idx::gamtilde_LL_12) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::gamtilde_LL_12))/dx[2];
        amrex::Real dDgamtilde_LLL_220 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::gamtilde_LL_22) - 1.0/12.0*state_fab(i + 2, j, k, Idx::gamtilde_LL_22) - 2.0/3.0*state_fab(i - 1, j, k, Idx::gamtilde_LL_22) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::gamtilde_LL_22))/dx[0];
        amrex::Real dDgamtilde_LLL_221 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::gamtilde_LL_22) - 1.0/12.0*state_fab(i, j + 2, k, Idx::gamtilde_LL_22) - 2.0/3.0*state_fab(i, j - 1, k, Idx::gamtilde_LL_22) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::gamtilde_LL_22))/dx[1];
        amrex::Real dDgamtilde_LLL_222 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::gamtilde_LL_22) - 1.0/12.0*state_fab(i, j, k + 2, Idx::gamtilde_LL_22) - 2.0/3.0*state_fab(i, j, k - 1, Idx::gamtilde_LL_22) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::gamtilde_LL_22))/dx[2];

        amrex::Real dDDgamtilde_LLLL_0000 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::gamtilde_LL_00) - 1.0/12.0*state_fab(i + 2, j, k, Idx::gamtilde_LL_00) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::gamtilde_LL_00) - 1.0/12.0*state_fab(i - 2, j, k, Idx::gamtilde_LL_00) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_00))/std::pow(dx[0], 2);
        amrex::Real dDDgamtilde_LLLL_0001 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::gamtilde_LL_00) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::gamtilde_LL_00) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::gamtilde_LL_00) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::gamtilde_LL_00) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::gamtilde_LL_00) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::gamtilde_LL_00) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::gamtilde_LL_00) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::gamtilde_LL_00) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::gamtilde_LL_00) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::gamtilde_LL_00) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::gamtilde_LL_00) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::gamtilde_LL_00) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::gamtilde_LL_00) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::gamtilde_LL_00) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::gamtilde_LL_00) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::gamtilde_LL_00))/(dx[0]*dx[1]);
        amrex::Real dDDgamtilde_LLLL_0002 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::gamtilde_LL_00) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::gamtilde_LL_00) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::gamtilde_LL_00) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::gamtilde_LL_00) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::gamtilde_LL_00) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::gamtilde_LL_00) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::gamtilde_LL_00) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::gamtilde_LL_00) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::gamtilde_LL_00) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::gamtilde_LL_00) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::gamtilde_LL_00) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::gamtilde_LL_00) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::gamtilde_LL_00) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::gamtilde_LL_00) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::gamtilde_LL_00) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::gamtilde_LL_00))/(dx[0]*dx[2]);
        amrex::Real dDDgamtilde_LLLL_0011 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::gamtilde_LL_00) - 1.0/12.0*state_fab(i, j + 2, k, Idx::gamtilde_LL_00) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::gamtilde_LL_00) - 1.0/12.0*state_fab(i, j - 2, k, Idx::gamtilde_LL_00) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_00))/std::pow(dx[1], 2);
        amrex::Real dDDgamtilde_LLLL_0012 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::gamtilde_LL_00) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::gamtilde_LL_00) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::gamtilde_LL_00) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::gamtilde_LL_00) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::gamtilde_LL_00) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::gamtilde_LL_00) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::gamtilde_LL_00) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::gamtilde_LL_00) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::gamtilde_LL_00) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::gamtilde_LL_00) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::gamtilde_LL_00) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::gamtilde_LL_00) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::gamtilde_LL_00) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::gamtilde_LL_00) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::gamtilde_LL_00) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::gamtilde_LL_00))/(dx[1]*dx[2]);
        amrex::Real dDDgamtilde_LLLL_0022 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::gamtilde_LL_00) - 1.0/12.0*state_fab(i, j, k + 2, Idx::gamtilde_LL_00) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::gamtilde_LL_00) - 1.0/12.0*state_fab(i, j, k - 2, Idx::gamtilde_LL_00) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_00))/std::pow(dx[2], 2);
        amrex::Real dDDgamtilde_LLLL_0100 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::gamtilde_LL_01) - 1.0/12.0*state_fab(i + 2, j, k, Idx::gamtilde_LL_01) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::gamtilde_LL_01) - 1.0/12.0*state_fab(i - 2, j, k, Idx::gamtilde_LL_01) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_01))/std::pow(dx[0], 2);
        amrex::Real dDDgamtilde_LLLL_0101 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::gamtilde_LL_01) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::gamtilde_LL_01) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::gamtilde_LL_01) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::gamtilde_LL_01) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::gamtilde_LL_01) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::gamtilde_LL_01) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::gamtilde_LL_01) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::gamtilde_LL_01) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::gamtilde_LL_01) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::gamtilde_LL_01) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::gamtilde_LL_01) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::gamtilde_LL_01) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::gamtilde_LL_01) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::gamtilde_LL_01) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::gamtilde_LL_01) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::gamtilde_LL_01))/(dx[0]*dx[1]);
        amrex::Real dDDgamtilde_LLLL_0102 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::gamtilde_LL_01) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::gamtilde_LL_01) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::gamtilde_LL_01) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::gamtilde_LL_01) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::gamtilde_LL_01) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::gamtilde_LL_01) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::gamtilde_LL_01) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::gamtilde_LL_01) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::gamtilde_LL_01) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::gamtilde_LL_01) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::gamtilde_LL_01) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::gamtilde_LL_01) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::gamtilde_LL_01) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::gamtilde_LL_01) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::gamtilde_LL_01) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::gamtilde_LL_01))/(dx[0]*dx[2]);
        amrex::Real dDDgamtilde_LLLL_0111 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::gamtilde_LL_01) - 1.0/12.0*state_fab(i, j + 2, k, Idx::gamtilde_LL_01) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::gamtilde_LL_01) - 1.0/12.0*state_fab(i, j - 2, k, Idx::gamtilde_LL_01) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_01))/std::pow(dx[1], 2);
        amrex::Real dDDgamtilde_LLLL_0112 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::gamtilde_LL_01) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::gamtilde_LL_01) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::gamtilde_LL_01) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::gamtilde_LL_01) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::gamtilde_LL_01) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::gamtilde_LL_01) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::gamtilde_LL_01) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::gamtilde_LL_01) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::gamtilde_LL_01) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::gamtilde_LL_01) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::gamtilde_LL_01) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::gamtilde_LL_01) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::gamtilde_LL_01) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::gamtilde_LL_01) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::gamtilde_LL_01) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::gamtilde_LL_01))/(dx[1]*dx[2]);
        amrex::Real dDDgamtilde_LLLL_0122 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::gamtilde_LL_01) - 1.0/12.0*state_fab(i, j, k + 2, Idx::gamtilde_LL_01) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::gamtilde_LL_01) - 1.0/12.0*state_fab(i, j, k - 2, Idx::gamtilde_LL_01) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_01))/std::pow(dx[2], 2);
        amrex::Real dDDgamtilde_LLLL_0200 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::gamtilde_LL_02) - 1.0/12.0*state_fab(i + 2, j, k, Idx::gamtilde_LL_02) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::gamtilde_LL_02) - 1.0/12.0*state_fab(i - 2, j, k, Idx::gamtilde_LL_02) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_02))/std::pow(dx[0], 2);
        amrex::Real dDDgamtilde_LLLL_0201 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::gamtilde_LL_02) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::gamtilde_LL_02) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::gamtilde_LL_02) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::gamtilde_LL_02) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::gamtilde_LL_02) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::gamtilde_LL_02) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::gamtilde_LL_02) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::gamtilde_LL_02) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::gamtilde_LL_02) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::gamtilde_LL_02) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::gamtilde_LL_02) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::gamtilde_LL_02) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::gamtilde_LL_02) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::gamtilde_LL_02) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::gamtilde_LL_02) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::gamtilde_LL_02))/(dx[0]*dx[1]);
        amrex::Real dDDgamtilde_LLLL_0202 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::gamtilde_LL_02) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::gamtilde_LL_02) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::gamtilde_LL_02) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::gamtilde_LL_02) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::gamtilde_LL_02) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::gamtilde_LL_02) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::gamtilde_LL_02) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::gamtilde_LL_02) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::gamtilde_LL_02) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::gamtilde_LL_02) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::gamtilde_LL_02) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::gamtilde_LL_02) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::gamtilde_LL_02) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::gamtilde_LL_02) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::gamtilde_LL_02) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::gamtilde_LL_02))/(dx[0]*dx[2]);
        amrex::Real dDDgamtilde_LLLL_0211 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::gamtilde_LL_02) - 1.0/12.0*state_fab(i, j + 2, k, Idx::gamtilde_LL_02) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::gamtilde_LL_02) - 1.0/12.0*state_fab(i, j - 2, k, Idx::gamtilde_LL_02) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_02))/std::pow(dx[1], 2);
        amrex::Real dDDgamtilde_LLLL_0212 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::gamtilde_LL_02) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::gamtilde_LL_02) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::gamtilde_LL_02) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::gamtilde_LL_02) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::gamtilde_LL_02) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::gamtilde_LL_02) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::gamtilde_LL_02) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::gamtilde_LL_02) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::gamtilde_LL_02) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::gamtilde_LL_02) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::gamtilde_LL_02) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::gamtilde_LL_02) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::gamtilde_LL_02) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::gamtilde_LL_02) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::gamtilde_LL_02) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::gamtilde_LL_02))/(dx[1]*dx[2]);
        amrex::Real dDDgamtilde_LLLL_0222 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::gamtilde_LL_02) - 1.0/12.0*state_fab(i, j, k + 2, Idx::gamtilde_LL_02) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::gamtilde_LL_02) - 1.0/12.0*state_fab(i, j, k - 2, Idx::gamtilde_LL_02) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_02))/std::pow(dx[2], 2);
        amrex::Real dDDgamtilde_LLLL_1100 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::gamtilde_LL_11) - 1.0/12.0*state_fab(i + 2, j, k, Idx::gamtilde_LL_11) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::gamtilde_LL_11) - 1.0/12.0*state_fab(i - 2, j, k, Idx::gamtilde_LL_11) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_11))/std::pow(dx[0], 2);
        amrex::Real dDDgamtilde_LLLL_1101 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::gamtilde_LL_11) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::gamtilde_LL_11) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::gamtilde_LL_11) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::gamtilde_LL_11) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::gamtilde_LL_11) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::gamtilde_LL_11) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::gamtilde_LL_11) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::gamtilde_LL_11) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::gamtilde_LL_11) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::gamtilde_LL_11) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::gamtilde_LL_11) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::gamtilde_LL_11) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::gamtilde_LL_11) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::gamtilde_LL_11) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::gamtilde_LL_11) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::gamtilde_LL_11))/(dx[0]*dx[1]);
        amrex::Real dDDgamtilde_LLLL_1102 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::gamtilde_LL_11) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::gamtilde_LL_11) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::gamtilde_LL_11) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::gamtilde_LL_11) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::gamtilde_LL_11) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::gamtilde_LL_11) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::gamtilde_LL_11) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::gamtilde_LL_11) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::gamtilde_LL_11) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::gamtilde_LL_11) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::gamtilde_LL_11) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::gamtilde_LL_11) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::gamtilde_LL_11) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::gamtilde_LL_11) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::gamtilde_LL_11) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::gamtilde_LL_11))/(dx[0]*dx[2]);
        amrex::Real dDDgamtilde_LLLL_1111 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::gamtilde_LL_11) - 1.0/12.0*state_fab(i, j + 2, k, Idx::gamtilde_LL_11) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::gamtilde_LL_11) - 1.0/12.0*state_fab(i, j - 2, k, Idx::gamtilde_LL_11) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_11))/std::pow(dx[1], 2);
        amrex::Real dDDgamtilde_LLLL_1112 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::gamtilde_LL_11) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::gamtilde_LL_11) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::gamtilde_LL_11) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::gamtilde_LL_11) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::gamtilde_LL_11) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::gamtilde_LL_11) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::gamtilde_LL_11) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::gamtilde_LL_11) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::gamtilde_LL_11) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::gamtilde_LL_11) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::gamtilde_LL_11) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::gamtilde_LL_11) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::gamtilde_LL_11) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::gamtilde_LL_11) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::gamtilde_LL_11) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::gamtilde_LL_11))/(dx[1]*dx[2]);
        amrex::Real dDDgamtilde_LLLL_1122 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::gamtilde_LL_11) - 1.0/12.0*state_fab(i, j, k + 2, Idx::gamtilde_LL_11) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::gamtilde_LL_11) - 1.0/12.0*state_fab(i, j, k - 2, Idx::gamtilde_LL_11) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_11))/std::pow(dx[2], 2);
        amrex::Real dDDgamtilde_LLLL_1200 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::gamtilde_LL_12) - 1.0/12.0*state_fab(i + 2, j, k, Idx::gamtilde_LL_12) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::gamtilde_LL_12) - 1.0/12.0*state_fab(i - 2, j, k, Idx::gamtilde_LL_12) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_12))/std::pow(dx[0], 2);
        amrex::Real dDDgamtilde_LLLL_1201 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::gamtilde_LL_12) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::gamtilde_LL_12) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::gamtilde_LL_12) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::gamtilde_LL_12) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::gamtilde_LL_12) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::gamtilde_LL_12) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::gamtilde_LL_12) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::gamtilde_LL_12) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::gamtilde_LL_12) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::gamtilde_LL_12) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::gamtilde_LL_12) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::gamtilde_LL_12) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::gamtilde_LL_12) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::gamtilde_LL_12) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::gamtilde_LL_12) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::gamtilde_LL_12))/(dx[0]*dx[1]);
        amrex::Real dDDgamtilde_LLLL_1202 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::gamtilde_LL_12) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::gamtilde_LL_12) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::gamtilde_LL_12) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::gamtilde_LL_12) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::gamtilde_LL_12) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::gamtilde_LL_12) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::gamtilde_LL_12) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::gamtilde_LL_12) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::gamtilde_LL_12) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::gamtilde_LL_12) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::gamtilde_LL_12) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::gamtilde_LL_12) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::gamtilde_LL_12) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::gamtilde_LL_12) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::gamtilde_LL_12) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::gamtilde_LL_12))/(dx[0]*dx[2]);
        amrex::Real dDDgamtilde_LLLL_1211 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::gamtilde_LL_12) - 1.0/12.0*state_fab(i, j + 2, k, Idx::gamtilde_LL_12) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::gamtilde_LL_12) - 1.0/12.0*state_fab(i, j - 2, k, Idx::gamtilde_LL_12) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_12))/std::pow(dx[1], 2);
        amrex::Real dDDgamtilde_LLLL_1212 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::gamtilde_LL_12) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::gamtilde_LL_12) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::gamtilde_LL_12) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::gamtilde_LL_12) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::gamtilde_LL_12) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::gamtilde_LL_12) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::gamtilde_LL_12) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::gamtilde_LL_12) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::gamtilde_LL_12) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::gamtilde_LL_12) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::gamtilde_LL_12) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::gamtilde_LL_12) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::gamtilde_LL_12) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::gamtilde_LL_12) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::gamtilde_LL_12) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::gamtilde_LL_12))/(dx[1]*dx[2]);
        amrex::Real dDDgamtilde_LLLL_1222 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::gamtilde_LL_12) - 1.0/12.0*state_fab(i, j, k + 2, Idx::gamtilde_LL_12) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::gamtilde_LL_12) - 1.0/12.0*state_fab(i, j, k - 2, Idx::gamtilde_LL_12) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_12))/std::pow(dx[2], 2);
        amrex::Real dDDgamtilde_LLLL_2200 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::gamtilde_LL_22) - 1.0/12.0*state_fab(i + 2, j, k, Idx::gamtilde_LL_22) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::gamtilde_LL_22) - 1.0/12.0*state_fab(i - 2, j, k, Idx::gamtilde_LL_22) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_22))/std::pow(dx[0], 2);
        amrex::Real dDDgamtilde_LLLL_2201 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::gamtilde_LL_22) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::gamtilde_LL_22) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::gamtilde_LL_22) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::gamtilde_LL_22) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::gamtilde_LL_22) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::gamtilde_LL_22) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::gamtilde_LL_22) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::gamtilde_LL_22) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::gamtilde_LL_22) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::gamtilde_LL_22) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::gamtilde_LL_22) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::gamtilde_LL_22) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::gamtilde_LL_22) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::gamtilde_LL_22) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::gamtilde_LL_22) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::gamtilde_LL_22))/(dx[0]*dx[1]);
        amrex::Real dDDgamtilde_LLLL_2202 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::gamtilde_LL_22) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::gamtilde_LL_22) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::gamtilde_LL_22) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::gamtilde_LL_22) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::gamtilde_LL_22) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::gamtilde_LL_22) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::gamtilde_LL_22) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::gamtilde_LL_22) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::gamtilde_LL_22) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::gamtilde_LL_22) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::gamtilde_LL_22) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::gamtilde_LL_22) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::gamtilde_LL_22) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::gamtilde_LL_22) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::gamtilde_LL_22) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::gamtilde_LL_22))/(dx[0]*dx[2]);
        amrex::Real dDDgamtilde_LLLL_2211 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::gamtilde_LL_22) - 1.0/12.0*state_fab(i, j + 2, k, Idx::gamtilde_LL_22) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::gamtilde_LL_22) - 1.0/12.0*state_fab(i, j - 2, k, Idx::gamtilde_LL_22) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_22))/std::pow(dx[1], 2);
        amrex::Real dDDgamtilde_LLLL_2212 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::gamtilde_LL_22) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::gamtilde_LL_22) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::gamtilde_LL_22) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::gamtilde_LL_22) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::gamtilde_LL_22) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::gamtilde_LL_22) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::gamtilde_LL_22) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::gamtilde_LL_22) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::gamtilde_LL_22) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::gamtilde_LL_22) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::gamtilde_LL_22) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::gamtilde_LL_22) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::gamtilde_LL_22) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::gamtilde_LL_22) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::gamtilde_LL_22) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::gamtilde_LL_22))/(dx[1]*dx[2]);
        amrex::Real dDDgamtilde_LLLL_2222 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::gamtilde_LL_22) - 1.0/12.0*state_fab(i, j, k + 2, Idx::gamtilde_LL_22) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::gamtilde_LL_22) - 1.0/12.0*state_fab(i, j, k - 2, Idx::gamtilde_LL_22) - 5.0/2.0*state_fab(i, j, k, Idx::gamtilde_LL_22))/std::pow(dx[2], 2);

        amrex::Real dDKhat_L_0 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::Khat) - 1.0/12.0*state_fab(i + 2, j, k, Idx::Khat) - 2.0/3.0*state_fab(i - 1, j, k, Idx::Khat) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::Khat))/dx[0];
        amrex::Real dDKhat_L_1 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::Khat) - 1.0/12.0*state_fab(i, j + 2, k, Idx::Khat) - 2.0/3.0*state_fab(i, j - 1, k, Idx::Khat) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::Khat))/dx[1];
        amrex::Real dDKhat_L_2 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::Khat) - 1.0/12.0*state_fab(i, j, k + 2, Idx::Khat) - 2.0/3.0*state_fab(i, j, k - 1, Idx::Khat) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::Khat))/dx[2];

        amrex::Real dDAtilde_LLL_000 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::Atilde_LL_00) - 1.0/12.0*state_fab(i + 2, j, k, Idx::Atilde_LL_00) - 2.0/3.0*state_fab(i - 1, j, k, Idx::Atilde_LL_00) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::Atilde_LL_00))/dx[0];
        amrex::Real dDAtilde_LLL_001 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::Atilde_LL_00) - 1.0/12.0*state_fab(i, j + 2, k, Idx::Atilde_LL_00) - 2.0/3.0*state_fab(i, j - 1, k, Idx::Atilde_LL_00) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::Atilde_LL_00))/dx[1];
        amrex::Real dDAtilde_LLL_002 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::Atilde_LL_00) - 1.0/12.0*state_fab(i, j, k + 2, Idx::Atilde_LL_00) - 2.0/3.0*state_fab(i, j, k - 1, Idx::Atilde_LL_00) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::Atilde_LL_00))/dx[2];
        amrex::Real dDAtilde_LLL_010 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::Atilde_LL_01) - 1.0/12.0*state_fab(i + 2, j, k, Idx::Atilde_LL_01) - 2.0/3.0*state_fab(i - 1, j, k, Idx::Atilde_LL_01) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::Atilde_LL_01))/dx[0];
        amrex::Real dDAtilde_LLL_011 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::Atilde_LL_01) - 1.0/12.0*state_fab(i, j + 2, k, Idx::Atilde_LL_01) - 2.0/3.0*state_fab(i, j - 1, k, Idx::Atilde_LL_01) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::Atilde_LL_01))/dx[1];
        amrex::Real dDAtilde_LLL_012 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::Atilde_LL_01) - 1.0/12.0*state_fab(i, j, k + 2, Idx::Atilde_LL_01) - 2.0/3.0*state_fab(i, j, k - 1, Idx::Atilde_LL_01) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::Atilde_LL_01))/dx[2];
        amrex::Real dDAtilde_LLL_020 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::Atilde_LL_02) - 1.0/12.0*state_fab(i + 2, j, k, Idx::Atilde_LL_02) - 2.0/3.0*state_fab(i - 1, j, k, Idx::Atilde_LL_02) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::Atilde_LL_02))/dx[0];
        amrex::Real dDAtilde_LLL_021 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::Atilde_LL_02) - 1.0/12.0*state_fab(i, j + 2, k, Idx::Atilde_LL_02) - 2.0/3.0*state_fab(i, j - 1, k, Idx::Atilde_LL_02) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::Atilde_LL_02))/dx[1];
        amrex::Real dDAtilde_LLL_022 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::Atilde_LL_02) - 1.0/12.0*state_fab(i, j, k + 2, Idx::Atilde_LL_02) - 2.0/3.0*state_fab(i, j, k - 1, Idx::Atilde_LL_02) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::Atilde_LL_02))/dx[2];
        amrex::Real dDAtilde_LLL_110 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::Atilde_LL_11) - 1.0/12.0*state_fab(i + 2, j, k, Idx::Atilde_LL_11) - 2.0/3.0*state_fab(i - 1, j, k, Idx::Atilde_LL_11) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::Atilde_LL_11))/dx[0];
        amrex::Real dDAtilde_LLL_111 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::Atilde_LL_11) - 1.0/12.0*state_fab(i, j + 2, k, Idx::Atilde_LL_11) - 2.0/3.0*state_fab(i, j - 1, k, Idx::Atilde_LL_11) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::Atilde_LL_11))/dx[1];
        amrex::Real dDAtilde_LLL_112 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::Atilde_LL_11) - 1.0/12.0*state_fab(i, j, k + 2, Idx::Atilde_LL_11) - 2.0/3.0*state_fab(i, j, k - 1, Idx::Atilde_LL_11) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::Atilde_LL_11))/dx[2];
        amrex::Real dDAtilde_LLL_120 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::Atilde_LL_12) - 1.0/12.0*state_fab(i + 2, j, k, Idx::Atilde_LL_12) - 2.0/3.0*state_fab(i - 1, j, k, Idx::Atilde_LL_12) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::Atilde_LL_12))/dx[0];
        amrex::Real dDAtilde_LLL_121 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::Atilde_LL_12) - 1.0/12.0*state_fab(i, j + 2, k, Idx::Atilde_LL_12) - 2.0/3.0*state_fab(i, j - 1, k, Idx::Atilde_LL_12) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::Atilde_LL_12))/dx[1];
        amrex::Real dDAtilde_LLL_122 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::Atilde_LL_12) - 1.0/12.0*state_fab(i, j, k + 2, Idx::Atilde_LL_12) - 2.0/3.0*state_fab(i, j, k - 1, Idx::Atilde_LL_12) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::Atilde_LL_12))/dx[2];
        amrex::Real dDAtilde_LLL_220 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::Atilde_LL_22) - 1.0/12.0*state_fab(i + 2, j, k, Idx::Atilde_LL_22) - 2.0/3.0*state_fab(i - 1, j, k, Idx::Atilde_LL_22) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::Atilde_LL_22))/dx[0];
        amrex::Real dDAtilde_LLL_221 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::Atilde_LL_22) - 1.0/12.0*state_fab(i, j + 2, k, Idx::Atilde_LL_22) - 2.0/3.0*state_fab(i, j - 1, k, Idx::Atilde_LL_22) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::Atilde_LL_22))/dx[1];
        amrex::Real dDAtilde_LLL_222 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::Atilde_LL_22) - 1.0/12.0*state_fab(i, j, k + 2, Idx::Atilde_LL_22) - 2.0/3.0*state_fab(i, j, k - 1, Idx::Atilde_LL_22) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::Atilde_LL_22))/dx[2];

        amrex::Real dDtheta_L_0 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::theta) - 1.0/12.0*state_fab(i + 2, j, k, Idx::theta) - 2.0/3.0*state_fab(i - 1, j, k, Idx::theta) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::theta))/dx[0];
        amrex::Real dDtheta_L_1 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::theta) - 1.0/12.0*state_fab(i, j + 2, k, Idx::theta) - 2.0/3.0*state_fab(i, j - 1, k, Idx::theta) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::theta))/dx[1];
        amrex::Real dDtheta_L_2 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::theta) - 1.0/12.0*state_fab(i, j, k + 2, Idx::theta) - 2.0/3.0*state_fab(i, j, k - 1, Idx::theta) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::theta))/dx[2];

        amrex::Real invgamtilde_UU_00 = (gamtilde_LL_11*gamtilde_LL_22 - std::pow(gamtilde_LL_12, 2))/(gamtilde_LL_00*gamtilde_LL_11*gamtilde_LL_22 - gamtilde_LL_00*std::pow(gamtilde_LL_12, 2) - std::pow(gamtilde_LL_01, 2)*gamtilde_LL_22 + 2*gamtilde_LL_01*gamtilde_LL_02*gamtilde_LL_12 - std::pow(gamtilde_LL_02, 2)*gamtilde_LL_11);
        amrex::Real invgamtilde_UU_01 = (-gamtilde_LL_01*gamtilde_LL_22 + gamtilde_LL_02*gamtilde_LL_12)/(gamtilde_LL_00*gamtilde_LL_11*gamtilde_LL_22 - gamtilde_LL_00*std::pow(gamtilde_LL_12, 2) - std::pow(gamtilde_LL_01, 2)*gamtilde_LL_22 + 2*gamtilde_LL_01*gamtilde_LL_02*gamtilde_LL_12 - std::pow(gamtilde_LL_02, 2)*gamtilde_LL_11);
        amrex::Real invgamtilde_UU_02 = (gamtilde_LL_01*gamtilde_LL_12 - gamtilde_LL_02*gamtilde_LL_11)/(gamtilde_LL_00*gamtilde_LL_11*gamtilde_LL_22 - gamtilde_LL_00*std::pow(gamtilde_LL_12, 2) - std::pow(gamtilde_LL_01, 2)*gamtilde_LL_22 + 2*gamtilde_LL_01*gamtilde_LL_02*gamtilde_LL_12 - std::pow(gamtilde_LL_02, 2)*gamtilde_LL_11);
        amrex::Real invgamtilde_UU_11 = (gamtilde_LL_00*gamtilde_LL_22 - std::pow(gamtilde_LL_02, 2))/(gamtilde_LL_00*gamtilde_LL_11*gamtilde_LL_22 - gamtilde_LL_00*std::pow(gamtilde_LL_12, 2) - std::pow(gamtilde_LL_01, 2)*gamtilde_LL_22 + 2*gamtilde_LL_01*gamtilde_LL_02*gamtilde_LL_12 - std::pow(gamtilde_LL_02, 2)*gamtilde_LL_11);
        amrex::Real invgamtilde_UU_12 = (-gamtilde_LL_00*gamtilde_LL_12 + gamtilde_LL_01*gamtilde_LL_02)/(gamtilde_LL_00*gamtilde_LL_11*gamtilde_LL_22 - gamtilde_LL_00*std::pow(gamtilde_LL_12, 2) - std::pow(gamtilde_LL_01, 2)*gamtilde_LL_22 + 2*gamtilde_LL_01*gamtilde_LL_02*gamtilde_LL_12 - std::pow(gamtilde_LL_02, 2)*gamtilde_LL_11);
        amrex::Real invgamtilde_UU_22 = (gamtilde_LL_00*gamtilde_LL_11 - std::pow(gamtilde_LL_01, 2))/(gamtilde_LL_00*gamtilde_LL_11*gamtilde_LL_22 - gamtilde_LL_00*std::pow(gamtilde_LL_12, 2) - std::pow(gamtilde_LL_01, 2)*gamtilde_LL_22 + 2*gamtilde_LL_01*gamtilde_LL_02*gamtilde_LL_12 - std::pow(gamtilde_LL_02, 2)*gamtilde_LL_11);

        amrex::Real Christoffelgamtilde_LLL_000 = 0.5*dDgamtilde_LLL_000;
        amrex::Real Christoffelgamtilde_LLL_001 = 0.5*dDgamtilde_LLL_001;
        amrex::Real Christoffelgamtilde_LLL_002 = 0.5*dDgamtilde_LLL_002;
        amrex::Real Christoffelgamtilde_LLL_011 = 1.0*dDgamtilde_LLL_011 - 0.5*dDgamtilde_LLL_110;
        amrex::Real Christoffelgamtilde_LLL_012 = 0.5*dDgamtilde_LLL_012 + 0.5*dDgamtilde_LLL_021 - 0.5*dDgamtilde_LLL_120;
        amrex::Real Christoffelgamtilde_LLL_022 = 1.0*dDgamtilde_LLL_022 - 0.5*dDgamtilde_LLL_220;
        amrex::Real Christoffelgamtilde_LLL_100 = -0.5*dDgamtilde_LLL_001 + 1.0*dDgamtilde_LLL_010;
        amrex::Real Christoffelgamtilde_LLL_101 = 0.5*dDgamtilde_LLL_110;
        amrex::Real Christoffelgamtilde_LLL_102 = 0.5*dDgamtilde_LLL_012 - 0.5*dDgamtilde_LLL_021 + 0.5*dDgamtilde_LLL_120;
        amrex::Real Christoffelgamtilde_LLL_111 = 0.5*dDgamtilde_LLL_111;
        amrex::Real Christoffelgamtilde_LLL_112 = 0.5*dDgamtilde_LLL_112;
        amrex::Real Christoffelgamtilde_LLL_122 = 1.0*dDgamtilde_LLL_122 - 0.5*dDgamtilde_LLL_221;
        amrex::Real Christoffelgamtilde_LLL_200 = -0.5*dDgamtilde_LLL_002 + 1.0*dDgamtilde_LLL_020;
        amrex::Real Christoffelgamtilde_LLL_201 = -0.5*dDgamtilde_LLL_012 + 0.5*dDgamtilde_LLL_021 + 0.5*dDgamtilde_LLL_120;
        amrex::Real Christoffelgamtilde_LLL_202 = 0.5*dDgamtilde_LLL_220;
        amrex::Real Christoffelgamtilde_LLL_211 = -0.5*dDgamtilde_LLL_112 + 1.0*dDgamtilde_LLL_121;
        amrex::Real Christoffelgamtilde_LLL_212 = 0.5*dDgamtilde_LLL_221;
        amrex::Real Christoffelgamtilde_LLL_222 = 0.5*dDgamtilde_LLL_222;

        amrex::Real Christoffelgamtilde_ULL_000 = Christoffelgamtilde_LLL_000*invgamtilde_UU_00 + Christoffelgamtilde_LLL_100*invgamtilde_UU_01 + Christoffelgamtilde_LLL_200*invgamtilde_UU_02;
        amrex::Real Christoffelgamtilde_ULL_001 = Christoffelgamtilde_LLL_001*invgamtilde_UU_00 + Christoffelgamtilde_LLL_101*invgamtilde_UU_01 + Christoffelgamtilde_LLL_201*invgamtilde_UU_02;
        amrex::Real Christoffelgamtilde_ULL_002 = Christoffelgamtilde_LLL_002*invgamtilde_UU_00 + Christoffelgamtilde_LLL_102*invgamtilde_UU_01 + Christoffelgamtilde_LLL_202*invgamtilde_UU_02;
        amrex::Real Christoffelgamtilde_ULL_011 = Christoffelgamtilde_LLL_011*invgamtilde_UU_00 + Christoffelgamtilde_LLL_111*invgamtilde_UU_01 + Christoffelgamtilde_LLL_211*invgamtilde_UU_02;
        amrex::Real Christoffelgamtilde_ULL_012 = Christoffelgamtilde_LLL_012*invgamtilde_UU_00 + Christoffelgamtilde_LLL_112*invgamtilde_UU_01 + Christoffelgamtilde_LLL_212*invgamtilde_UU_02;
        amrex::Real Christoffelgamtilde_ULL_022 = Christoffelgamtilde_LLL_022*invgamtilde_UU_00 + Christoffelgamtilde_LLL_122*invgamtilde_UU_01 + Christoffelgamtilde_LLL_222*invgamtilde_UU_02;
        amrex::Real Christoffelgamtilde_ULL_100 = Christoffelgamtilde_LLL_000*invgamtilde_UU_01 + Christoffelgamtilde_LLL_100*invgamtilde_UU_11 + Christoffelgamtilde_LLL_200*invgamtilde_UU_12;
        amrex::Real Christoffelgamtilde_ULL_101 = Christoffelgamtilde_LLL_001*invgamtilde_UU_01 + Christoffelgamtilde_LLL_101*invgamtilde_UU_11 + Christoffelgamtilde_LLL_201*invgamtilde_UU_12;
        amrex::Real Christoffelgamtilde_ULL_102 = Christoffelgamtilde_LLL_002*invgamtilde_UU_01 + Christoffelgamtilde_LLL_102*invgamtilde_UU_11 + Christoffelgamtilde_LLL_202*invgamtilde_UU_12;
        amrex::Real Christoffelgamtilde_ULL_111 = Christoffelgamtilde_LLL_011*invgamtilde_UU_01 + Christoffelgamtilde_LLL_111*invgamtilde_UU_11 + Christoffelgamtilde_LLL_211*invgamtilde_UU_12;
        amrex::Real Christoffelgamtilde_ULL_112 = Christoffelgamtilde_LLL_012*invgamtilde_UU_01 + Christoffelgamtilde_LLL_112*invgamtilde_UU_11 + Christoffelgamtilde_LLL_212*invgamtilde_UU_12;
        amrex::Real Christoffelgamtilde_ULL_122 = Christoffelgamtilde_LLL_022*invgamtilde_UU_01 + Christoffelgamtilde_LLL_122*invgamtilde_UU_11 + Christoffelgamtilde_LLL_222*invgamtilde_UU_12;
        amrex::Real Christoffelgamtilde_ULL_200 = Christoffelgamtilde_LLL_000*invgamtilde_UU_02 + Christoffelgamtilde_LLL_100*invgamtilde_UU_12 + Christoffelgamtilde_LLL_200*invgamtilde_UU_22;
        amrex::Real Christoffelgamtilde_ULL_201 = Christoffelgamtilde_LLL_001*invgamtilde_UU_02 + Christoffelgamtilde_LLL_101*invgamtilde_UU_12 + Christoffelgamtilde_LLL_201*invgamtilde_UU_22;
        amrex::Real Christoffelgamtilde_ULL_202 = Christoffelgamtilde_LLL_002*invgamtilde_UU_02 + Christoffelgamtilde_LLL_102*invgamtilde_UU_12 + Christoffelgamtilde_LLL_202*invgamtilde_UU_22;
        amrex::Real Christoffelgamtilde_ULL_211 = Christoffelgamtilde_LLL_011*invgamtilde_UU_02 + Christoffelgamtilde_LLL_111*invgamtilde_UU_12 + Christoffelgamtilde_LLL_211*invgamtilde_UU_22;
        amrex::Real Christoffelgamtilde_ULL_212 = Christoffelgamtilde_LLL_012*invgamtilde_UU_02 + Christoffelgamtilde_LLL_112*invgamtilde_UU_12 + Christoffelgamtilde_LLL_212*invgamtilde_UU_22;
        amrex::Real Christoffelgamtilde_ULL_222 = Christoffelgamtilde_LLL_022*invgamtilde_UU_02 + Christoffelgamtilde_LLL_122*invgamtilde_UU_12 + Christoffelgamtilde_LLL_222*invgamtilde_UU_22;

        amrex::Real gam_LL_00 = gamtilde_LL_00/chi;
        amrex::Real gam_LL_01 = gamtilde_LL_01/chi;
        amrex::Real gam_LL_02 = gamtilde_LL_02/chi;
        amrex::Real gam_LL_10 = gamtilde_LL_01/chi;
        amrex::Real gam_LL_11 = gamtilde_LL_11/chi;
        amrex::Real gam_LL_12 = gamtilde_LL_12/chi;
        amrex::Real gam_LL_20 = gamtilde_LL_02/chi;
        amrex::Real gam_LL_21 = gamtilde_LL_12/chi;
        amrex::Real gam_LL_22 = gamtilde_LL_22/chi;

        amrex::Real invgam_UU_00 = chi*invgamtilde_UU_00;
        amrex::Real invgam_UU_01 = chi*invgamtilde_UU_01;
        amrex::Real invgam_UU_02 = chi*invgamtilde_UU_02;
        amrex::Real invgam_UU_10 = chi*invgamtilde_UU_01;
        amrex::Real invgam_UU_11 = chi*invgamtilde_UU_11;
        amrex::Real invgam_UU_12 = chi*invgamtilde_UU_12;
        amrex::Real invgam_UU_20 = chi*invgamtilde_UU_02;
        amrex::Real invgam_UU_21 = chi*invgamtilde_UU_12;
        amrex::Real invgam_UU_22 = chi*invgamtilde_UU_22;

        amrex::Real ConformalChristoffelchigamtilde_ULL_000 = Christoffelgamtilde_ULL_000 + 0.5*dDchi_L_0*gamtilde_LL_00*invgamtilde_UU_00/chi - 1.0*dDchi_L_0/chi + 0.5*dDchi_L_1*gamtilde_LL_00*invgamtilde_UU_01/chi + 0.5*dDchi_L_2*gamtilde_LL_00*invgamtilde_UU_02/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_001 = Christoffelgamtilde_ULL_001 + 0.5*dDchi_L_0*gamtilde_LL_01*invgamtilde_UU_00/chi + 0.5*dDchi_L_1*gamtilde_LL_01*invgamtilde_UU_01/chi - 0.5*dDchi_L_1/chi + 0.5*dDchi_L_2*gamtilde_LL_01*invgamtilde_UU_02/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_002 = Christoffelgamtilde_ULL_002 + 0.5*dDchi_L_0*gamtilde_LL_02*invgamtilde_UU_00/chi + 0.5*dDchi_L_1*gamtilde_LL_02*invgamtilde_UU_01/chi + 0.5*dDchi_L_2*gamtilde_LL_02*invgamtilde_UU_02/chi - 0.5*dDchi_L_2/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_011 = Christoffelgamtilde_ULL_011 + 0.5*dDchi_L_0*gamtilde_LL_11*invgamtilde_UU_00/chi + 0.5*dDchi_L_1*gamtilde_LL_11*invgamtilde_UU_01/chi + 0.5*dDchi_L_2*gamtilde_LL_11*invgamtilde_UU_02/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_012 = Christoffelgamtilde_ULL_012 + 0.5*dDchi_L_0*gamtilde_LL_12*invgamtilde_UU_00/chi + 0.5*dDchi_L_1*gamtilde_LL_12*invgamtilde_UU_01/chi + 0.5*dDchi_L_2*gamtilde_LL_12*invgamtilde_UU_02/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_022 = Christoffelgamtilde_ULL_022 + 0.5*dDchi_L_0*gamtilde_LL_22*invgamtilde_UU_00/chi + 0.5*dDchi_L_1*gamtilde_LL_22*invgamtilde_UU_01/chi + 0.5*dDchi_L_2*gamtilde_LL_22*invgamtilde_UU_02/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_100 = Christoffelgamtilde_ULL_100 + 0.5*dDchi_L_0*gamtilde_LL_00*invgamtilde_UU_01/chi + 0.5*dDchi_L_1*gamtilde_LL_00*invgamtilde_UU_11/chi + 0.5*dDchi_L_2*gamtilde_LL_00*invgamtilde_UU_12/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_101 = Christoffelgamtilde_ULL_101 + 0.5*dDchi_L_0*gamtilde_LL_01*invgamtilde_UU_01/chi - 0.5*dDchi_L_0/chi + 0.5*dDchi_L_1*gamtilde_LL_01*invgamtilde_UU_11/chi + 0.5*dDchi_L_2*gamtilde_LL_01*invgamtilde_UU_12/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_102 = Christoffelgamtilde_ULL_102 + 0.5*dDchi_L_0*gamtilde_LL_02*invgamtilde_UU_01/chi + 0.5*dDchi_L_1*gamtilde_LL_02*invgamtilde_UU_11/chi + 0.5*dDchi_L_2*gamtilde_LL_02*invgamtilde_UU_12/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_111 = Christoffelgamtilde_ULL_111 + 0.5*dDchi_L_0*gamtilde_LL_11*invgamtilde_UU_01/chi + 0.5*dDchi_L_1*gamtilde_LL_11*invgamtilde_UU_11/chi - 1.0*dDchi_L_1/chi + 0.5*dDchi_L_2*gamtilde_LL_11*invgamtilde_UU_12/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_112 = Christoffelgamtilde_ULL_112 + 0.5*dDchi_L_0*gamtilde_LL_12*invgamtilde_UU_01/chi + 0.5*dDchi_L_1*gamtilde_LL_12*invgamtilde_UU_11/chi + 0.5*dDchi_L_2*gamtilde_LL_12*invgamtilde_UU_12/chi - 0.5*dDchi_L_2/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_122 = Christoffelgamtilde_ULL_122 + 0.5*dDchi_L_0*gamtilde_LL_22*invgamtilde_UU_01/chi + 0.5*dDchi_L_1*gamtilde_LL_22*invgamtilde_UU_11/chi + 0.5*dDchi_L_2*gamtilde_LL_22*invgamtilde_UU_12/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_200 = Christoffelgamtilde_ULL_200 + 0.5*dDchi_L_0*gamtilde_LL_00*invgamtilde_UU_02/chi + 0.5*dDchi_L_1*gamtilde_LL_00*invgamtilde_UU_12/chi + 0.5*dDchi_L_2*gamtilde_LL_00*invgamtilde_UU_22/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_201 = Christoffelgamtilde_ULL_201 + 0.5*dDchi_L_0*gamtilde_LL_01*invgamtilde_UU_02/chi + 0.5*dDchi_L_1*gamtilde_LL_01*invgamtilde_UU_12/chi + 0.5*dDchi_L_2*gamtilde_LL_01*invgamtilde_UU_22/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_202 = Christoffelgamtilde_ULL_202 + 0.5*dDchi_L_0*gamtilde_LL_02*invgamtilde_UU_02/chi - 0.5*dDchi_L_0/chi + 0.5*dDchi_L_1*gamtilde_LL_02*invgamtilde_UU_12/chi + 0.5*dDchi_L_2*gamtilde_LL_02*invgamtilde_UU_22/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_211 = Christoffelgamtilde_ULL_211 + 0.5*dDchi_L_0*gamtilde_LL_11*invgamtilde_UU_02/chi + 0.5*dDchi_L_1*gamtilde_LL_11*invgamtilde_UU_12/chi + 0.5*dDchi_L_2*gamtilde_LL_11*invgamtilde_UU_22/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_212 = Christoffelgamtilde_ULL_212 + 0.5*dDchi_L_0*gamtilde_LL_12*invgamtilde_UU_02/chi + 0.5*dDchi_L_1*gamtilde_LL_12*invgamtilde_UU_12/chi - 0.5*dDchi_L_1/chi + 0.5*dDchi_L_2*gamtilde_LL_12*invgamtilde_UU_22/chi;
        amrex::Real ConformalChristoffelchigamtilde_ULL_222 = Christoffelgamtilde_ULL_222 + 0.5*dDchi_L_0*gamtilde_LL_22*invgamtilde_UU_02/chi + 0.5*dDchi_L_1*gamtilde_LL_22*invgamtilde_UU_12/chi + 0.5*dDchi_L_2*gamtilde_LL_22*invgamtilde_UU_22/chi - 1.0*dDchi_L_2/chi;

        amrex::Real ConformalChristoffelchigamtilde_LLL_000 = ConformalChristoffelchigamtilde_ULL_000*gamtilde_LL_00/chi + ConformalChristoffelchigamtilde_ULL_100*gamtilde_LL_01/chi + ConformalChristoffelchigamtilde_ULL_200*gamtilde_LL_02/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_001 = ConformalChristoffelchigamtilde_ULL_001*gamtilde_LL_00/chi + ConformalChristoffelchigamtilde_ULL_101*gamtilde_LL_01/chi + ConformalChristoffelchigamtilde_ULL_201*gamtilde_LL_02/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_002 = ConformalChristoffelchigamtilde_ULL_002*gamtilde_LL_00/chi + ConformalChristoffelchigamtilde_ULL_102*gamtilde_LL_01/chi + ConformalChristoffelchigamtilde_ULL_202*gamtilde_LL_02/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_011 = ConformalChristoffelchigamtilde_ULL_011*gamtilde_LL_00/chi + ConformalChristoffelchigamtilde_ULL_111*gamtilde_LL_01/chi + ConformalChristoffelchigamtilde_ULL_211*gamtilde_LL_02/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_012 = ConformalChristoffelchigamtilde_ULL_012*gamtilde_LL_00/chi + ConformalChristoffelchigamtilde_ULL_112*gamtilde_LL_01/chi + ConformalChristoffelchigamtilde_ULL_212*gamtilde_LL_02/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_022 = ConformalChristoffelchigamtilde_ULL_022*gamtilde_LL_00/chi + ConformalChristoffelchigamtilde_ULL_122*gamtilde_LL_01/chi + ConformalChristoffelchigamtilde_ULL_222*gamtilde_LL_02/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_100 = ConformalChristoffelchigamtilde_ULL_000*gamtilde_LL_01/chi + ConformalChristoffelchigamtilde_ULL_100*gamtilde_LL_11/chi + ConformalChristoffelchigamtilde_ULL_200*gamtilde_LL_12/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_101 = ConformalChristoffelchigamtilde_ULL_001*gamtilde_LL_01/chi + ConformalChristoffelchigamtilde_ULL_101*gamtilde_LL_11/chi + ConformalChristoffelchigamtilde_ULL_201*gamtilde_LL_12/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_102 = ConformalChristoffelchigamtilde_ULL_002*gamtilde_LL_01/chi + ConformalChristoffelchigamtilde_ULL_102*gamtilde_LL_11/chi + ConformalChristoffelchigamtilde_ULL_202*gamtilde_LL_12/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_111 = ConformalChristoffelchigamtilde_ULL_011*gamtilde_LL_01/chi + ConformalChristoffelchigamtilde_ULL_111*gamtilde_LL_11/chi + ConformalChristoffelchigamtilde_ULL_211*gamtilde_LL_12/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_112 = ConformalChristoffelchigamtilde_ULL_012*gamtilde_LL_01/chi + ConformalChristoffelchigamtilde_ULL_112*gamtilde_LL_11/chi + ConformalChristoffelchigamtilde_ULL_212*gamtilde_LL_12/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_122 = ConformalChristoffelchigamtilde_ULL_022*gamtilde_LL_01/chi + ConformalChristoffelchigamtilde_ULL_122*gamtilde_LL_11/chi + ConformalChristoffelchigamtilde_ULL_222*gamtilde_LL_12/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_200 = ConformalChristoffelchigamtilde_ULL_000*gamtilde_LL_02/chi + ConformalChristoffelchigamtilde_ULL_100*gamtilde_LL_12/chi + ConformalChristoffelchigamtilde_ULL_200*gamtilde_LL_22/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_201 = ConformalChristoffelchigamtilde_ULL_001*gamtilde_LL_02/chi + ConformalChristoffelchigamtilde_ULL_101*gamtilde_LL_12/chi + ConformalChristoffelchigamtilde_ULL_201*gamtilde_LL_22/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_202 = ConformalChristoffelchigamtilde_ULL_002*gamtilde_LL_02/chi + ConformalChristoffelchigamtilde_ULL_102*gamtilde_LL_12/chi + ConformalChristoffelchigamtilde_ULL_202*gamtilde_LL_22/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_211 = ConformalChristoffelchigamtilde_ULL_011*gamtilde_LL_02/chi + ConformalChristoffelchigamtilde_ULL_111*gamtilde_LL_12/chi + ConformalChristoffelchigamtilde_ULL_211*gamtilde_LL_22/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_212 = ConformalChristoffelchigamtilde_ULL_012*gamtilde_LL_02/chi + ConformalChristoffelchigamtilde_ULL_112*gamtilde_LL_12/chi + ConformalChristoffelchigamtilde_ULL_212*gamtilde_LL_22/chi;
        amrex::Real ConformalChristoffelchigamtilde_LLL_222 = ConformalChristoffelchigamtilde_ULL_022*gamtilde_LL_02/chi + ConformalChristoffelchigamtilde_ULL_122*gamtilde_LL_12/chi + ConformalChristoffelchigamtilde_ULL_222*gamtilde_LL_22/chi;

        amrex::Real dDDgam_LLLL_0000 = dDDgamtilde_LLLL_0000/chi + dDDchi_LL_00*gamtilde_LL_00/std::pow(chi, 2) - 2*dDchi_L_0*dDgamtilde_LLL_000/std::pow(chi, 2) + 2*std::pow(dDchi_L_0, 2)*gamtilde_LL_00/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0001 = dDDgamtilde_LLLL_0001/chi + dDDchi_LL_01*gamtilde_LL_00/std::pow(chi, 2) - dDchi_L_0*dDgamtilde_LLL_001/std::pow(chi, 2) - dDchi_L_1*dDgamtilde_LLL_000/std::pow(chi, 2) + 2*dDchi_L_0*dDchi_L_1*gamtilde_LL_00/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0002 = dDDgamtilde_LLLL_0002/chi + dDDchi_LL_02*gamtilde_LL_00/std::pow(chi, 2) - dDchi_L_0*dDgamtilde_LLL_002/std::pow(chi, 2) - dDchi_L_2*dDgamtilde_LLL_000/std::pow(chi, 2) + 2*dDchi_L_0*dDchi_L_2*gamtilde_LL_00/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0011 = dDDgamtilde_LLLL_0011/chi + dDDchi_LL_11*gamtilde_LL_00/std::pow(chi, 2) - 2*dDchi_L_1*dDgamtilde_LLL_001/std::pow(chi, 2) + 2*std::pow(dDchi_L_1, 2)*gamtilde_LL_00/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0012 = dDDgamtilde_LLLL_0012/chi + dDDchi_LL_12*gamtilde_LL_00/std::pow(chi, 2) - dDchi_L_1*dDgamtilde_LLL_002/std::pow(chi, 2) - dDchi_L_2*dDgamtilde_LLL_001/std::pow(chi, 2) + 2*dDchi_L_1*dDchi_L_2*gamtilde_LL_00/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0022 = dDDgamtilde_LLLL_0022/chi + dDDchi_LL_22*gamtilde_LL_00/std::pow(chi, 2) - 2*dDchi_L_2*dDgamtilde_LLL_002/std::pow(chi, 2) + 2*std::pow(dDchi_L_2, 2)*gamtilde_LL_00/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0100 = dDDgamtilde_LLLL_0100/chi + dDDchi_LL_00*gamtilde_LL_01/std::pow(chi, 2) - 2*dDchi_L_0*dDgamtilde_LLL_010/std::pow(chi, 2) + 2*std::pow(dDchi_L_0, 2)*gamtilde_LL_01/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0101 = dDDgamtilde_LLLL_0101/chi + dDDchi_LL_01*gamtilde_LL_01/std::pow(chi, 2) - dDchi_L_0*dDgamtilde_LLL_011/std::pow(chi, 2) - dDchi_L_1*dDgamtilde_LLL_010/std::pow(chi, 2) + 2*dDchi_L_0*dDchi_L_1*gamtilde_LL_01/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0102 = dDDgamtilde_LLLL_0102/chi + dDDchi_LL_02*gamtilde_LL_01/std::pow(chi, 2) - dDchi_L_0*dDgamtilde_LLL_012/std::pow(chi, 2) - dDchi_L_2*dDgamtilde_LLL_010/std::pow(chi, 2) + 2*dDchi_L_0*dDchi_L_2*gamtilde_LL_01/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0111 = dDDgamtilde_LLLL_0111/chi + dDDchi_LL_11*gamtilde_LL_01/std::pow(chi, 2) - 2*dDchi_L_1*dDgamtilde_LLL_011/std::pow(chi, 2) + 2*std::pow(dDchi_L_1, 2)*gamtilde_LL_01/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0112 = dDDgamtilde_LLLL_0112/chi + dDDchi_LL_12*gamtilde_LL_01/std::pow(chi, 2) - dDchi_L_1*dDgamtilde_LLL_012/std::pow(chi, 2) - dDchi_L_2*dDgamtilde_LLL_011/std::pow(chi, 2) + 2*dDchi_L_1*dDchi_L_2*gamtilde_LL_01/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0122 = dDDgamtilde_LLLL_0122/chi + dDDchi_LL_22*gamtilde_LL_01/std::pow(chi, 2) - 2*dDchi_L_2*dDgamtilde_LLL_012/std::pow(chi, 2) + 2*std::pow(dDchi_L_2, 2)*gamtilde_LL_01/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0200 = dDDgamtilde_LLLL_0200/chi + dDDchi_LL_00*gamtilde_LL_02/std::pow(chi, 2) - 2*dDchi_L_0*dDgamtilde_LLL_020/std::pow(chi, 2) + 2*std::pow(dDchi_L_0, 2)*gamtilde_LL_02/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0201 = dDDgamtilde_LLLL_0201/chi + dDDchi_LL_01*gamtilde_LL_02/std::pow(chi, 2) - dDchi_L_0*dDgamtilde_LLL_021/std::pow(chi, 2) - dDchi_L_1*dDgamtilde_LLL_020/std::pow(chi, 2) + 2*dDchi_L_0*dDchi_L_1*gamtilde_LL_02/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0202 = dDDgamtilde_LLLL_0202/chi + dDDchi_LL_02*gamtilde_LL_02/std::pow(chi, 2) - dDchi_L_0*dDgamtilde_LLL_022/std::pow(chi, 2) - dDchi_L_2*dDgamtilde_LLL_020/std::pow(chi, 2) + 2*dDchi_L_0*dDchi_L_2*gamtilde_LL_02/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0211 = dDDgamtilde_LLLL_0211/chi + dDDchi_LL_11*gamtilde_LL_02/std::pow(chi, 2) - 2*dDchi_L_1*dDgamtilde_LLL_021/std::pow(chi, 2) + 2*std::pow(dDchi_L_1, 2)*gamtilde_LL_02/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0212 = dDDgamtilde_LLLL_0212/chi + dDDchi_LL_12*gamtilde_LL_02/std::pow(chi, 2) - dDchi_L_1*dDgamtilde_LLL_022/std::pow(chi, 2) - dDchi_L_2*dDgamtilde_LLL_021/std::pow(chi, 2) + 2*dDchi_L_1*dDchi_L_2*gamtilde_LL_02/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_0222 = dDDgamtilde_LLLL_0222/chi + dDDchi_LL_22*gamtilde_LL_02/std::pow(chi, 2) - 2*dDchi_L_2*dDgamtilde_LLL_022/std::pow(chi, 2) + 2*std::pow(dDchi_L_2, 2)*gamtilde_LL_02/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_1100 = dDDgamtilde_LLLL_1100/chi + dDDchi_LL_00*gamtilde_LL_11/std::pow(chi, 2) - 2*dDchi_L_0*dDgamtilde_LLL_110/std::pow(chi, 2) + 2*std::pow(dDchi_L_0, 2)*gamtilde_LL_11/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_1101 = dDDgamtilde_LLLL_1101/chi + dDDchi_LL_01*gamtilde_LL_11/std::pow(chi, 2) - dDchi_L_0*dDgamtilde_LLL_111/std::pow(chi, 2) - dDchi_L_1*dDgamtilde_LLL_110/std::pow(chi, 2) + 2*dDchi_L_0*dDchi_L_1*gamtilde_LL_11/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_1102 = dDDgamtilde_LLLL_1102/chi + dDDchi_LL_02*gamtilde_LL_11/std::pow(chi, 2) - dDchi_L_0*dDgamtilde_LLL_112/std::pow(chi, 2) - dDchi_L_2*dDgamtilde_LLL_110/std::pow(chi, 2) + 2*dDchi_L_0*dDchi_L_2*gamtilde_LL_11/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_1111 = dDDgamtilde_LLLL_1111/chi + dDDchi_LL_11*gamtilde_LL_11/std::pow(chi, 2) - 2*dDchi_L_1*dDgamtilde_LLL_111/std::pow(chi, 2) + 2*std::pow(dDchi_L_1, 2)*gamtilde_LL_11/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_1112 = dDDgamtilde_LLLL_1112/chi + dDDchi_LL_12*gamtilde_LL_11/std::pow(chi, 2) - dDchi_L_1*dDgamtilde_LLL_112/std::pow(chi, 2) - dDchi_L_2*dDgamtilde_LLL_111/std::pow(chi, 2) + 2*dDchi_L_1*dDchi_L_2*gamtilde_LL_11/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_1122 = dDDgamtilde_LLLL_1122/chi + dDDchi_LL_22*gamtilde_LL_11/std::pow(chi, 2) - 2*dDchi_L_2*dDgamtilde_LLL_112/std::pow(chi, 2) + 2*std::pow(dDchi_L_2, 2)*gamtilde_LL_11/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_1200 = dDDgamtilde_LLLL_1200/chi + dDDchi_LL_00*gamtilde_LL_12/std::pow(chi, 2) - 2*dDchi_L_0*dDgamtilde_LLL_120/std::pow(chi, 2) + 2*std::pow(dDchi_L_0, 2)*gamtilde_LL_12/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_1201 = dDDgamtilde_LLLL_1201/chi + dDDchi_LL_01*gamtilde_LL_12/std::pow(chi, 2) - dDchi_L_0*dDgamtilde_LLL_121/std::pow(chi, 2) - dDchi_L_1*dDgamtilde_LLL_120/std::pow(chi, 2) + 2*dDchi_L_0*dDchi_L_1*gamtilde_LL_12/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_1202 = dDDgamtilde_LLLL_1202/chi + dDDchi_LL_02*gamtilde_LL_12/std::pow(chi, 2) - dDchi_L_0*dDgamtilde_LLL_122/std::pow(chi, 2) - dDchi_L_2*dDgamtilde_LLL_120/std::pow(chi, 2) + 2*dDchi_L_0*dDchi_L_2*gamtilde_LL_12/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_1211 = dDDgamtilde_LLLL_1211/chi + dDDchi_LL_11*gamtilde_LL_12/std::pow(chi, 2) - 2*dDchi_L_1*dDgamtilde_LLL_121/std::pow(chi, 2) + 2*std::pow(dDchi_L_1, 2)*gamtilde_LL_12/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_1212 = dDDgamtilde_LLLL_1212/chi + dDDchi_LL_12*gamtilde_LL_12/std::pow(chi, 2) - dDchi_L_1*dDgamtilde_LLL_122/std::pow(chi, 2) - dDchi_L_2*dDgamtilde_LLL_121/std::pow(chi, 2) + 2*dDchi_L_1*dDchi_L_2*gamtilde_LL_12/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_1222 = dDDgamtilde_LLLL_1222/chi + dDDchi_LL_22*gamtilde_LL_12/std::pow(chi, 2) - 2*dDchi_L_2*dDgamtilde_LLL_122/std::pow(chi, 2) + 2*std::pow(dDchi_L_2, 2)*gamtilde_LL_12/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_2200 = dDDgamtilde_LLLL_2200/chi + dDDchi_LL_00*gamtilde_LL_22/std::pow(chi, 2) - 2*dDchi_L_0*dDgamtilde_LLL_220/std::pow(chi, 2) + 2*std::pow(dDchi_L_0, 2)*gamtilde_LL_22/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_2201 = dDDgamtilde_LLLL_2201/chi + dDDchi_LL_01*gamtilde_LL_22/std::pow(chi, 2) - dDchi_L_0*dDgamtilde_LLL_221/std::pow(chi, 2) - dDchi_L_1*dDgamtilde_LLL_220/std::pow(chi, 2) + 2*dDchi_L_0*dDchi_L_1*gamtilde_LL_22/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_2202 = dDDgamtilde_LLLL_2202/chi + dDDchi_LL_02*gamtilde_LL_22/std::pow(chi, 2) - dDchi_L_0*dDgamtilde_LLL_222/std::pow(chi, 2) - dDchi_L_2*dDgamtilde_LLL_220/std::pow(chi, 2) + 2*dDchi_L_0*dDchi_L_2*gamtilde_LL_22/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_2211 = dDDgamtilde_LLLL_2211/chi + dDDchi_LL_11*gamtilde_LL_22/std::pow(chi, 2) - 2*dDchi_L_1*dDgamtilde_LLL_221/std::pow(chi, 2) + 2*std::pow(dDchi_L_1, 2)*gamtilde_LL_22/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_2212 = dDDgamtilde_LLLL_2212/chi + dDDchi_LL_12*gamtilde_LL_22/std::pow(chi, 2) - dDchi_L_1*dDgamtilde_LLL_222/std::pow(chi, 2) - dDchi_L_2*dDgamtilde_LLL_221/std::pow(chi, 2) + 2*dDchi_L_1*dDchi_L_2*gamtilde_LL_22/std::pow(chi, 3);
        amrex::Real dDDgam_LLLL_2222 = dDDgamtilde_LLLL_2222/chi + dDDchi_LL_22*gamtilde_LL_22/std::pow(chi, 2) - 2*dDchi_L_2*dDgamtilde_LLL_222/std::pow(chi, 2) + 2*std::pow(dDchi_L_2, 2)*gamtilde_LL_22/std::pow(chi, 3);

        amrex::Real Riem_LLLL_0101 = ConformalChristoffelchigamtilde_LLL_001*ConformalChristoffelchigamtilde_ULL_001 - ConformalChristoffelchigamtilde_LLL_011*ConformalChristoffelchigamtilde_ULL_000 + ConformalChristoffelchigamtilde_LLL_101*ConformalChristoffelchigamtilde_ULL_101 - ConformalChristoffelchigamtilde_LLL_111*ConformalChristoffelchigamtilde_ULL_100 + ConformalChristoffelchigamtilde_LLL_201*ConformalChristoffelchigamtilde_ULL_201 - ConformalChristoffelchigamtilde_LLL_211*ConformalChristoffelchigamtilde_ULL_200 - 0.5*dDDgam_LLLL_0011 + 1.0*dDDgam_LLLL_0101 - 0.5*dDDgam_LLLL_1100;
        amrex::Real Riem_LLLL_0102 = ConformalChristoffelchigamtilde_LLL_001*ConformalChristoffelchigamtilde_ULL_002 - ConformalChristoffelchigamtilde_LLL_012*ConformalChristoffelchigamtilde_ULL_000 + ConformalChristoffelchigamtilde_LLL_101*ConformalChristoffelchigamtilde_ULL_102 - ConformalChristoffelchigamtilde_LLL_112*ConformalChristoffelchigamtilde_ULL_100 + ConformalChristoffelchigamtilde_LLL_201*ConformalChristoffelchigamtilde_ULL_202 - ConformalChristoffelchigamtilde_LLL_212*ConformalChristoffelchigamtilde_ULL_200 - 0.5*dDDgam_LLLL_0012 + 0.5*dDDgam_LLLL_0102 + 0.5*dDDgam_LLLL_0201 - 0.5*dDDgam_LLLL_1200;
        amrex::Real Riem_LLLL_0112 = ConformalChristoffelchigamtilde_LLL_011*ConformalChristoffelchigamtilde_ULL_002 - ConformalChristoffelchigamtilde_LLL_012*ConformalChristoffelchigamtilde_ULL_001 + ConformalChristoffelchigamtilde_LLL_111*ConformalChristoffelchigamtilde_ULL_102 - ConformalChristoffelchigamtilde_LLL_112*ConformalChristoffelchigamtilde_ULL_101 + ConformalChristoffelchigamtilde_LLL_211*ConformalChristoffelchigamtilde_ULL_202 - ConformalChristoffelchigamtilde_LLL_212*ConformalChristoffelchigamtilde_ULL_201 - 0.5*dDDgam_LLLL_0112 + 0.5*dDDgam_LLLL_0211 + 0.5*dDDgam_LLLL_1102 - 0.5*dDDgam_LLLL_1201;
        amrex::Real Riem_LLLL_0201 = ConformalChristoffelchigamtilde_LLL_002*ConformalChristoffelchigamtilde_ULL_001 - ConformalChristoffelchigamtilde_LLL_012*ConformalChristoffelchigamtilde_ULL_000 + ConformalChristoffelchigamtilde_LLL_102*ConformalChristoffelchigamtilde_ULL_101 - ConformalChristoffelchigamtilde_LLL_112*ConformalChristoffelchigamtilde_ULL_100 + ConformalChristoffelchigamtilde_LLL_202*ConformalChristoffelchigamtilde_ULL_201 - ConformalChristoffelchigamtilde_LLL_212*ConformalChristoffelchigamtilde_ULL_200 - 0.5*dDDgam_LLLL_0012 + 0.5*dDDgam_LLLL_0102 + 0.5*dDDgam_LLLL_0201 - 0.5*dDDgam_LLLL_1200;
        amrex::Real Riem_LLLL_0202 = ConformalChristoffelchigamtilde_LLL_002*ConformalChristoffelchigamtilde_ULL_002 - ConformalChristoffelchigamtilde_LLL_022*ConformalChristoffelchigamtilde_ULL_000 + ConformalChristoffelchigamtilde_LLL_102*ConformalChristoffelchigamtilde_ULL_102 - ConformalChristoffelchigamtilde_LLL_122*ConformalChristoffelchigamtilde_ULL_100 + ConformalChristoffelchigamtilde_LLL_202*ConformalChristoffelchigamtilde_ULL_202 - ConformalChristoffelchigamtilde_LLL_222*ConformalChristoffelchigamtilde_ULL_200 - 0.5*dDDgam_LLLL_0022 + 1.0*dDDgam_LLLL_0202 - 0.5*dDDgam_LLLL_2200;
        amrex::Real Riem_LLLL_0212 = ConformalChristoffelchigamtilde_LLL_012*ConformalChristoffelchigamtilde_ULL_002 - ConformalChristoffelchigamtilde_LLL_022*ConformalChristoffelchigamtilde_ULL_001 + ConformalChristoffelchigamtilde_LLL_112*ConformalChristoffelchigamtilde_ULL_102 - ConformalChristoffelchigamtilde_LLL_122*ConformalChristoffelchigamtilde_ULL_101 + ConformalChristoffelchigamtilde_LLL_212*ConformalChristoffelchigamtilde_ULL_202 - ConformalChristoffelchigamtilde_LLL_222*ConformalChristoffelchigamtilde_ULL_201 - 0.5*dDDgam_LLLL_0122 + 0.5*dDDgam_LLLL_0212 + 0.5*dDDgam_LLLL_1202 - 0.5*dDDgam_LLLL_2201;
        amrex::Real Riem_LLLL_1201 = ConformalChristoffelchigamtilde_LLL_002*ConformalChristoffelchigamtilde_ULL_011 - ConformalChristoffelchigamtilde_LLL_012*ConformalChristoffelchigamtilde_ULL_001 + ConformalChristoffelchigamtilde_LLL_102*ConformalChristoffelchigamtilde_ULL_111 - ConformalChristoffelchigamtilde_LLL_112*ConformalChristoffelchigamtilde_ULL_101 + ConformalChristoffelchigamtilde_LLL_202*ConformalChristoffelchigamtilde_ULL_211 - ConformalChristoffelchigamtilde_LLL_212*ConformalChristoffelchigamtilde_ULL_201 - 0.5*dDDgam_LLLL_0112 + 0.5*dDDgam_LLLL_0211 + 0.5*dDDgam_LLLL_1102 - 0.5*dDDgam_LLLL_1201;
        amrex::Real Riem_LLLL_1202 = ConformalChristoffelchigamtilde_LLL_002*ConformalChristoffelchigamtilde_ULL_012 - ConformalChristoffelchigamtilde_LLL_022*ConformalChristoffelchigamtilde_ULL_001 + ConformalChristoffelchigamtilde_LLL_102*ConformalChristoffelchigamtilde_ULL_112 - ConformalChristoffelchigamtilde_LLL_122*ConformalChristoffelchigamtilde_ULL_101 + ConformalChristoffelchigamtilde_LLL_202*ConformalChristoffelchigamtilde_ULL_212 - ConformalChristoffelchigamtilde_LLL_222*ConformalChristoffelchigamtilde_ULL_201 - 0.5*dDDgam_LLLL_0122 + 0.5*dDDgam_LLLL_0212 + 0.5*dDDgam_LLLL_1202 - 0.5*dDDgam_LLLL_2201;
        amrex::Real Riem_LLLL_1212 = ConformalChristoffelchigamtilde_LLL_012*ConformalChristoffelchigamtilde_ULL_012 - ConformalChristoffelchigamtilde_LLL_022*ConformalChristoffelchigamtilde_ULL_011 + ConformalChristoffelchigamtilde_LLL_112*ConformalChristoffelchigamtilde_ULL_112 - ConformalChristoffelchigamtilde_LLL_122*ConformalChristoffelchigamtilde_ULL_111 + ConformalChristoffelchigamtilde_LLL_212*ConformalChristoffelchigamtilde_ULL_212 - ConformalChristoffelchigamtilde_LLL_222*ConformalChristoffelchigamtilde_ULL_211 - 0.5*dDDgam_LLLL_1122 + 1.0*dDDgam_LLLL_1212 - 0.5*dDDgam_LLLL_2211;

        amrex::Real Ricci_LL_00 = 3*Riem_LLLL_0101*invgam_UU_11 + 3*Riem_LLLL_0102*invgam_UU_12 + 3*Riem_LLLL_0201*invgam_UU_21 + 3*Riem_LLLL_0202*invgam_UU_22;
        amrex::Real Ricci_LL_01 = -3*Riem_LLLL_0101*invgam_UU_10 + 3*Riem_LLLL_0112*invgam_UU_12 - 3*Riem_LLLL_0201*invgam_UU_20 + 3*Riem_LLLL_0212*invgam_UU_22;
        amrex::Real Ricci_LL_02 = -3*Riem_LLLL_0102*invgam_UU_10 - 3*Riem_LLLL_0112*invgam_UU_11 - 3*Riem_LLLL_0202*invgam_UU_20 - 3*Riem_LLLL_0212*invgam_UU_21;
        amrex::Real Ricci_LL_11 = 3*Riem_LLLL_0101*invgam_UU_00 - 3*Riem_LLLL_0112*invgam_UU_02 - 3*Riem_LLLL_1201*invgam_UU_20 + 3*Riem_LLLL_1212*invgam_UU_22;
        amrex::Real Ricci_LL_12 = 3*Riem_LLLL_0102*invgam_UU_00 + 3*Riem_LLLL_0112*invgam_UU_01 - 3*Riem_LLLL_1202*invgam_UU_20 - 3*Riem_LLLL_1212*invgam_UU_21;
        amrex::Real Ricci_LL_22 = 3*Riem_LLLL_0202*invgam_UU_00 + 3*Riem_LLLL_0212*invgam_UU_01 + 3*Riem_LLLL_1202*invgam_UU_10 + 3*Riem_LLLL_1212*invgam_UU_11;

        amrex::Real K_LL_00 = (Atilde_LL_00 + 0.33333333333333331*gamtilde_LL_00*(Khat + 2*theta))/chi;
        amrex::Real K_LL_01 = (Atilde_LL_01 + 0.33333333333333331*gamtilde_LL_01*(Khat + 2*theta))/chi;
        amrex::Real K_LL_02 = (Atilde_LL_02 + 0.33333333333333331*gamtilde_LL_02*(Khat + 2*theta))/chi;
        amrex::Real K_LL_11 = (Atilde_LL_11 + 0.33333333333333331*gamtilde_LL_11*(Khat + 2*theta))/chi;
        amrex::Real K_LL_12 = (Atilde_LL_12 + 0.33333333333333331*gamtilde_LL_12*(Khat + 2*theta))/chi;
        amrex::Real K_LL_22 = (Atilde_LL_22 + 0.33333333333333331*gamtilde_LL_22*(Khat + 2*theta))/chi;

        amrex::Real Ksclr = Khat + 2*theta;

        /*amrex::Real detgam = gam_LL_00*gam_LL_11*gam_LL_22 - gam_LL_00*std::pow(gam_LL_12, 2) - std::pow(gam_LL_01, 2)*gam_LL_22 + 2*gam_LL_01*gam_LL_02*gam_LL_12 - std::pow(gam_LL_02, 2)*gam_LL_11;

        amrex::Real x = x0;

        amrex::Real y = x1;

        amrex::Real z = x2;

        amrex::Real rcoord = std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));

        amrex::Real v0_U_0 = -y;
        amrex::Real v0_U_1 = x;
        amrex::Real v0_U_2 = 0;

        amrex::Real v1_U_0 = x;
        amrex::Real v1_U_1 = y;
        amrex::Real v1_U_2 = z;

        amrex::Real v2_U_0 = std::sqrt(detgam)*invgam_UU_00*v0_U_1*v1_U_2 - std::sqrt(detgam)*invgam_UU_00*v0_U_2*v1_U_1 - std::sqrt(detgam)*invgam_UU_01*v0_U_0*v1_U_2 + std::sqrt(detgam)*invgam_UU_01*v0_U_2*v1_U_0 + std::sqrt(detgam)*invgam_UU_02*v0_U_0*v1_U_1 - std::sqrt(detgam)*invgam_UU_02*v0_U_1*v1_U_0;
        amrex::Real v2_U_1 = std::sqrt(detgam)*invgam_UU_10*v0_U_1*v1_U_2 - std::sqrt(detgam)*invgam_UU_10*v0_U_2*v1_U_1 - std::sqrt(detgam)*invgam_UU_11*v0_U_0*v1_U_2 + std::sqrt(detgam)*invgam_UU_11*v0_U_2*v1_U_0 + std::sqrt(detgam)*invgam_UU_12*v0_U_0*v1_U_1 - std::sqrt(detgam)*invgam_UU_12*v0_U_1*v1_U_0;
        amrex::Real v2_U_2 = std::sqrt(detgam)*invgam_UU_20*v0_U_1*v1_U_2 - std::sqrt(detgam)*invgam_UU_20*v0_U_2*v1_U_1 - std::sqrt(detgam)*invgam_UU_21*v0_U_0*v1_U_2 + std::sqrt(detgam)*invgam_UU_21*v0_U_2*v1_U_0 + std::sqrt(detgam)*invgam_UU_22*v0_U_0*v1_U_1 - std::sqrt(detgam)*invgam_UU_22*v0_U_1*v1_U_0;

        amrex::Real omega00 = gam_LL_00*std::pow(v0_U_0, 2) + gam_LL_01*v0_U_0*v0_U_1 + gam_LL_02*v0_U_0*v0_U_2 + gam_LL_10*v0_U_0*v0_U_1 + gam_LL_11*std::pow(v0_U_1, 2) + gam_LL_12*v0_U_1*v0_U_2 + gam_LL_20*v0_U_0*v0_U_2 + gam_LL_21*v0_U_1*v0_U_2 + gam_LL_22*std::pow(v0_U_2, 2);

        amrex::Real e0_U_0 = v0_U_0/omega00;
        amrex::Real e0_U_1 = v0_U_1/omega00;
        amrex::Real e0_U_2 = v0_U_2/omega00;

        amrex::Real omega01 = e0_U_0*gam_LL_00*v1_U_0 + e0_U_0*gam_LL_01*v1_U_1 + e0_U_0*gam_LL_02*v1_U_2 + e0_U_1*gam_LL_10*v1_U_0 + e0_U_1*gam_LL_11*v1_U_1 + e0_U_1*gam_LL_12*v1_U_2 + e0_U_2*gam_LL_20*v1_U_0 + e0_U_2*gam_LL_21*v1_U_1 + e0_U_2*gam_LL_22*v1_U_2;

        amrex::Real vp1_U_0 = -e0_U_0*omega01 + v1_U_0;
        amrex::Real vp1_U_1 = -e0_U_1*omega01 + v1_U_1;
        amrex::Real vp1_U_2 = -e0_U_2*omega01 + v1_U_2;

        amrex::Real omega11 = gam_LL_00*std::pow(vp1_U_0, 2) + gam_LL_01*vp1_U_0*vp1_U_1 + gam_LL_02*vp1_U_0*vp1_U_2 + gam_LL_10*vp1_U_0*vp1_U_1 + gam_LL_11*std::pow(vp1_U_1, 2) + gam_LL_12*vp1_U_1*vp1_U_2 + gam_LL_20*vp1_U_0*vp1_U_2 + gam_LL_21*vp1_U_1*vp1_U_2 + gam_LL_22*std::pow(vp1_U_2, 2);

        amrex::Real e1_U_0 = vp1_U_0/omega11;
        amrex::Real e1_U_1 = vp1_U_1/omega11;
        amrex::Real e1_U_2 = vp1_U_2/omega11;

        amrex::Real omega02 = e0_U_0*gam_LL_00*v2_U_0 + e0_U_0*gam_LL_01*v2_U_1 + e0_U_0*gam_LL_02*v2_U_2 + e0_U_1*gam_LL_10*v2_U_0 + e0_U_1*gam_LL_11*v2_U_1 + e0_U_1*gam_LL_12*v2_U_2 + e0_U_2*gam_LL_20*v2_U_0 + e0_U_2*gam_LL_21*v2_U_1 + e0_U_2*gam_LL_22*v2_U_2;

        amrex::Real omega12 = e1_U_0*gam_LL_00*v2_U_0 + e1_U_0*gam_LL_01*v2_U_1 + e1_U_0*gam_LL_02*v2_U_2 + e1_U_1*gam_LL_10*v2_U_0 + e1_U_1*gam_LL_11*v2_U_1 + e1_U_1*gam_LL_12*v2_U_2 + e1_U_2*gam_LL_20*v2_U_0 + e1_U_2*gam_LL_21*v2_U_1 + e1_U_2*gam_LL_22*v2_U_2;

        amrex::Real vp2_U_0 = -e0_U_0*omega02 - e1_U_0*omega12 + v2_U_0;
        amrex::Real vp2_U_1 = -e0_U_1*omega02 - e1_U_1*omega12 + v2_U_1;
        amrex::Real vp2_U_2 = -e0_U_2*omega02 - e1_U_2*omega12 + v2_U_2;

        amrex::Real omega22 = gam_LL_00*std::pow(vp2_U_0, 2) + gam_LL_01*vp2_U_0*vp2_U_1 + gam_LL_02*vp2_U_0*vp2_U_2 + gam_LL_10*vp2_U_0*vp2_U_1 + gam_LL_11*std::pow(vp2_U_1, 2) + gam_LL_12*vp2_U_1*vp2_U_2 + gam_LL_20*vp2_U_0*vp2_U_2 + gam_LL_21*vp2_U_1*vp2_U_2 + gam_LL_22*std::pow(vp2_U_2, 2) + (std::fabs(gam_LL_00*std::pow(vp2_U_0, 2)) < 1.0000000000000001e-5) + (std::fabs(gam_LL_11*std::pow(vp2_U_1, 2)) < 1.0000000000000001e-5) + (std::fabs(gam_LL_22*std::pow(vp2_U_2, 2)) < 1.0000000000000001e-5) + (std::fabs(gam_LL_01*vp2_U_0*vp2_U_1) < 1.0000000000000001e-5) + (std::fabs(gam_LL_02*vp2_U_0*vp2_U_2) < 1.0000000000000001e-5) + (std::fabs(gam_LL_10*vp2_U_0*vp2_U_1) < 1.0000000000000001e-5) + (std::fabs(gam_LL_12*vp2_U_1*vp2_U_2) < 1.0000000000000001e-5) + (std::fabs(gam_LL_20*vp2_U_0*vp2_U_2) < 1.0000000000000001e-5) + (std::fabs(gam_LL_21*vp2_U_1*vp2_U_2) < 1.0000000000000001e-5);

        amrex::Real e2_U_0 = vp2_U_0/omega22;
        amrex::Real e2_U_1 = vp2_U_1/omega22;
        amrex::Real e2_U_2 = vp2_U_2/omega22;

        amrex::Real u_U_0 = -beta_U_0/alpha;
        amrex::Real u_U_1 = -beta_U_1/alpha;
        amrex::Real u_U_2 = -beta_U_2/alpha;
        amrex::Real u_U_3 = 1.0/alpha;

        amrex::Real r_U_0 = e1_U_0;
        amrex::Real r_U_1 = e1_U_1;
        amrex::Real r_U_2 = e1_U_2;
        amrex::Real r_U_3 = 0;

        amrex::Real theta_U_0 = e2_U_0;
        amrex::Real theta_U_1 = e2_U_1;
        amrex::Real theta_U_2 = e2_U_2;
        amrex::Real theta_U_3 = 0;

        amrex::Real phi_U_0 = e0_U_0;
        amrex::Real phi_U_1 = e0_U_1;
        amrex::Real phi_U_2 = e0_U_2;
        amrex::Real phi_U_3 = 0;

        amrex::Real l_U_0 = (1.0/2.0)*M_SQRT2*(r_U_0 + u_U_0);
        amrex::Real l_U_1 = (1.0/2.0)*M_SQRT2*(r_U_1 + u_U_1);
        amrex::Real l_U_2 = (1.0/2.0)*M_SQRT2*(r_U_2 + u_U_2);
        amrex::Real l_U_3 = (1.0/2.0)*M_SQRT2*(r_U_3 + u_U_3);

        amrex::Real n_U_0 = (1.0/2.0)*M_SQRT2*(-r_U_0 + u_U_0);
        amrex::Real n_U_1 = (1.0/2.0)*M_SQRT2*(-r_U_1 + u_U_1);
        amrex::Real n_U_2 = (1.0/2.0)*M_SQRT2*(-r_U_2 + u_U_2);
        amrex::Real n_U_3 = (1.0/2.0)*M_SQRT2*(-r_U_3 + u_U_3);

        amrex::Real mre_U_0 = (1.0/2.0)*M_SQRT2*theta_U_0;
        amrex::Real mre_U_1 = (1.0/2.0)*M_SQRT2*theta_U_1;
        amrex::Real mre_U_2 = (1.0/2.0)*M_SQRT2*theta_U_2;
        amrex::Real mre_U_3 = (1.0/2.0)*M_SQRT2*theta_U_3;

        amrex::Real mim_U_0 = (1.0/2.0)*M_SQRT2*phi_U_0;
        amrex::Real mim_U_1 = (1.0/2.0)*M_SQRT2*phi_U_1;
        amrex::Real mim_U_2 = (1.0/2.0)*M_SQRT2*phi_U_2;
        amrex::Real mim_U_3 = (1.0/2.0)*M_SQRT2*phi_U_3;

        amrex::Real nmnmre_UUUU_0000 = -std::pow(mim_U_0, 2)*std::pow(n_U_0, 2) + std::pow(mre_U_0, 2)*std::pow(n_U_0, 2);
        amrex::Real nmnmre_UUUU_0001 = -mim_U_0*mim_U_1*std::pow(n_U_0, 2) + mre_U_0*mre_U_1*std::pow(n_U_0, 2);
        amrex::Real nmnmre_UUUU_0002 = -mim_U_0*mim_U_2*std::pow(n_U_0, 2) + mre_U_0*mre_U_2*std::pow(n_U_0, 2);
        amrex::Real nmnmre_UUUU_0003 = -mim_U_0*mim_U_3*std::pow(n_U_0, 2) + mre_U_0*mre_U_3*std::pow(n_U_0, 2);
        amrex::Real nmnmre_UUUU_0010 = -std::pow(mim_U_0, 2)*n_U_0*n_U_1 + std::pow(mre_U_0, 2)*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_0011 = -mim_U_0*mim_U_1*n_U_0*n_U_1 + mre_U_0*mre_U_1*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_0012 = -mim_U_0*mim_U_2*n_U_0*n_U_1 + mre_U_0*mre_U_2*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_0013 = -mim_U_0*mim_U_3*n_U_0*n_U_1 + mre_U_0*mre_U_3*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_0020 = -std::pow(mim_U_0, 2)*n_U_0*n_U_2 + std::pow(mre_U_0, 2)*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_0021 = -mim_U_0*mim_U_1*n_U_0*n_U_2 + mre_U_0*mre_U_1*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_0022 = -mim_U_0*mim_U_2*n_U_0*n_U_2 + mre_U_0*mre_U_2*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_0023 = -mim_U_0*mim_U_3*n_U_0*n_U_2 + mre_U_0*mre_U_3*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_0030 = -std::pow(mim_U_0, 2)*n_U_0*n_U_3 + std::pow(mre_U_0, 2)*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_0031 = -mim_U_0*mim_U_1*n_U_0*n_U_3 + mre_U_0*mre_U_1*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_0032 = -mim_U_0*mim_U_2*n_U_0*n_U_3 + mre_U_0*mre_U_2*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_0033 = -mim_U_0*mim_U_3*n_U_0*n_U_3 + mre_U_0*mre_U_3*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_0101 = -std::pow(mim_U_1, 2)*std::pow(n_U_0, 2) + std::pow(mre_U_1, 2)*std::pow(n_U_0, 2);
        amrex::Real nmnmre_UUUU_0102 = -mim_U_1*mim_U_2*std::pow(n_U_0, 2) + mre_U_1*mre_U_2*std::pow(n_U_0, 2);
        amrex::Real nmnmre_UUUU_0103 = -mim_U_1*mim_U_3*std::pow(n_U_0, 2) + mre_U_1*mre_U_3*std::pow(n_U_0, 2);
        amrex::Real nmnmre_UUUU_0111 = -std::pow(mim_U_1, 2)*n_U_0*n_U_1 + std::pow(mre_U_1, 2)*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_0112 = -mim_U_1*mim_U_2*n_U_0*n_U_1 + mre_U_1*mre_U_2*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_0113 = -mim_U_1*mim_U_3*n_U_0*n_U_1 + mre_U_1*mre_U_3*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_0121 = -std::pow(mim_U_1, 2)*n_U_0*n_U_2 + std::pow(mre_U_1, 2)*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_0122 = -mim_U_1*mim_U_2*n_U_0*n_U_2 + mre_U_1*mre_U_2*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_0123 = -mim_U_1*mim_U_3*n_U_0*n_U_2 + mre_U_1*mre_U_3*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_0131 = -std::pow(mim_U_1, 2)*n_U_0*n_U_3 + std::pow(mre_U_1, 2)*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_0132 = -mim_U_1*mim_U_2*n_U_0*n_U_3 + mre_U_1*mre_U_2*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_0133 = -mim_U_1*mim_U_3*n_U_0*n_U_3 + mre_U_1*mre_U_3*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_0202 = -std::pow(mim_U_2, 2)*std::pow(n_U_0, 2) + std::pow(mre_U_2, 2)*std::pow(n_U_0, 2);
        amrex::Real nmnmre_UUUU_0203 = -mim_U_2*mim_U_3*std::pow(n_U_0, 2) + mre_U_2*mre_U_3*std::pow(n_U_0, 2);
        amrex::Real nmnmre_UUUU_0212 = -std::pow(mim_U_2, 2)*n_U_0*n_U_1 + std::pow(mre_U_2, 2)*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_0213 = -mim_U_2*mim_U_3*n_U_0*n_U_1 + mre_U_2*mre_U_3*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_0222 = -std::pow(mim_U_2, 2)*n_U_0*n_U_2 + std::pow(mre_U_2, 2)*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_0223 = -mim_U_2*mim_U_3*n_U_0*n_U_2 + mre_U_2*mre_U_3*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_0232 = -std::pow(mim_U_2, 2)*n_U_0*n_U_3 + std::pow(mre_U_2, 2)*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_0233 = -mim_U_2*mim_U_3*n_U_0*n_U_3 + mre_U_2*mre_U_3*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_0303 = -std::pow(mim_U_3, 2)*std::pow(n_U_0, 2) + std::pow(mre_U_3, 2)*std::pow(n_U_0, 2);
        amrex::Real nmnmre_UUUU_0313 = -std::pow(mim_U_3, 2)*n_U_0*n_U_1 + std::pow(mre_U_3, 2)*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_0323 = -std::pow(mim_U_3, 2)*n_U_0*n_U_2 + std::pow(mre_U_3, 2)*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_0333 = -std::pow(mim_U_3, 2)*n_U_0*n_U_3 + std::pow(mre_U_3, 2)*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_1000 = -std::pow(mim_U_0, 2)*n_U_0*n_U_1 + std::pow(mre_U_0, 2)*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_1001 = -mim_U_0*mim_U_1*n_U_0*n_U_1 + mre_U_0*mre_U_1*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_1002 = -mim_U_0*mim_U_2*n_U_0*n_U_1 + mre_U_0*mre_U_2*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_1003 = -mim_U_0*mim_U_3*n_U_0*n_U_1 + mre_U_0*mre_U_3*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_1010 = -std::pow(mim_U_0, 2)*std::pow(n_U_1, 2) + std::pow(mre_U_0, 2)*std::pow(n_U_1, 2);
        amrex::Real nmnmre_UUUU_1011 = -mim_U_0*mim_U_1*std::pow(n_U_1, 2) + mre_U_0*mre_U_1*std::pow(n_U_1, 2);
        amrex::Real nmnmre_UUUU_1012 = -mim_U_0*mim_U_2*std::pow(n_U_1, 2) + mre_U_0*mre_U_2*std::pow(n_U_1, 2);
        amrex::Real nmnmre_UUUU_1013 = -mim_U_0*mim_U_3*std::pow(n_U_1, 2) + mre_U_0*mre_U_3*std::pow(n_U_1, 2);
        amrex::Real nmnmre_UUUU_1020 = -std::pow(mim_U_0, 2)*n_U_1*n_U_2 + std::pow(mre_U_0, 2)*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_1021 = -mim_U_0*mim_U_1*n_U_1*n_U_2 + mre_U_0*mre_U_1*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_1022 = -mim_U_0*mim_U_2*n_U_1*n_U_2 + mre_U_0*mre_U_2*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_1023 = -mim_U_0*mim_U_3*n_U_1*n_U_2 + mre_U_0*mre_U_3*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_1030 = -std::pow(mim_U_0, 2)*n_U_1*n_U_3 + std::pow(mre_U_0, 2)*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_1031 = -mim_U_0*mim_U_1*n_U_1*n_U_3 + mre_U_0*mre_U_1*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_1032 = -mim_U_0*mim_U_2*n_U_1*n_U_3 + mre_U_0*mre_U_2*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_1033 = -mim_U_0*mim_U_3*n_U_1*n_U_3 + mre_U_0*mre_U_3*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_1101 = -std::pow(mim_U_1, 2)*n_U_0*n_U_1 + std::pow(mre_U_1, 2)*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_1102 = -mim_U_1*mim_U_2*n_U_0*n_U_1 + mre_U_1*mre_U_2*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_1103 = -mim_U_1*mim_U_3*n_U_0*n_U_1 + mre_U_1*mre_U_3*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_1111 = -std::pow(mim_U_1, 2)*std::pow(n_U_1, 2) + std::pow(mre_U_1, 2)*std::pow(n_U_1, 2);
        amrex::Real nmnmre_UUUU_1112 = -mim_U_1*mim_U_2*std::pow(n_U_1, 2) + mre_U_1*mre_U_2*std::pow(n_U_1, 2);
        amrex::Real nmnmre_UUUU_1113 = -mim_U_1*mim_U_3*std::pow(n_U_1, 2) + mre_U_1*mre_U_3*std::pow(n_U_1, 2);
        amrex::Real nmnmre_UUUU_1121 = -std::pow(mim_U_1, 2)*n_U_1*n_U_2 + std::pow(mre_U_1, 2)*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_1122 = -mim_U_1*mim_U_2*n_U_1*n_U_2 + mre_U_1*mre_U_2*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_1123 = -mim_U_1*mim_U_3*n_U_1*n_U_2 + mre_U_1*mre_U_3*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_1131 = -std::pow(mim_U_1, 2)*n_U_1*n_U_3 + std::pow(mre_U_1, 2)*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_1132 = -mim_U_1*mim_U_2*n_U_1*n_U_3 + mre_U_1*mre_U_2*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_1133 = -mim_U_1*mim_U_3*n_U_1*n_U_3 + mre_U_1*mre_U_3*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_1202 = -std::pow(mim_U_2, 2)*n_U_0*n_U_1 + std::pow(mre_U_2, 2)*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_1203 = -mim_U_2*mim_U_3*n_U_0*n_U_1 + mre_U_2*mre_U_3*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_1212 = -std::pow(mim_U_2, 2)*std::pow(n_U_1, 2) + std::pow(mre_U_2, 2)*std::pow(n_U_1, 2);
        amrex::Real nmnmre_UUUU_1213 = -mim_U_2*mim_U_3*std::pow(n_U_1, 2) + mre_U_2*mre_U_3*std::pow(n_U_1, 2);
        amrex::Real nmnmre_UUUU_1222 = -std::pow(mim_U_2, 2)*n_U_1*n_U_2 + std::pow(mre_U_2, 2)*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_1223 = -mim_U_2*mim_U_3*n_U_1*n_U_2 + mre_U_2*mre_U_3*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_1232 = -std::pow(mim_U_2, 2)*n_U_1*n_U_3 + std::pow(mre_U_2, 2)*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_1233 = -mim_U_2*mim_U_3*n_U_1*n_U_3 + mre_U_2*mre_U_3*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_1303 = -std::pow(mim_U_3, 2)*n_U_0*n_U_1 + std::pow(mre_U_3, 2)*n_U_0*n_U_1;
        amrex::Real nmnmre_UUUU_1313 = -std::pow(mim_U_3, 2)*std::pow(n_U_1, 2) + std::pow(mre_U_3, 2)*std::pow(n_U_1, 2);
        amrex::Real nmnmre_UUUU_1323 = -std::pow(mim_U_3, 2)*n_U_1*n_U_2 + std::pow(mre_U_3, 2)*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_1333 = -std::pow(mim_U_3, 2)*n_U_1*n_U_3 + std::pow(mre_U_3, 2)*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_2000 = -std::pow(mim_U_0, 2)*n_U_0*n_U_2 + std::pow(mre_U_0, 2)*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_2001 = -mim_U_0*mim_U_1*n_U_0*n_U_2 + mre_U_0*mre_U_1*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_2002 = -mim_U_0*mim_U_2*n_U_0*n_U_2 + mre_U_0*mre_U_2*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_2003 = -mim_U_0*mim_U_3*n_U_0*n_U_2 + mre_U_0*mre_U_3*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_2010 = -std::pow(mim_U_0, 2)*n_U_1*n_U_2 + std::pow(mre_U_0, 2)*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_2011 = -mim_U_0*mim_U_1*n_U_1*n_U_2 + mre_U_0*mre_U_1*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_2012 = -mim_U_0*mim_U_2*n_U_1*n_U_2 + mre_U_0*mre_U_2*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_2013 = -mim_U_0*mim_U_3*n_U_1*n_U_2 + mre_U_0*mre_U_3*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_2020 = -std::pow(mim_U_0, 2)*std::pow(n_U_2, 2) + std::pow(mre_U_0, 2)*std::pow(n_U_2, 2);
        amrex::Real nmnmre_UUUU_2021 = -mim_U_0*mim_U_1*std::pow(n_U_2, 2) + mre_U_0*mre_U_1*std::pow(n_U_2, 2);
        amrex::Real nmnmre_UUUU_2022 = -mim_U_0*mim_U_2*std::pow(n_U_2, 2) + mre_U_0*mre_U_2*std::pow(n_U_2, 2);
        amrex::Real nmnmre_UUUU_2023 = -mim_U_0*mim_U_3*std::pow(n_U_2, 2) + mre_U_0*mre_U_3*std::pow(n_U_2, 2);
        amrex::Real nmnmre_UUUU_2030 = -std::pow(mim_U_0, 2)*n_U_2*n_U_3 + std::pow(mre_U_0, 2)*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_2031 = -mim_U_0*mim_U_1*n_U_2*n_U_3 + mre_U_0*mre_U_1*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_2032 = -mim_U_0*mim_U_2*n_U_2*n_U_3 + mre_U_0*mre_U_2*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_2033 = -mim_U_0*mim_U_3*n_U_2*n_U_3 + mre_U_0*mre_U_3*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_2101 = -std::pow(mim_U_1, 2)*n_U_0*n_U_2 + std::pow(mre_U_1, 2)*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_2102 = -mim_U_1*mim_U_2*n_U_0*n_U_2 + mre_U_1*mre_U_2*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_2103 = -mim_U_1*mim_U_3*n_U_0*n_U_2 + mre_U_1*mre_U_3*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_2111 = -std::pow(mim_U_1, 2)*n_U_1*n_U_2 + std::pow(mre_U_1, 2)*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_2112 = -mim_U_1*mim_U_2*n_U_1*n_U_2 + mre_U_1*mre_U_2*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_2113 = -mim_U_1*mim_U_3*n_U_1*n_U_2 + mre_U_1*mre_U_3*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_2121 = -std::pow(mim_U_1, 2)*std::pow(n_U_2, 2) + std::pow(mre_U_1, 2)*std::pow(n_U_2, 2);
        amrex::Real nmnmre_UUUU_2122 = -mim_U_1*mim_U_2*std::pow(n_U_2, 2) + mre_U_1*mre_U_2*std::pow(n_U_2, 2);
        amrex::Real nmnmre_UUUU_2123 = -mim_U_1*mim_U_3*std::pow(n_U_2, 2) + mre_U_1*mre_U_3*std::pow(n_U_2, 2);
        amrex::Real nmnmre_UUUU_2131 = -std::pow(mim_U_1, 2)*n_U_2*n_U_3 + std::pow(mre_U_1, 2)*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_2132 = -mim_U_1*mim_U_2*n_U_2*n_U_3 + mre_U_1*mre_U_2*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_2133 = -mim_U_1*mim_U_3*n_U_2*n_U_3 + mre_U_1*mre_U_3*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_2202 = -std::pow(mim_U_2, 2)*n_U_0*n_U_2 + std::pow(mre_U_2, 2)*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_2203 = -mim_U_2*mim_U_3*n_U_0*n_U_2 + mre_U_2*mre_U_3*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_2212 = -std::pow(mim_U_2, 2)*n_U_1*n_U_2 + std::pow(mre_U_2, 2)*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_2213 = -mim_U_2*mim_U_3*n_U_1*n_U_2 + mre_U_2*mre_U_3*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_2222 = -std::pow(mim_U_2, 2)*std::pow(n_U_2, 2) + std::pow(mre_U_2, 2)*std::pow(n_U_2, 2);
        amrex::Real nmnmre_UUUU_2223 = -mim_U_2*mim_U_3*std::pow(n_U_2, 2) + mre_U_2*mre_U_3*std::pow(n_U_2, 2);
        amrex::Real nmnmre_UUUU_2232 = -std::pow(mim_U_2, 2)*n_U_2*n_U_3 + std::pow(mre_U_2, 2)*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_2233 = -mim_U_2*mim_U_3*n_U_2*n_U_3 + mre_U_2*mre_U_3*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_2303 = -std::pow(mim_U_3, 2)*n_U_0*n_U_2 + std::pow(mre_U_3, 2)*n_U_0*n_U_2;
        amrex::Real nmnmre_UUUU_2313 = -std::pow(mim_U_3, 2)*n_U_1*n_U_2 + std::pow(mre_U_3, 2)*n_U_1*n_U_2;
        amrex::Real nmnmre_UUUU_2323 = -std::pow(mim_U_3, 2)*std::pow(n_U_2, 2) + std::pow(mre_U_3, 2)*std::pow(n_U_2, 2);
        amrex::Real nmnmre_UUUU_2333 = -std::pow(mim_U_3, 2)*n_U_2*n_U_3 + std::pow(mre_U_3, 2)*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_3000 = -std::pow(mim_U_0, 2)*n_U_0*n_U_3 + std::pow(mre_U_0, 2)*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_3001 = -mim_U_0*mim_U_1*n_U_0*n_U_3 + mre_U_0*mre_U_1*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_3002 = -mim_U_0*mim_U_2*n_U_0*n_U_3 + mre_U_0*mre_U_2*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_3003 = -mim_U_0*mim_U_3*n_U_0*n_U_3 + mre_U_0*mre_U_3*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_3010 = -std::pow(mim_U_0, 2)*n_U_1*n_U_3 + std::pow(mre_U_0, 2)*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_3011 = -mim_U_0*mim_U_1*n_U_1*n_U_3 + mre_U_0*mre_U_1*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_3012 = -mim_U_0*mim_U_2*n_U_1*n_U_3 + mre_U_0*mre_U_2*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_3013 = -mim_U_0*mim_U_3*n_U_1*n_U_3 + mre_U_0*mre_U_3*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_3020 = -std::pow(mim_U_0, 2)*n_U_2*n_U_3 + std::pow(mre_U_0, 2)*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_3021 = -mim_U_0*mim_U_1*n_U_2*n_U_3 + mre_U_0*mre_U_1*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_3022 = -mim_U_0*mim_U_2*n_U_2*n_U_3 + mre_U_0*mre_U_2*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_3023 = -mim_U_0*mim_U_3*n_U_2*n_U_3 + mre_U_0*mre_U_3*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_3030 = -std::pow(mim_U_0, 2)*std::pow(n_U_3, 2) + std::pow(mre_U_0, 2)*std::pow(n_U_3, 2);
        amrex::Real nmnmre_UUUU_3031 = -mim_U_0*mim_U_1*std::pow(n_U_3, 2) + mre_U_0*mre_U_1*std::pow(n_U_3, 2);
        amrex::Real nmnmre_UUUU_3032 = -mim_U_0*mim_U_2*std::pow(n_U_3, 2) + mre_U_0*mre_U_2*std::pow(n_U_3, 2);
        amrex::Real nmnmre_UUUU_3033 = -mim_U_0*mim_U_3*std::pow(n_U_3, 2) + mre_U_0*mre_U_3*std::pow(n_U_3, 2);
        amrex::Real nmnmre_UUUU_3101 = -std::pow(mim_U_1, 2)*n_U_0*n_U_3 + std::pow(mre_U_1, 2)*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_3102 = -mim_U_1*mim_U_2*n_U_0*n_U_3 + mre_U_1*mre_U_2*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_3103 = -mim_U_1*mim_U_3*n_U_0*n_U_3 + mre_U_1*mre_U_3*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_3111 = -std::pow(mim_U_1, 2)*n_U_1*n_U_3 + std::pow(mre_U_1, 2)*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_3112 = -mim_U_1*mim_U_2*n_U_1*n_U_3 + mre_U_1*mre_U_2*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_3113 = -mim_U_1*mim_U_3*n_U_1*n_U_3 + mre_U_1*mre_U_3*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_3121 = -std::pow(mim_U_1, 2)*n_U_2*n_U_3 + std::pow(mre_U_1, 2)*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_3122 = -mim_U_1*mim_U_2*n_U_2*n_U_3 + mre_U_1*mre_U_2*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_3123 = -mim_U_1*mim_U_3*n_U_2*n_U_3 + mre_U_1*mre_U_3*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_3131 = -std::pow(mim_U_1, 2)*std::pow(n_U_3, 2) + std::pow(mre_U_1, 2)*std::pow(n_U_3, 2);
        amrex::Real nmnmre_UUUU_3132 = -mim_U_1*mim_U_2*std::pow(n_U_3, 2) + mre_U_1*mre_U_2*std::pow(n_U_3, 2);
        amrex::Real nmnmre_UUUU_3133 = -mim_U_1*mim_U_3*std::pow(n_U_3, 2) + mre_U_1*mre_U_3*std::pow(n_U_3, 2);
        amrex::Real nmnmre_UUUU_3202 = -std::pow(mim_U_2, 2)*n_U_0*n_U_3 + std::pow(mre_U_2, 2)*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_3203 = -mim_U_2*mim_U_3*n_U_0*n_U_3 + mre_U_2*mre_U_3*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_3212 = -std::pow(mim_U_2, 2)*n_U_1*n_U_3 + std::pow(mre_U_2, 2)*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_3213 = -mim_U_2*mim_U_3*n_U_1*n_U_3 + mre_U_2*mre_U_3*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_3222 = -std::pow(mim_U_2, 2)*n_U_2*n_U_3 + std::pow(mre_U_2, 2)*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_3223 = -mim_U_2*mim_U_3*n_U_2*n_U_3 + mre_U_2*mre_U_3*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_3232 = -std::pow(mim_U_2, 2)*std::pow(n_U_3, 2) + std::pow(mre_U_2, 2)*std::pow(n_U_3, 2);
        amrex::Real nmnmre_UUUU_3233 = -mim_U_2*mim_U_3*std::pow(n_U_3, 2) + mre_U_2*mre_U_3*std::pow(n_U_3, 2);
        amrex::Real nmnmre_UUUU_3303 = -std::pow(mim_U_3, 2)*n_U_0*n_U_3 + std::pow(mre_U_3, 2)*n_U_0*n_U_3;
        amrex::Real nmnmre_UUUU_3313 = -std::pow(mim_U_3, 2)*n_U_1*n_U_3 + std::pow(mre_U_3, 2)*n_U_1*n_U_3;
        amrex::Real nmnmre_UUUU_3323 = -std::pow(mim_U_3, 2)*n_U_2*n_U_3 + std::pow(mre_U_3, 2)*n_U_2*n_U_3;
        amrex::Real nmnmre_UUUU_3333 = -std::pow(mim_U_3, 2)*std::pow(n_U_3, 2) + std::pow(mre_U_3, 2)*std::pow(n_U_3, 2);

        amrex::Real nmnmim_UUUU_0000 = -2*mim_U_0*mre_U_0*std::pow(n_U_0, 2);
        amrex::Real nmnmim_UUUU_0001 = -mim_U_0*mre_U_1*std::pow(n_U_0, 2) - mim_U_1*mre_U_0*std::pow(n_U_0, 2);
        amrex::Real nmnmim_UUUU_0002 = -mim_U_0*mre_U_2*std::pow(n_U_0, 2) - mim_U_2*mre_U_0*std::pow(n_U_0, 2);
        amrex::Real nmnmim_UUUU_0003 = -mim_U_0*mre_U_3*std::pow(n_U_0, 2) - mim_U_3*mre_U_0*std::pow(n_U_0, 2);
        amrex::Real nmnmim_UUUU_0010 = -2*mim_U_0*mre_U_0*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_0011 = -mim_U_0*mre_U_1*n_U_0*n_U_1 - mim_U_1*mre_U_0*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_0012 = -mim_U_0*mre_U_2*n_U_0*n_U_1 - mim_U_2*mre_U_0*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_0013 = -mim_U_0*mre_U_3*n_U_0*n_U_1 - mim_U_3*mre_U_0*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_0020 = -2*mim_U_0*mre_U_0*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_0021 = -mim_U_0*mre_U_1*n_U_0*n_U_2 - mim_U_1*mre_U_0*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_0022 = -mim_U_0*mre_U_2*n_U_0*n_U_2 - mim_U_2*mre_U_0*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_0023 = -mim_U_0*mre_U_3*n_U_0*n_U_2 - mim_U_3*mre_U_0*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_0030 = -2*mim_U_0*mre_U_0*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_0031 = -mim_U_0*mre_U_1*n_U_0*n_U_3 - mim_U_1*mre_U_0*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_0032 = -mim_U_0*mre_U_2*n_U_0*n_U_3 - mim_U_2*mre_U_0*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_0033 = -mim_U_0*mre_U_3*n_U_0*n_U_3 - mim_U_3*mre_U_0*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_0101 = -2*mim_U_1*mre_U_1*std::pow(n_U_0, 2);
        amrex::Real nmnmim_UUUU_0102 = -mim_U_1*mre_U_2*std::pow(n_U_0, 2) - mim_U_2*mre_U_1*std::pow(n_U_0, 2);
        amrex::Real nmnmim_UUUU_0103 = -mim_U_1*mre_U_3*std::pow(n_U_0, 2) - mim_U_3*mre_U_1*std::pow(n_U_0, 2);
        amrex::Real nmnmim_UUUU_0111 = -2*mim_U_1*mre_U_1*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_0112 = -mim_U_1*mre_U_2*n_U_0*n_U_1 - mim_U_2*mre_U_1*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_0113 = -mim_U_1*mre_U_3*n_U_0*n_U_1 - mim_U_3*mre_U_1*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_0121 = -2*mim_U_1*mre_U_1*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_0122 = -mim_U_1*mre_U_2*n_U_0*n_U_2 - mim_U_2*mre_U_1*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_0123 = -mim_U_1*mre_U_3*n_U_0*n_U_2 - mim_U_3*mre_U_1*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_0131 = -2*mim_U_1*mre_U_1*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_0132 = -mim_U_1*mre_U_2*n_U_0*n_U_3 - mim_U_2*mre_U_1*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_0133 = -mim_U_1*mre_U_3*n_U_0*n_U_3 - mim_U_3*mre_U_1*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_0202 = -2*mim_U_2*mre_U_2*std::pow(n_U_0, 2);
        amrex::Real nmnmim_UUUU_0203 = -mim_U_2*mre_U_3*std::pow(n_U_0, 2) - mim_U_3*mre_U_2*std::pow(n_U_0, 2);
        amrex::Real nmnmim_UUUU_0212 = -2*mim_U_2*mre_U_2*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_0213 = -mim_U_2*mre_U_3*n_U_0*n_U_1 - mim_U_3*mre_U_2*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_0222 = -2*mim_U_2*mre_U_2*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_0223 = -mim_U_2*mre_U_3*n_U_0*n_U_2 - mim_U_3*mre_U_2*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_0232 = -2*mim_U_2*mre_U_2*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_0233 = -mim_U_2*mre_U_3*n_U_0*n_U_3 - mim_U_3*mre_U_2*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_0303 = -2*mim_U_3*mre_U_3*std::pow(n_U_0, 2);
        amrex::Real nmnmim_UUUU_0313 = -2*mim_U_3*mre_U_3*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_0323 = -2*mim_U_3*mre_U_3*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_0333 = -2*mim_U_3*mre_U_3*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_1000 = -2*mim_U_0*mre_U_0*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_1001 = -mim_U_0*mre_U_1*n_U_0*n_U_1 - mim_U_1*mre_U_0*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_1002 = -mim_U_0*mre_U_2*n_U_0*n_U_1 - mim_U_2*mre_U_0*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_1003 = -mim_U_0*mre_U_3*n_U_0*n_U_1 - mim_U_3*mre_U_0*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_1010 = -2*mim_U_0*mre_U_0*std::pow(n_U_1, 2);
        amrex::Real nmnmim_UUUU_1011 = -mim_U_0*mre_U_1*std::pow(n_U_1, 2) - mim_U_1*mre_U_0*std::pow(n_U_1, 2);
        amrex::Real nmnmim_UUUU_1012 = -mim_U_0*mre_U_2*std::pow(n_U_1, 2) - mim_U_2*mre_U_0*std::pow(n_U_1, 2);
        amrex::Real nmnmim_UUUU_1013 = -mim_U_0*mre_U_3*std::pow(n_U_1, 2) - mim_U_3*mre_U_0*std::pow(n_U_1, 2);
        amrex::Real nmnmim_UUUU_1020 = -2*mim_U_0*mre_U_0*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_1021 = -mim_U_0*mre_U_1*n_U_1*n_U_2 - mim_U_1*mre_U_0*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_1022 = -mim_U_0*mre_U_2*n_U_1*n_U_2 - mim_U_2*mre_U_0*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_1023 = -mim_U_0*mre_U_3*n_U_1*n_U_2 - mim_U_3*mre_U_0*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_1030 = -2*mim_U_0*mre_U_0*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_1031 = -mim_U_0*mre_U_1*n_U_1*n_U_3 - mim_U_1*mre_U_0*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_1032 = -mim_U_0*mre_U_2*n_U_1*n_U_3 - mim_U_2*mre_U_0*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_1033 = -mim_U_0*mre_U_3*n_U_1*n_U_3 - mim_U_3*mre_U_0*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_1101 = -2*mim_U_1*mre_U_1*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_1102 = -mim_U_1*mre_U_2*n_U_0*n_U_1 - mim_U_2*mre_U_1*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_1103 = -mim_U_1*mre_U_3*n_U_0*n_U_1 - mim_U_3*mre_U_1*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_1111 = -2*mim_U_1*mre_U_1*std::pow(n_U_1, 2);
        amrex::Real nmnmim_UUUU_1112 = -mim_U_1*mre_U_2*std::pow(n_U_1, 2) - mim_U_2*mre_U_1*std::pow(n_U_1, 2);
        amrex::Real nmnmim_UUUU_1113 = -mim_U_1*mre_U_3*std::pow(n_U_1, 2) - mim_U_3*mre_U_1*std::pow(n_U_1, 2);
        amrex::Real nmnmim_UUUU_1121 = -2*mim_U_1*mre_U_1*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_1122 = -mim_U_1*mre_U_2*n_U_1*n_U_2 - mim_U_2*mre_U_1*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_1123 = -mim_U_1*mre_U_3*n_U_1*n_U_2 - mim_U_3*mre_U_1*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_1131 = -2*mim_U_1*mre_U_1*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_1132 = -mim_U_1*mre_U_2*n_U_1*n_U_3 - mim_U_2*mre_U_1*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_1133 = -mim_U_1*mre_U_3*n_U_1*n_U_3 - mim_U_3*mre_U_1*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_1202 = -2*mim_U_2*mre_U_2*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_1203 = -mim_U_2*mre_U_3*n_U_0*n_U_1 - mim_U_3*mre_U_2*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_1212 = -2*mim_U_2*mre_U_2*std::pow(n_U_1, 2);
        amrex::Real nmnmim_UUUU_1213 = -mim_U_2*mre_U_3*std::pow(n_U_1, 2) - mim_U_3*mre_U_2*std::pow(n_U_1, 2);
        amrex::Real nmnmim_UUUU_1222 = -2*mim_U_2*mre_U_2*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_1223 = -mim_U_2*mre_U_3*n_U_1*n_U_2 - mim_U_3*mre_U_2*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_1232 = -2*mim_U_2*mre_U_2*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_1233 = -mim_U_2*mre_U_3*n_U_1*n_U_3 - mim_U_3*mre_U_2*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_1303 = -2*mim_U_3*mre_U_3*n_U_0*n_U_1;
        amrex::Real nmnmim_UUUU_1313 = -2*mim_U_3*mre_U_3*std::pow(n_U_1, 2);
        amrex::Real nmnmim_UUUU_1323 = -2*mim_U_3*mre_U_3*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_1333 = -2*mim_U_3*mre_U_3*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_2000 = -2*mim_U_0*mre_U_0*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_2001 = -mim_U_0*mre_U_1*n_U_0*n_U_2 - mim_U_1*mre_U_0*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_2002 = -mim_U_0*mre_U_2*n_U_0*n_U_2 - mim_U_2*mre_U_0*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_2003 = -mim_U_0*mre_U_3*n_U_0*n_U_2 - mim_U_3*mre_U_0*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_2010 = -2*mim_U_0*mre_U_0*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_2011 = -mim_U_0*mre_U_1*n_U_1*n_U_2 - mim_U_1*mre_U_0*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_2012 = -mim_U_0*mre_U_2*n_U_1*n_U_2 - mim_U_2*mre_U_0*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_2013 = -mim_U_0*mre_U_3*n_U_1*n_U_2 - mim_U_3*mre_U_0*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_2020 = -2*mim_U_0*mre_U_0*std::pow(n_U_2, 2);
        amrex::Real nmnmim_UUUU_2021 = -mim_U_0*mre_U_1*std::pow(n_U_2, 2) - mim_U_1*mre_U_0*std::pow(n_U_2, 2);
        amrex::Real nmnmim_UUUU_2022 = -mim_U_0*mre_U_2*std::pow(n_U_2, 2) - mim_U_2*mre_U_0*std::pow(n_U_2, 2);
        amrex::Real nmnmim_UUUU_2023 = -mim_U_0*mre_U_3*std::pow(n_U_2, 2) - mim_U_3*mre_U_0*std::pow(n_U_2, 2);
        amrex::Real nmnmim_UUUU_2030 = -2*mim_U_0*mre_U_0*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_2031 = -mim_U_0*mre_U_1*n_U_2*n_U_3 - mim_U_1*mre_U_0*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_2032 = -mim_U_0*mre_U_2*n_U_2*n_U_3 - mim_U_2*mre_U_0*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_2033 = -mim_U_0*mre_U_3*n_U_2*n_U_3 - mim_U_3*mre_U_0*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_2101 = -2*mim_U_1*mre_U_1*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_2102 = -mim_U_1*mre_U_2*n_U_0*n_U_2 - mim_U_2*mre_U_1*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_2103 = -mim_U_1*mre_U_3*n_U_0*n_U_2 - mim_U_3*mre_U_1*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_2111 = -2*mim_U_1*mre_U_1*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_2112 = -mim_U_1*mre_U_2*n_U_1*n_U_2 - mim_U_2*mre_U_1*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_2113 = -mim_U_1*mre_U_3*n_U_1*n_U_2 - mim_U_3*mre_U_1*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_2121 = -2*mim_U_1*mre_U_1*std::pow(n_U_2, 2);
        amrex::Real nmnmim_UUUU_2122 = -mim_U_1*mre_U_2*std::pow(n_U_2, 2) - mim_U_2*mre_U_1*std::pow(n_U_2, 2);
        amrex::Real nmnmim_UUUU_2123 = -mim_U_1*mre_U_3*std::pow(n_U_2, 2) - mim_U_3*mre_U_1*std::pow(n_U_2, 2);
        amrex::Real nmnmim_UUUU_2131 = -2*mim_U_1*mre_U_1*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_2132 = -mim_U_1*mre_U_2*n_U_2*n_U_3 - mim_U_2*mre_U_1*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_2133 = -mim_U_1*mre_U_3*n_U_2*n_U_3 - mim_U_3*mre_U_1*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_2202 = -2*mim_U_2*mre_U_2*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_2203 = -mim_U_2*mre_U_3*n_U_0*n_U_2 - mim_U_3*mre_U_2*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_2212 = -2*mim_U_2*mre_U_2*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_2213 = -mim_U_2*mre_U_3*n_U_1*n_U_2 - mim_U_3*mre_U_2*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_2222 = -2*mim_U_2*mre_U_2*std::pow(n_U_2, 2);
        amrex::Real nmnmim_UUUU_2223 = -mim_U_2*mre_U_3*std::pow(n_U_2, 2) - mim_U_3*mre_U_2*std::pow(n_U_2, 2);
        amrex::Real nmnmim_UUUU_2232 = -2*mim_U_2*mre_U_2*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_2233 = -mim_U_2*mre_U_3*n_U_2*n_U_3 - mim_U_3*mre_U_2*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_2303 = -2*mim_U_3*mre_U_3*n_U_0*n_U_2;
        amrex::Real nmnmim_UUUU_2313 = -2*mim_U_3*mre_U_3*n_U_1*n_U_2;
        amrex::Real nmnmim_UUUU_2323 = -2*mim_U_3*mre_U_3*std::pow(n_U_2, 2);
        amrex::Real nmnmim_UUUU_2333 = -2*mim_U_3*mre_U_3*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_3000 = -2*mim_U_0*mre_U_0*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_3001 = -mim_U_0*mre_U_1*n_U_0*n_U_3 - mim_U_1*mre_U_0*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_3002 = -mim_U_0*mre_U_2*n_U_0*n_U_3 - mim_U_2*mre_U_0*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_3003 = -mim_U_0*mre_U_3*n_U_0*n_U_3 - mim_U_3*mre_U_0*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_3010 = -2*mim_U_0*mre_U_0*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_3011 = -mim_U_0*mre_U_1*n_U_1*n_U_3 - mim_U_1*mre_U_0*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_3012 = -mim_U_0*mre_U_2*n_U_1*n_U_3 - mim_U_2*mre_U_0*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_3013 = -mim_U_0*mre_U_3*n_U_1*n_U_3 - mim_U_3*mre_U_0*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_3020 = -2*mim_U_0*mre_U_0*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_3021 = -mim_U_0*mre_U_1*n_U_2*n_U_3 - mim_U_1*mre_U_0*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_3022 = -mim_U_0*mre_U_2*n_U_2*n_U_3 - mim_U_2*mre_U_0*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_3023 = -mim_U_0*mre_U_3*n_U_2*n_U_3 - mim_U_3*mre_U_0*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_3030 = -2*mim_U_0*mre_U_0*std::pow(n_U_3, 2);
        amrex::Real nmnmim_UUUU_3031 = -mim_U_0*mre_U_1*std::pow(n_U_3, 2) - mim_U_1*mre_U_0*std::pow(n_U_3, 2);
        amrex::Real nmnmim_UUUU_3032 = -mim_U_0*mre_U_2*std::pow(n_U_3, 2) - mim_U_2*mre_U_0*std::pow(n_U_3, 2);
        amrex::Real nmnmim_UUUU_3033 = -mim_U_0*mre_U_3*std::pow(n_U_3, 2) - mim_U_3*mre_U_0*std::pow(n_U_3, 2);
        amrex::Real nmnmim_UUUU_3101 = -2*mim_U_1*mre_U_1*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_3102 = -mim_U_1*mre_U_2*n_U_0*n_U_3 - mim_U_2*mre_U_1*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_3103 = -mim_U_1*mre_U_3*n_U_0*n_U_3 - mim_U_3*mre_U_1*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_3111 = -2*mim_U_1*mre_U_1*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_3112 = -mim_U_1*mre_U_2*n_U_1*n_U_3 - mim_U_2*mre_U_1*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_3113 = -mim_U_1*mre_U_3*n_U_1*n_U_3 - mim_U_3*mre_U_1*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_3121 = -2*mim_U_1*mre_U_1*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_3122 = -mim_U_1*mre_U_2*n_U_2*n_U_3 - mim_U_2*mre_U_1*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_3123 = -mim_U_1*mre_U_3*n_U_2*n_U_3 - mim_U_3*mre_U_1*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_3131 = -2*mim_U_1*mre_U_1*std::pow(n_U_3, 2);
        amrex::Real nmnmim_UUUU_3132 = -mim_U_1*mre_U_2*std::pow(n_U_3, 2) - mim_U_2*mre_U_1*std::pow(n_U_3, 2);
        amrex::Real nmnmim_UUUU_3133 = -mim_U_1*mre_U_3*std::pow(n_U_3, 2) - mim_U_3*mre_U_1*std::pow(n_U_3, 2);
        amrex::Real nmnmim_UUUU_3202 = -2*mim_U_2*mre_U_2*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_3203 = -mim_U_2*mre_U_3*n_U_0*n_U_3 - mim_U_3*mre_U_2*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_3212 = -2*mim_U_2*mre_U_2*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_3213 = -mim_U_2*mre_U_3*n_U_1*n_U_3 - mim_U_3*mre_U_2*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_3222 = -2*mim_U_2*mre_U_2*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_3223 = -mim_U_2*mre_U_3*n_U_2*n_U_3 - mim_U_3*mre_U_2*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_3232 = -2*mim_U_2*mre_U_2*std::pow(n_U_3, 2);
        amrex::Real nmnmim_UUUU_3233 = -mim_U_2*mre_U_3*std::pow(n_U_3, 2) - mim_U_3*mre_U_2*std::pow(n_U_3, 2);
        amrex::Real nmnmim_UUUU_3303 = -2*mim_U_3*mre_U_3*n_U_0*n_U_3;
        amrex::Real nmnmim_UUUU_3313 = -2*mim_U_3*mre_U_3*n_U_1*n_U_3;
        amrex::Real nmnmim_UUUU_3323 = -2*mim_U_3*mre_U_3*n_U_2*n_U_3;
        amrex::Real nmnmim_UUUU_3333 = -2*mim_U_3*mre_U_3*std::pow(n_U_3, 2);

        amrex::Real rank4Part_LLLL_0101 = K_LL_00*K_LL_11 - std::pow(K_LL_01, 2) + Riem_LLLL_0101;
        amrex::Real rank4Part_LLLL_0102 = K_LL_00*K_LL_12 - K_LL_01*K_LL_02 + Riem_LLLL_0102;
        amrex::Real rank4Part_LLLL_0112 = K_LL_01*K_LL_12 - K_LL_02*K_LL_11 + Riem_LLLL_0112;
        amrex::Real rank4Part_LLLL_0201 = K_LL_00*K_LL_12 - K_LL_01*K_LL_02 + Riem_LLLL_0201;
        amrex::Real rank4Part_LLLL_0202 = K_LL_00*K_LL_22 - std::pow(K_LL_02, 2) + Riem_LLLL_0202;
        amrex::Real rank4Part_LLLL_0212 = K_LL_01*K_LL_22 - K_LL_02*K_LL_12 + Riem_LLLL_0212;
        amrex::Real rank4Part_LLLL_1201 = K_LL_01*K_LL_12 - K_LL_02*K_LL_11 + Riem_LLLL_1201;
        amrex::Real rank4Part_LLLL_1202 = K_LL_01*K_LL_22 - K_LL_02*K_LL_12 + Riem_LLLL_1202;
        amrex::Real rank4Part_LLLL_1212 = K_LL_11*K_LL_22 - std::pow(K_LL_12, 2) + Riem_LLLL_1212;

        amrex::Real dDK_LLL_000 = (-K_LL_00*dDchi_L_0 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_000 + dDAtilde_LLL_000 + 0.33333333333333331*gamtilde_LL_00*(dDKhat_L_0 + 2*dDtheta_L_0))/chi;
        amrex::Real dDK_LLL_001 = (-K_LL_00*dDchi_L_1 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_001 + dDAtilde_LLL_001 + 0.33333333333333331*gamtilde_LL_00*(dDKhat_L_1 + 2*dDtheta_L_1))/chi;
        amrex::Real dDK_LLL_002 = (-K_LL_00*dDchi_L_2 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_002 + dDAtilde_LLL_002 + 0.33333333333333331*gamtilde_LL_00*(dDKhat_L_2 + 2*dDtheta_L_2))/chi;
        amrex::Real dDK_LLL_010 = (-K_LL_01*dDchi_L_0 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_010 + dDAtilde_LLL_010 + 0.33333333333333331*gamtilde_LL_01*(dDKhat_L_0 + 2*dDtheta_L_0))/chi;
        amrex::Real dDK_LLL_011 = (-K_LL_01*dDchi_L_1 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_011 + dDAtilde_LLL_011 + 0.33333333333333331*gamtilde_LL_01*(dDKhat_L_1 + 2*dDtheta_L_1))/chi;
        amrex::Real dDK_LLL_012 = (-K_LL_01*dDchi_L_2 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_012 + dDAtilde_LLL_012 + 0.33333333333333331*gamtilde_LL_01*(dDKhat_L_2 + 2*dDtheta_L_2))/chi;
        amrex::Real dDK_LLL_020 = (-K_LL_02*dDchi_L_0 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_020 + dDAtilde_LLL_020 + 0.33333333333333331*gamtilde_LL_02*(dDKhat_L_0 + 2*dDtheta_L_0))/chi;
        amrex::Real dDK_LLL_021 = (-K_LL_02*dDchi_L_1 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_021 + dDAtilde_LLL_021 + 0.33333333333333331*gamtilde_LL_02*(dDKhat_L_1 + 2*dDtheta_L_1))/chi;
        amrex::Real dDK_LLL_022 = (-K_LL_02*dDchi_L_2 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_022 + dDAtilde_LLL_022 + 0.33333333333333331*gamtilde_LL_02*(dDKhat_L_2 + 2*dDtheta_L_2))/chi;
        amrex::Real dDK_LLL_110 = (-K_LL_11*dDchi_L_0 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_110 + dDAtilde_LLL_110 + 0.33333333333333331*gamtilde_LL_11*(dDKhat_L_0 + 2*dDtheta_L_0))/chi;
        amrex::Real dDK_LLL_111 = (-K_LL_11*dDchi_L_1 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_111 + dDAtilde_LLL_111 + 0.33333333333333331*gamtilde_LL_11*(dDKhat_L_1 + 2*dDtheta_L_1))/chi;
        amrex::Real dDK_LLL_112 = (-K_LL_11*dDchi_L_2 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_112 + dDAtilde_LLL_112 + 0.33333333333333331*gamtilde_LL_11*(dDKhat_L_2 + 2*dDtheta_L_2))/chi;
        amrex::Real dDK_LLL_120 = (-K_LL_12*dDchi_L_0 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_120 + dDAtilde_LLL_120 + 0.33333333333333331*gamtilde_LL_12*(dDKhat_L_0 + 2*dDtheta_L_0))/chi;
        amrex::Real dDK_LLL_121 = (-K_LL_12*dDchi_L_1 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_121 + dDAtilde_LLL_121 + 0.33333333333333331*gamtilde_LL_12*(dDKhat_L_1 + 2*dDtheta_L_1))/chi;
        amrex::Real dDK_LLL_122 = (-K_LL_12*dDchi_L_2 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_122 + dDAtilde_LLL_122 + 0.33333333333333331*gamtilde_LL_12*(dDKhat_L_2 + 2*dDtheta_L_2))/chi;
        amrex::Real dDK_LLL_220 = (-K_LL_22*dDchi_L_0 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_220 + dDAtilde_LLL_220 + 0.33333333333333331*gamtilde_LL_22*(dDKhat_L_0 + 2*dDtheta_L_0))/chi;
        amrex::Real dDK_LLL_221 = (-K_LL_22*dDchi_L_1 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_221 + dDAtilde_LLL_221 + 0.33333333333333331*gamtilde_LL_22*(dDKhat_L_1 + 2*dDtheta_L_1))/chi;
        amrex::Real dDK_LLL_222 = (-K_LL_22*dDchi_L_2 + 0.33333333333333331*Ksclr*dDgamtilde_LLL_222 + dDAtilde_LLL_222 + 0.33333333333333331*gamtilde_LL_22*(dDKhat_L_2 + 2*dDtheta_L_2))/chi;

        amrex::Real rank3Part_LLL_001 = -4*ConformalChristoffelchigamtilde_ULL_000*K_LL_01 + 4*ConformalChristoffelchigamtilde_ULL_001*K_LL_00 - 4*ConformalChristoffelchigamtilde_ULL_100*K_LL_11 + 4*ConformalChristoffelchigamtilde_ULL_101*K_LL_01 - 4*ConformalChristoffelchigamtilde_ULL_200*K_LL_12 + 4*ConformalChristoffelchigamtilde_ULL_201*K_LL_02 - 4*dDK_LLL_001 + 4*dDK_LLL_010;
        amrex::Real rank3Part_LLL_002 = -4*ConformalChristoffelchigamtilde_ULL_000*K_LL_02 + 4*ConformalChristoffelchigamtilde_ULL_002*K_LL_00 - 4*ConformalChristoffelchigamtilde_ULL_100*K_LL_12 + 4*ConformalChristoffelchigamtilde_ULL_102*K_LL_01 - 4*ConformalChristoffelchigamtilde_ULL_200*K_LL_22 + 4*ConformalChristoffelchigamtilde_ULL_202*K_LL_02 - 4*dDK_LLL_002 + 4*dDK_LLL_020;
        amrex::Real rank3Part_LLL_012 = -4*ConformalChristoffelchigamtilde_ULL_001*K_LL_02 + 4*ConformalChristoffelchigamtilde_ULL_002*K_LL_01 - 4*ConformalChristoffelchigamtilde_ULL_101*K_LL_12 + 4*ConformalChristoffelchigamtilde_ULL_102*K_LL_11 - 4*ConformalChristoffelchigamtilde_ULL_201*K_LL_22 + 4*ConformalChristoffelchigamtilde_ULL_202*K_LL_12 - 4*dDK_LLL_012 + 4*dDK_LLL_021;
        amrex::Real rank3Part_LLL_101 = -4*ConformalChristoffelchigamtilde_ULL_001*K_LL_01 + 4*ConformalChristoffelchigamtilde_ULL_011*K_LL_00 - 4*ConformalChristoffelchigamtilde_ULL_101*K_LL_11 + 4*ConformalChristoffelchigamtilde_ULL_111*K_LL_01 - 4*ConformalChristoffelchigamtilde_ULL_201*K_LL_12 + 4*ConformalChristoffelchigamtilde_ULL_211*K_LL_02 - 4*dDK_LLL_011 + 4*dDK_LLL_110;
        amrex::Real rank3Part_LLL_102 = -4*ConformalChristoffelchigamtilde_ULL_001*K_LL_02 + 4*ConformalChristoffelchigamtilde_ULL_012*K_LL_00 - 4*ConformalChristoffelchigamtilde_ULL_101*K_LL_12 + 4*ConformalChristoffelchigamtilde_ULL_112*K_LL_01 - 4*ConformalChristoffelchigamtilde_ULL_201*K_LL_22 + 4*ConformalChristoffelchigamtilde_ULL_212*K_LL_02 - 4*dDK_LLL_012 + 4*dDK_LLL_120;
        amrex::Real rank3Part_LLL_112 = -4*ConformalChristoffelchigamtilde_ULL_011*K_LL_02 + 4*ConformalChristoffelchigamtilde_ULL_012*K_LL_01 - 4*ConformalChristoffelchigamtilde_ULL_111*K_LL_12 + 4*ConformalChristoffelchigamtilde_ULL_112*K_LL_11 - 4*ConformalChristoffelchigamtilde_ULL_211*K_LL_22 + 4*ConformalChristoffelchigamtilde_ULL_212*K_LL_12 - 4*dDK_LLL_112 + 4*dDK_LLL_121;
        amrex::Real rank3Part_LLL_201 = -4*ConformalChristoffelchigamtilde_ULL_002*K_LL_01 + 4*ConformalChristoffelchigamtilde_ULL_012*K_LL_00 - 4*ConformalChristoffelchigamtilde_ULL_102*K_LL_11 + 4*ConformalChristoffelchigamtilde_ULL_112*K_LL_01 - 4*ConformalChristoffelchigamtilde_ULL_202*K_LL_12 + 4*ConformalChristoffelchigamtilde_ULL_212*K_LL_02 - 4*dDK_LLL_021 + 4*dDK_LLL_120;
        amrex::Real rank3Part_LLL_202 = -4*ConformalChristoffelchigamtilde_ULL_002*K_LL_02 + 4*ConformalChristoffelchigamtilde_ULL_022*K_LL_00 - 4*ConformalChristoffelchigamtilde_ULL_102*K_LL_12 + 4*ConformalChristoffelchigamtilde_ULL_122*K_LL_01 - 4*ConformalChristoffelchigamtilde_ULL_202*K_LL_22 + 4*ConformalChristoffelchigamtilde_ULL_222*K_LL_02 - 4*dDK_LLL_022 + 4*dDK_LLL_220;
        amrex::Real rank3Part_LLL_212 = -4*ConformalChristoffelchigamtilde_ULL_012*K_LL_02 + 4*ConformalChristoffelchigamtilde_ULL_022*K_LL_01 - 4*ConformalChristoffelchigamtilde_ULL_112*K_LL_12 + 4*ConformalChristoffelchigamtilde_ULL_122*K_LL_11 - 4*ConformalChristoffelchigamtilde_ULL_212*K_LL_22 + 4*ConformalChristoffelchigamtilde_ULL_222*K_LL_12 - 4*dDK_LLL_122 + 4*dDK_LLL_221;

        amrex::Real rank2Part_LL_00 = -12*K_LL_00*K_LL_02*invgam_UU_20 + 4*K_LL_00*Ksclr - 12*K_LL_01*K_LL_02*invgam_UU_21 - 12*std::pow(K_LL_02, 2)*invgam_UU_22 + 4*Ricci_LL_00;
        amrex::Real rank2Part_LL_01 = -12*K_LL_01*K_LL_02*invgam_UU_20 + 4*K_LL_01*Ksclr - 12*K_LL_02*K_LL_11*invgam_UU_21 - 12*K_LL_02*K_LL_12*invgam_UU_22 + 4*Ricci_LL_01;
        amrex::Real rank2Part_LL_02 = -12*std::pow(K_LL_02, 2)*invgam_UU_20 - 12*K_LL_02*K_LL_12*invgam_UU_21 - 12*K_LL_02*K_LL_22*invgam_UU_22 + 4*K_LL_02*Ksclr + 4*Ricci_LL_02;
        amrex::Real rank2Part_LL_11 = -12*K_LL_01*K_LL_12*invgam_UU_20 - 12*K_LL_11*K_LL_12*invgam_UU_21 + 4*K_LL_11*Ksclr - 12*std::pow(K_LL_12, 2)*invgam_UU_22 + 4*Ricci_LL_11;
        amrex::Real rank2Part_LL_12 = -12*K_LL_02*K_LL_12*invgam_UU_20 - 12*std::pow(K_LL_12, 2)*invgam_UU_21 - 12*K_LL_12*K_LL_22*invgam_UU_22 + 4*K_LL_12*Ksclr + 4*Ricci_LL_12;
        amrex::Real rank2Part_LL_22 = -12*K_LL_02*K_LL_22*invgam_UU_20 - 12*K_LL_12*K_LL_22*invgam_UU_21 - 12*std::pow(K_LL_22, 2)*invgam_UU_22 + 4*K_LL_22*Ksclr + 4*Ricci_LL_22;

        amrex::Real Psi4re = (-nmnmre_UUUU_0011*rank4Part_LLLL_0101 - nmnmre_UUUU_0012*rank4Part_LLLL_0201 - nmnmre_UUUU_0021*rank4Part_LLLL_0102 - nmnmre_UUUU_0022*rank4Part_LLLL_0202 + nmnmre_UUUU_0101*rank4Part_LLLL_0101 + nmnmre_UUUU_0102*rank4Part_LLLL_0102 + nmnmre_UUUU_0102*rank4Part_LLLL_0201 + nmnmre_UUUU_0112*rank4Part_LLLL_0112 - nmnmre_UUUU_0121*rank4Part_LLLL_0112 - nmnmre_UUUU_0122*rank4Part_LLLL_0212 + nmnmre_UUUU_0202*rank4Part_LLLL_0202 + nmnmre_UUUU_0212*rank4Part_LLLL_0212 - nmnmre_UUUU_1001*rank4Part_LLLL_0101 - nmnmre_UUUU_1002*rank4Part_LLLL_0102 + nmnmre_UUUU_1010*rank4Part_LLLL_0101 - nmnmre_UUUU_1012*rank4Part_LLLL_0112 - nmnmre_UUUU_1012*rank4Part_LLLL_1201 + nmnmre_UUUU_1020*rank4Part_LLLL_0102 + nmnmre_UUUU_1021*rank4Part_LLLL_0112 - nmnmre_UUUU_1022*rank4Part_LLLL_1202 + nmnmre_UUUU_1102*rank4Part_LLLL_1201 - nmnmre_UUUU_1122*rank4Part_LLLL_1212 + nmnmre_UUUU_1202*rank4Part_LLLL_1202 + nmnmre_UUUU_1212*rank4Part_LLLL_1212 - nmnmre_UUUU_2001*rank4Part_LLLL_0201 - nmnmre_UUUU_2002*rank4Part_LLLL_0202 + nmnmre_UUUU_2010*rank4Part_LLLL_0201 + nmnmre_UUUU_2011*rank4Part_LLLL_1201 - nmnmre_UUUU_2012*rank4Part_LLLL_0212 + nmnmre_UUUU_2020*rank4Part_LLLL_0202 + nmnmre_UUUU_2021*rank4Part_LLLL_0212 + nmnmre_UUUU_2021*rank4Part_LLLL_1202 - nmnmre_UUUU_2101*rank4Part_LLLL_1201 - nmnmre_UUUU_2102*rank4Part_LLLL_1202 - nmnmre_UUUU_2112*rank4Part_LLLL_1212 + nmnmre_UUUU_2121*rank4Part_LLLL_1212 + 0.25*rank2Part_LL_00*(-nmnmre_UUUU_0033 + nmnmre_UUUU_0303 - nmnmre_UUUU_3003 + nmnmre_UUUU_3030) + 0.25*rank2Part_LL_01*(-nmnmre_UUUU_0133 + nmnmre_UUUU_1303 + nmnmre_UUUU_3031 - nmnmre_UUUU_3103) + 0.25*rank2Part_LL_01*(nmnmre_UUUU_0313 - nmnmre_UUUU_1033 - nmnmre_UUUU_3013 + nmnmre_UUUU_3031) + 0.25*rank2Part_LL_02*(-nmnmre_UUUU_0233 + nmnmre_UUUU_2303 + nmnmre_UUUU_3032 - nmnmre_UUUU_3203) + 0.25*rank2Part_LL_02*(nmnmre_UUUU_0323 - nmnmre_UUUU_2033 - nmnmre_UUUU_3023 + nmnmre_UUUU_3032) + 0.25*rank2Part_LL_11*(-nmnmre_UUUU_1133 + nmnmre_UUUU_1313 - nmnmre_UUUU_3113 + nmnmre_UUUU_3131) + 0.25*rank2Part_LL_12*(-nmnmre_UUUU_1233 + nmnmre_UUUU_2313 + nmnmre_UUUU_3132 - nmnmre_UUUU_3213) + 0.25*rank2Part_LL_12*(nmnmre_UUUU_1323 - nmnmre_UUUU_2133 - nmnmre_UUUU_3123 + nmnmre_UUUU_3132) + 0.25*rank2Part_LL_22*(-nmnmre_UUUU_2233 + nmnmre_UUUU_2323 - nmnmre_UUUU_3223 + nmnmre_UUUU_3232) - 0.5*rank3Part_LLL_001*(-nmnmre_UUUU_0013 + nmnmre_UUUU_3010) + 0.5*rank3Part_LLL_001*(-nmnmre_UUUU_0103 + nmnmre_UUUU_3001) - 0.5*rank3Part_LLL_002*(-nmnmre_UUUU_0023 + nmnmre_UUUU_3020) + 0.5*rank3Part_LLL_002*(-nmnmre_UUUU_0203 + nmnmre_UUUU_3002) - 0.5*rank3Part_LLL_012*(-nmnmre_UUUU_0123 + nmnmre_UUUU_3021) + 0.5*rank3Part_LLL_012*(-nmnmre_UUUU_0213 + nmnmre_UUUU_3012) - 0.5*rank3Part_LLL_101*(-nmnmre_UUUU_1013 + nmnmre_UUUU_3011) + 0.5*rank3Part_LLL_101*(-nmnmre_UUUU_1103 + nmnmre_UUUU_3101) - 0.5*rank3Part_LLL_102*(-nmnmre_UUUU_1023 + nmnmre_UUUU_3021) + 0.5*rank3Part_LLL_102*(-nmnmre_UUUU_1203 + nmnmre_UUUU_3102) - 0.5*rank3Part_LLL_112*(-nmnmre_UUUU_1123 + nmnmre_UUUU_3121) + 0.5*rank3Part_LLL_112*(-nmnmre_UUUU_1213 + nmnmre_UUUU_3112) - 0.5*rank3Part_LLL_201*(-nmnmre_UUUU_2013 + nmnmre_UUUU_3012) + 0.5*rank3Part_LLL_201*(-nmnmre_UUUU_2103 + nmnmre_UUUU_3102) - 0.5*rank3Part_LLL_202*(-nmnmre_UUUU_2023 + nmnmre_UUUU_3022) + 0.5*rank3Part_LLL_202*(-nmnmre_UUUU_2203 + nmnmre_UUUU_3202) - 0.5*rank3Part_LLL_212*(-nmnmre_UUUU_2123 + nmnmre_UUUU_3122) + 0.5*rank3Part_LLL_212*(-nmnmre_UUUU_2213 + nmnmre_UUUU_3212))*(rcoord > 0);
        
        diag(i, j, k, Diag::Psi4re) = Psi4re;

        amrex::Real Psi4im = (-nmnmim_UUUU_0011*rank4Part_LLLL_0101 - nmnmim_UUUU_0012*rank4Part_LLLL_0201 - nmnmim_UUUU_0021*rank4Part_LLLL_0102 - nmnmim_UUUU_0022*rank4Part_LLLL_0202 + nmnmim_UUUU_0101*rank4Part_LLLL_0101 + nmnmim_UUUU_0102*rank4Part_LLLL_0102 + nmnmim_UUUU_0102*rank4Part_LLLL_0201 + nmnmim_UUUU_0112*rank4Part_LLLL_0112 - nmnmim_UUUU_0121*rank4Part_LLLL_0112 - nmnmim_UUUU_0122*rank4Part_LLLL_0212 + nmnmim_UUUU_0202*rank4Part_LLLL_0202 + nmnmim_UUUU_0212*rank4Part_LLLL_0212 - nmnmim_UUUU_1001*rank4Part_LLLL_0101 - nmnmim_UUUU_1002*rank4Part_LLLL_0102 + nmnmim_UUUU_1010*rank4Part_LLLL_0101 - nmnmim_UUUU_1012*rank4Part_LLLL_0112 - nmnmim_UUUU_1012*rank4Part_LLLL_1201 + nmnmim_UUUU_1020*rank4Part_LLLL_0102 + nmnmim_UUUU_1021*rank4Part_LLLL_0112 - nmnmim_UUUU_1022*rank4Part_LLLL_1202 + nmnmim_UUUU_1102*rank4Part_LLLL_1201 - nmnmim_UUUU_1122*rank4Part_LLLL_1212 + nmnmim_UUUU_1202*rank4Part_LLLL_1202 + nmnmim_UUUU_1212*rank4Part_LLLL_1212 - nmnmim_UUUU_2001*rank4Part_LLLL_0201 - nmnmim_UUUU_2002*rank4Part_LLLL_0202 + nmnmim_UUUU_2010*rank4Part_LLLL_0201 + nmnmim_UUUU_2011*rank4Part_LLLL_1201 - nmnmim_UUUU_2012*rank4Part_LLLL_0212 + nmnmim_UUUU_2020*rank4Part_LLLL_0202 + nmnmim_UUUU_2021*rank4Part_LLLL_0212 + nmnmim_UUUU_2021*rank4Part_LLLL_1202 - nmnmim_UUUU_2101*rank4Part_LLLL_1201 - nmnmim_UUUU_2102*rank4Part_LLLL_1202 - nmnmim_UUUU_2112*rank4Part_LLLL_1212 + nmnmim_UUUU_2121*rank4Part_LLLL_1212 + 0.25*rank2Part_LL_00*(-nmnmim_UUUU_0033 + nmnmim_UUUU_0303 - nmnmim_UUUU_3003 + nmnmim_UUUU_3030) + 0.25*rank2Part_LL_01*(-nmnmim_UUUU_0133 + nmnmim_UUUU_1303 + nmnmim_UUUU_3031 - nmnmim_UUUU_3103) + 0.25*rank2Part_LL_01*(nmnmim_UUUU_0313 - nmnmim_UUUU_1033 - nmnmim_UUUU_3013 + nmnmim_UUUU_3031) + 0.25*rank2Part_LL_02*(-nmnmim_UUUU_0233 + nmnmim_UUUU_2303 + nmnmim_UUUU_3032 - nmnmim_UUUU_3203) + 0.25*rank2Part_LL_02*(nmnmim_UUUU_0323 - nmnmim_UUUU_2033 - nmnmim_UUUU_3023 + nmnmim_UUUU_3032) + 0.25*rank2Part_LL_11*(-nmnmim_UUUU_1133 + nmnmim_UUUU_1313 - nmnmim_UUUU_3113 + nmnmim_UUUU_3131) + 0.25*rank2Part_LL_12*(-nmnmim_UUUU_1233 + nmnmim_UUUU_2313 + nmnmim_UUUU_3132 - nmnmim_UUUU_3213) + 0.25*rank2Part_LL_12*(nmnmim_UUUU_1323 - nmnmim_UUUU_2133 - nmnmim_UUUU_3123 + nmnmim_UUUU_3132) + 0.25*rank2Part_LL_22*(-nmnmim_UUUU_2233 + nmnmim_UUUU_2323 - nmnmim_UUUU_3223 + nmnmim_UUUU_3232) - 0.5*rank3Part_LLL_001*(-nmnmim_UUUU_0013 + nmnmim_UUUU_3010) + 0.5*rank3Part_LLL_001*(-nmnmim_UUUU_0103 + nmnmim_UUUU_3001) - 0.5*rank3Part_LLL_002*(-nmnmim_UUUU_0023 + nmnmim_UUUU_3020) + 0.5*rank3Part_LLL_002*(-nmnmim_UUUU_0203 + nmnmim_UUUU_3002) - 0.5*rank3Part_LLL_012*(-nmnmim_UUUU_0123 + nmnmim_UUUU_3021) + 0.5*rank3Part_LLL_012*(-nmnmim_UUUU_0213 + nmnmim_UUUU_3012) - 0.5*rank3Part_LLL_101*(-nmnmim_UUUU_1013 + nmnmim_UUUU_3011) + 0.5*rank3Part_LLL_101*(-nmnmim_UUUU_1103 + nmnmim_UUUU_3101) - 0.5*rank3Part_LLL_102*(-nmnmim_UUUU_1023 + nmnmim_UUUU_3021) + 0.5*rank3Part_LLL_102*(-nmnmim_UUUU_1203 + nmnmim_UUUU_3102) - 0.5*rank3Part_LLL_112*(-nmnmim_UUUU_1123 + nmnmim_UUUU_3121) + 0.5*rank3Part_LLL_112*(-nmnmim_UUUU_1213 + nmnmim_UUUU_3112) - 0.5*rank3Part_LLL_201*(-nmnmim_UUUU_2013 + nmnmim_UUUU_3012) + 0.5*rank3Part_LLL_201*(-nmnmim_UUUU_2103 + nmnmim_UUUU_3102) - 0.5*rank3Part_LLL_202*(-nmnmim_UUUU_2023 + nmnmim_UUUU_3022) + 0.5*rank3Part_LLL_202*(-nmnmim_UUUU_2203 + nmnmim_UUUU_3202) - 0.5*rank3Part_LLL_212*(-nmnmim_UUUU_2123 + nmnmim_UUUU_3122) + 0.5*rank3Part_LLL_212*(-nmnmim_UUUU_2213 + nmnmim_UUUU_3212))*(rcoord > 0);
        
        diag(i, j, k, Diag::Psi4im) = Psi4im;*/
        
        diag(i, j, k, Diag::AbsKhat) = Khat;
        diag(i, j, k, Diag::AbsKsclr) = Khat + 2*theta;
        diag(i, j, k, Diag::AbsAlpha) = std::fabs(1-alpha);
        
        
        //diag(i, j, k, Diag::rPsi4Real) = rcoord*std::fabs(Psi4re);
        //diag(i, j, k, Diag::rPsi4Imaginary) = rcoord*std::fabs(Psi4im);
        //diag(i, j, k, Diag::rModPsi4) = rcoord*std::sqrt(Psi4re*Psi4re + Psi4im*Psi4im);
        
}
#endif
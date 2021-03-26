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
        
        amrex::Real v = 5.0;

        amrex::Real lam = 0.25;

        amrex::Real Omega = 0.25;

        amrex::Real e = 0.3/std::sqrt(2);

        amrex::Real PhiReal = state_fab(i, j, k, Idx::PhiReal);

        amrex::Real PhiImaginary = state_fab(i, j, k, Idx::PhiImaginary);

        amrex::Real Chi = state_fab(i, j, k, Idx::Chi);

        amrex::Real Chi0 = state_fab(i, j, k, Idx::Chi0);

        amrex::Real Chi1 = state_fab(i, j, k, Idx::Chi1);

        amrex::Real A_L_0 = state_fab(i, j, k, Idx::A_L_0);
        amrex::Real A_L_1 = state_fab(i, j, k, Idx::A_L_1);
        amrex::Real A_L_2 = state_fab(i, j, k, Idx::A_L_2);

        amrex::Real AT = state_fab(i, j, k, Idx::AT);

        amrex::Real PiPhiReal = state_fab(i, j, k, Idx::PiPhiReal);

        amrex::Real PiPhiImaginary = state_fab(i, j, k, Idx::PiPhiImaginary);

        amrex::Real PiChi = state_fab(i, j, k, Idx::PiChi);

        amrex::Real PiChi0 = state_fab(i, j, k, Idx::PiChi0);

        amrex::Real PiChi1 = state_fab(i, j, k, Idx::PiChi1);

        amrex::Real PiA_L_0 = state_fab(i, j, k, Idx::PiA_L_0);
        amrex::Real PiA_L_1 = state_fab(i, j, k, Idx::PiA_L_1);
        amrex::Real PiA_L_2 = state_fab(i, j, k, Idx::PiA_L_2);

        amrex::Real PiAT = state_fab(i, j, k, Idx::PiAT);

        amrex::Real dDPhiReal_L_0 = ((1.0/2.0)*state_fab(i + 1, j, k, Idx::PhiReal) - 1.0/2.0*state_fab(i - 1, j, k, Idx::PhiReal))/dx[0];
        amrex::Real dDPhiReal_L_1 = ((1.0/2.0)*state_fab(i, j + 1, k, Idx::PhiReal) - 1.0/2.0*state_fab(i, j - 1, k, Idx::PhiReal))/dx[1];
        amrex::Real dDPhiReal_L_2 = ((1.0/2.0)*state_fab(i, j, k + 1, Idx::PhiReal) - 1.0/2.0*state_fab(i, j, k - 1, Idx::PhiReal))/dx[2];

        amrex::Real dDPhiImaginary_L_0 = ((1.0/2.0)*state_fab(i + 1, j, k, Idx::PhiImaginary) - 1.0/2.0*state_fab(i - 1, j, k, Idx::PhiImaginary))/dx[0];
        amrex::Real dDPhiImaginary_L_1 = ((1.0/2.0)*state_fab(i, j + 1, k, Idx::PhiImaginary) - 1.0/2.0*state_fab(i, j - 1, k, Idx::PhiImaginary))/dx[1];
        amrex::Real dDPhiImaginary_L_2 = ((1.0/2.0)*state_fab(i, j, k + 1, Idx::PhiImaginary) - 1.0/2.0*state_fab(i, j, k - 1, Idx::PhiImaginary))/dx[2];

        amrex::Real dDDPhiReal_LL_00 = (state_fab(i + 1, j, k, Idx::PhiReal) + state_fab(i - 1, j, k, Idx::PhiReal) - 2*state_fab(i, j, k, Idx::PhiReal))/std::pow(dx[0], 2);
        amrex::Real dDDPhiReal_LL_01 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::PhiReal) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::PhiReal) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::PhiReal) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::PhiReal))/(dx[0]*dx[1]);
        amrex::Real dDDPhiReal_LL_02 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::PhiReal) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::PhiReal) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::PhiReal) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::PhiReal))/(dx[0]*dx[2]);
        amrex::Real dDDPhiReal_LL_11 = (state_fab(i, j + 1, k, Idx::PhiReal) + state_fab(i, j - 1, k, Idx::PhiReal) - 2*state_fab(i, j, k, Idx::PhiReal))/std::pow(dx[1], 2);
        amrex::Real dDDPhiReal_LL_12 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::PhiReal) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::PhiReal) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::PhiReal) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::PhiReal))/(dx[1]*dx[2]);
        amrex::Real dDDPhiReal_LL_22 = (state_fab(i, j, k + 1, Idx::PhiReal) + state_fab(i, j, k - 1, Idx::PhiReal) - 2*state_fab(i, j, k, Idx::PhiReal))/std::pow(dx[2], 2);

        amrex::Real dDDPhiImaginary_LL_00 = (state_fab(i + 1, j, k, Idx::PhiImaginary) + state_fab(i - 1, j, k, Idx::PhiImaginary) - 2*state_fab(i, j, k, Idx::PhiImaginary))/std::pow(dx[0], 2);
        amrex::Real dDDPhiImaginary_LL_01 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::PhiImaginary) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::PhiImaginary) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::PhiImaginary) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::PhiImaginary))/(dx[0]*dx[1]);
        amrex::Real dDDPhiImaginary_LL_02 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::PhiImaginary) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::PhiImaginary) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::PhiImaginary) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::PhiImaginary))/(dx[0]*dx[2]);
        amrex::Real dDDPhiImaginary_LL_11 = (state_fab(i, j + 1, k, Idx::PhiImaginary) + state_fab(i, j - 1, k, Idx::PhiImaginary) - 2*state_fab(i, j, k, Idx::PhiImaginary))/std::pow(dx[1], 2);
        amrex::Real dDDPhiImaginary_LL_12 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::PhiImaginary) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::PhiImaginary) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::PhiImaginary) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::PhiImaginary))/(dx[1]*dx[2]);
        amrex::Real dDDPhiImaginary_LL_22 = (state_fab(i, j, k + 1, Idx::PhiImaginary) + state_fab(i, j, k - 1, Idx::PhiImaginary) - 2*state_fab(i, j, k, Idx::PhiImaginary))/std::pow(dx[2], 2);

        amrex::Real dDDChi_LL_00 = (state_fab(i + 1, j, k, Idx::Chi) + state_fab(i - 1, j, k, Idx::Chi) - 2*state_fab(i, j, k, Idx::Chi))/std::pow(dx[0], 2);
        amrex::Real dDDChi_LL_01 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::Chi) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::Chi) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::Chi) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::Chi))/(dx[0]*dx[1]);
        amrex::Real dDDChi_LL_02 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::Chi) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::Chi) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::Chi) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::Chi))/(dx[0]*dx[2]);
        amrex::Real dDDChi_LL_11 = (state_fab(i, j + 1, k, Idx::Chi) + state_fab(i, j - 1, k, Idx::Chi) - 2*state_fab(i, j, k, Idx::Chi))/std::pow(dx[1], 2);
        amrex::Real dDDChi_LL_12 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::Chi) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::Chi) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::Chi) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::Chi))/(dx[1]*dx[2]);
        amrex::Real dDDChi_LL_22 = (state_fab(i, j, k + 1, Idx::Chi) + state_fab(i, j, k - 1, Idx::Chi) - 2*state_fab(i, j, k, Idx::Chi))/std::pow(dx[2], 2);

        amrex::Real dDDChi0_LL_00 = (state_fab(i + 1, j, k, Idx::Chi0) + state_fab(i - 1, j, k, Idx::Chi0) - 2*state_fab(i, j, k, Idx::Chi0))/std::pow(dx[0], 2);
        amrex::Real dDDChi0_LL_01 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::Chi0) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::Chi0) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::Chi0) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::Chi0))/(dx[0]*dx[1]);
        amrex::Real dDDChi0_LL_02 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::Chi0) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::Chi0) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::Chi0) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::Chi0))/(dx[0]*dx[2]);
        amrex::Real dDDChi0_LL_11 = (state_fab(i, j + 1, k, Idx::Chi0) + state_fab(i, j - 1, k, Idx::Chi0) - 2*state_fab(i, j, k, Idx::Chi0))/std::pow(dx[1], 2);
        amrex::Real dDDChi0_LL_12 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::Chi0) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::Chi0) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::Chi0) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::Chi0))/(dx[1]*dx[2]);
        amrex::Real dDDChi0_LL_22 = (state_fab(i, j, k + 1, Idx::Chi0) + state_fab(i, j, k - 1, Idx::Chi0) - 2*state_fab(i, j, k, Idx::Chi0))/std::pow(dx[2], 2);

        amrex::Real dDDChi1_LL_00 = (state_fab(i + 1, j, k, Idx::Chi1) + state_fab(i - 1, j, k, Idx::Chi1) - 2*state_fab(i, j, k, Idx::Chi1))/std::pow(dx[0], 2);
        amrex::Real dDDChi1_LL_01 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::Chi1) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::Chi1) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::Chi1) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::Chi1))/(dx[0]*dx[1]);
        amrex::Real dDDChi1_LL_02 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::Chi1) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::Chi1) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::Chi1) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::Chi1))/(dx[0]*dx[2]);
        amrex::Real dDDChi1_LL_11 = (state_fab(i, j + 1, k, Idx::Chi1) + state_fab(i, j - 1, k, Idx::Chi1) - 2*state_fab(i, j, k, Idx::Chi1))/std::pow(dx[1], 2);
        amrex::Real dDDChi1_LL_12 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::Chi1) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::Chi1) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::Chi1) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::Chi1))/(dx[1]*dx[2]);
        amrex::Real dDDChi1_LL_22 = (state_fab(i, j, k + 1, Idx::Chi1) + state_fab(i, j, k - 1, Idx::Chi1) - 2*state_fab(i, j, k, Idx::Chi1))/std::pow(dx[2], 2);

        amrex::Real dDDA_LLL_000 = (state_fab(i + 1, j, k, Idx::A_L_0) + state_fab(i - 1, j, k, Idx::A_L_0) - 2*state_fab(i, j, k, Idx::A_L_0))/std::pow(dx[0], 2);
        amrex::Real dDDA_LLL_001 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::A_L_0) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::A_L_0) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::A_L_0) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::A_L_0))/(dx[0]*dx[1]);
        amrex::Real dDDA_LLL_002 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::A_L_0) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::A_L_0) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::A_L_0) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::A_L_0))/(dx[0]*dx[2]);
        amrex::Real dDDA_LLL_011 = (state_fab(i, j + 1, k, Idx::A_L_0) + state_fab(i, j - 1, k, Idx::A_L_0) - 2*state_fab(i, j, k, Idx::A_L_0))/std::pow(dx[1], 2);
        amrex::Real dDDA_LLL_012 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::A_L_0) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::A_L_0) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::A_L_0) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::A_L_0))/(dx[1]*dx[2]);
        amrex::Real dDDA_LLL_022 = (state_fab(i, j, k + 1, Idx::A_L_0) + state_fab(i, j, k - 1, Idx::A_L_0) - 2*state_fab(i, j, k, Idx::A_L_0))/std::pow(dx[2], 2);
        amrex::Real dDDA_LLL_100 = (state_fab(i + 1, j, k, Idx::A_L_1) + state_fab(i - 1, j, k, Idx::A_L_1) - 2*state_fab(i, j, k, Idx::A_L_1))/std::pow(dx[0], 2);
        amrex::Real dDDA_LLL_101 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::A_L_1) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::A_L_1) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::A_L_1) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::A_L_1))/(dx[0]*dx[1]);
        amrex::Real dDDA_LLL_102 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::A_L_1) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::A_L_1) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::A_L_1) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::A_L_1))/(dx[0]*dx[2]);
        amrex::Real dDDA_LLL_111 = (state_fab(i, j + 1, k, Idx::A_L_1) + state_fab(i, j - 1, k, Idx::A_L_1) - 2*state_fab(i, j, k, Idx::A_L_1))/std::pow(dx[1], 2);
        amrex::Real dDDA_LLL_112 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::A_L_1) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::A_L_1) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::A_L_1) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::A_L_1))/(dx[1]*dx[2]);
        amrex::Real dDDA_LLL_122 = (state_fab(i, j, k + 1, Idx::A_L_1) + state_fab(i, j, k - 1, Idx::A_L_1) - 2*state_fab(i, j, k, Idx::A_L_1))/std::pow(dx[2], 2);
        amrex::Real dDDA_LLL_200 = (state_fab(i + 1, j, k, Idx::A_L_2) + state_fab(i - 1, j, k, Idx::A_L_2) - 2*state_fab(i, j, k, Idx::A_L_2))/std::pow(dx[0], 2);
        amrex::Real dDDA_LLL_201 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::A_L_2) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::A_L_2) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::A_L_2) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::A_L_2))/(dx[0]*dx[1]);
        amrex::Real dDDA_LLL_202 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::A_L_2) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::A_L_2) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::A_L_2) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::A_L_2))/(dx[0]*dx[2]);
        amrex::Real dDDA_LLL_211 = (state_fab(i, j + 1, k, Idx::A_L_2) + state_fab(i, j - 1, k, Idx::A_L_2) - 2*state_fab(i, j, k, Idx::A_L_2))/std::pow(dx[1], 2);
        amrex::Real dDDA_LLL_212 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::A_L_2) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::A_L_2) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::A_L_2) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::A_L_2))/(dx[1]*dx[2]);
        amrex::Real dDDA_LLL_222 = (state_fab(i, j, k + 1, Idx::A_L_2) + state_fab(i, j, k - 1, Idx::A_L_2) - 2*state_fab(i, j, k, Idx::A_L_2))/std::pow(dx[2], 2);

        amrex::Real dDDAT_LL_00 = (state_fab(i + 1, j, k, Idx::AT) + state_fab(i - 1, j, k, Idx::AT) - 2*state_fab(i, j, k, Idx::AT))/std::pow(dx[0], 2);
        amrex::Real dDDAT_LL_01 = ((1.0/4.0)*state_fab(i + 1, j + 1, k, Idx::AT) - 1.0/4.0*state_fab(i + 1, j - 1, k, Idx::AT) - 1.0/4.0*state_fab(i - 1, j + 1, k, Idx::AT) + (1.0/4.0)*state_fab(i - 1, j - 1, k, Idx::AT))/(dx[0]*dx[1]);
        amrex::Real dDDAT_LL_02 = ((1.0/4.0)*state_fab(i + 1, j, k + 1, Idx::AT) - 1.0/4.0*state_fab(i + 1, j, k - 1, Idx::AT) - 1.0/4.0*state_fab(i - 1, j, k + 1, Idx::AT) + (1.0/4.0)*state_fab(i - 1, j, k - 1, Idx::AT))/(dx[0]*dx[2]);
        amrex::Real dDDAT_LL_11 = (state_fab(i, j + 1, k, Idx::AT) + state_fab(i, j - 1, k, Idx::AT) - 2*state_fab(i, j, k, Idx::AT))/std::pow(dx[1], 2);
        amrex::Real dDDAT_LL_12 = ((1.0/4.0)*state_fab(i, j + 1, k + 1, Idx::AT) - 1.0/4.0*state_fab(i, j + 1, k - 1, Idx::AT) - 1.0/4.0*state_fab(i, j - 1, k + 1, Idx::AT) + (1.0/4.0)*state_fab(i, j - 1, k - 1, Idx::AT))/(dx[1]*dx[2]);
        amrex::Real dDDAT_LL_22 = (state_fab(i, j, k + 1, Idx::AT) + state_fab(i, j, k - 1, Idx::AT) - 2*state_fab(i, j, k, Idx::AT))/std::pow(dx[2], 2);

        amrex::Real LapPhiReal = dDDPhiReal_LL_00 + dDDPhiReal_LL_11 + dDDPhiReal_LL_22;

        amrex::Real LapPhiImaginary = dDDPhiImaginary_LL_00 + dDDPhiImaginary_LL_11 + dDDPhiImaginary_LL_22;

        amrex::Real LapChi = dDDChi_LL_00 + dDDChi_LL_11 + dDDChi_LL_22;
        amrex::Real LapChi0 = dDDChi0_LL_00 + dDDChi0_LL_11 + dDDChi0_LL_22;
        amrex::Real LapChi1 = dDDChi1_LL_00 + dDDChi1_LL_11 + dDDChi1_LL_22;

        amrex::Real LapA_L_0 = dDDA_LLL_000 + dDDA_LLL_011 + dDDA_LLL_022;
        amrex::Real LapA_L_1 = dDDA_LLL_100 + dDDA_LLL_111 + dDDA_LLL_122;
        amrex::Real LapA_L_2 = dDDA_LLL_200 + dDDA_LLL_211 + dDDA_LLL_222;

        amrex::Real LapAT = dDDAT_LL_00 + dDDAT_LL_11 + dDDAT_LL_22;

        amrex::Real j_L_0 = -PhiImaginary*dDPhiReal_L_0 + PhiReal*dDPhiImaginary_L_0;
        amrex::Real j_L_1 = -PhiImaginary*dDPhiReal_L_1 + PhiReal*dDPhiImaginary_L_1;
        amrex::Real j_L_2 = -PhiImaginary*dDPhiReal_L_2 + PhiReal*dDPhiImaginary_L_2;

        amrex::Real jT = -PhiImaginary*PiPhiReal + PhiReal*PiPhiImaginary;

        rhs_fab(i, j, k, Idx::PhiReal) = PiPhiReal;

        rhs_fab(i, j, k, Idx::PhiImaginary) = PiPhiImaginary;

        rhs_fab(i, j, k, Idx::Chi) = PiChi;
        rhs_fab(i, j, k, Idx::Chi0) = PiChi0;
        rhs_fab(i, j, k, Idx::Chi1) = PiChi1;

        rhs_fab(i, j, k, Idx::A_L_0) = PiA_L_0;
        rhs_fab(i, j, k, Idx::A_L_1) = PiA_L_1;
        rhs_fab(i, j, k, Idx::A_L_2) = PiA_L_2;

        rhs_fab(i, j, k, Idx::AT) = PiAT;

        rhs_fab(i, j, k, Idx::PiPhiReal) = std::pow(AT, 2)*PhiReal*std::pow(e, 2) + 2*AT*PiPhiImaginary*e - std::pow(A_L_0, 2)*PhiReal*std::pow(e, 2) - 2*A_L_0*dDPhiImaginary_L_0*e - std::pow(A_L_1, 2)*PhiReal*std::pow(e, 2) - 2*A_L_1*dDPhiImaginary_L_1*e - std::pow(A_L_2, 2)*PhiReal*std::pow(e, 2) - 2*A_L_2*dDPhiImaginary_L_2*e + LapPhiReal - 2*PhiReal*lam*(std::pow(Chi, 2)+std::pow(Chi0, 2)+std::pow(Chi1, 2) + std::pow(PhiImaginary, 2) + std::pow(PhiReal, 2) - std::pow(v, 2));

        rhs_fab(i, j, k, Idx::PiPhiImaginary) = std::pow(AT, 2)*PhiImaginary*std::pow(e, 2) - 2*AT*PiPhiReal*e - std::pow(A_L_0, 2)*PhiImaginary*std::pow(e, 2) + 2*A_L_0*dDPhiReal_L_0*e - std::pow(A_L_1, 2)*PhiImaginary*std::pow(e, 2) + 2*A_L_1*dDPhiReal_L_1*e - std::pow(A_L_2, 2)*PhiImaginary*std::pow(e, 2) + 2*A_L_2*dDPhiReal_L_2*e + LapPhiImaginary - 2*PhiImaginary*lam*(std::pow(Chi, 2)+std::pow(Chi0, 2)+std::pow(Chi1, 2) + std::pow(PhiImaginary, 2) + std::pow(PhiReal, 2) - std::pow(v, 2));

        rhs_fab(i, j, k, Idx::PiChi) = -Chi*std::pow(Omega, 2) - 2*Chi*lam*(std::pow(Chi, 2)+std::pow(Chi0, 2)+std::pow(Chi1, 2) + std::pow(PhiImaginary, 2) + std::pow(PhiReal, 2) - std::pow(v, 2)) + LapChi;
        
        rhs_fab(i, j, k, Idx::PiChi0) = -Chi0*std::pow(Omega, 2) - 2*Chi0*lam*(std::pow(Chi, 2)+std::pow(Chi0, 2)+std::pow(Chi1, 2) + std::pow(PhiImaginary, 2) + std::pow(PhiReal, 2) - std::pow(v, 2)) + LapChi0;
        
        rhs_fab(i, j, k, Idx::PiChi1) = -Chi1*std::pow(Omega, 2) - 2*Chi1*lam*(std::pow(Chi, 2)+std::pow(Chi0, 2)+std::pow(Chi1, 2) + std::pow(PhiImaginary, 2) + std::pow(PhiReal, 2) - std::pow(v, 2)) + LapChi1;

        rhs_fab(i, j, k, Idx::PiA_L_0) = -2*A_L_0*std::pow(e, 2)*(std::pow(PhiImaginary, 2) + std::pow(PhiReal, 2)) + LapA_L_0 - 2*e*j_L_0;
        rhs_fab(i, j, k, Idx::PiA_L_1) = -2*A_L_1*std::pow(e, 2)*(std::pow(PhiImaginary, 2) + std::pow(PhiReal, 2)) + LapA_L_1 - 2*e*j_L_1;
        rhs_fab(i, j, k, Idx::PiA_L_2) = -2*A_L_2*std::pow(e, 2)*(std::pow(PhiImaginary, 2) + std::pow(PhiReal, 2)) + LapA_L_2 - 2*e*j_L_2;

        rhs_fab(i, j, k, Idx::PiAT) = -2*AT*std::pow(e, 2)*(std::pow(PhiImaginary, 2) + std::pow(PhiReal, 2)) + LapAT - 2*e*jT;
        
}
#endif
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

        amrex::Real du0 = ((2.0/3.0)*state_fab(i + 1, j, k, Idx::u) - 1.0/12.0*state_fab(i + 2, j, k, Idx::u) - 2.0/3.0*state_fab(i - 1, j, k, Idx::u) + (1.0/12.0)*state_fab(i - 2, j, k, Idx::u))/dx[0];
        amrex::Real du1 = ((2.0/3.0)*state_fab(i, j + 1, k, Idx::u) - 1.0/12.0*state_fab(i, j + 2, k, Idx::u) - 2.0/3.0*state_fab(i, j - 1, k, Idx::u) + (1.0/12.0)*state_fab(i, j - 2, k, Idx::u))/dx[1];
        amrex::Real du2 = ((2.0/3.0)*state_fab(i, j, k + 1, Idx::u) - 1.0/12.0*state_fab(i, j, k + 2, Idx::u) - 2.0/3.0*state_fab(i, j, k - 1, Idx::u) + (1.0/12.0)*state_fab(i, j, k - 2, Idx::u))/dx[2];
        amrex::Real ddu00 = ((4.0/3.0)*state_fab(i + 1, j, k, Idx::u) - 1.0/12.0*state_fab(i + 2, j, k, Idx::u) + (4.0/3.0)*state_fab(i - 1, j, k, Idx::u) - 1.0/12.0*state_fab(i - 2, j, k, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[0], 2);
        amrex::Real ddu01 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::u) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::u) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::u) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::u) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::u) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::u) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::u) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::u) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::u) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::u) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::u) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::u) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::u) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::u) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::u) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::u))/(dx[0]*dx[1]);
        amrex::Real ddu02 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::u) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::u) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::u) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::u) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::u) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::u) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::u) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::u) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::u) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::u) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::u) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::u) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::u) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::u) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::u) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::u))/(dx[0]*dx[2]);
        amrex::Real ddu10 = ((4.0/9.0)*state_fab(i + 1, j + 1, k, Idx::u) - 1.0/18.0*state_fab(i + 1, j + 2, k, Idx::u) - 4.0/9.0*state_fab(i + 1, j - 1, k, Idx::u) + (1.0/18.0)*state_fab(i + 1, j - 2, k, Idx::u) - 1.0/18.0*state_fab(i + 2, j + 1, k, Idx::u) + (1.0/144.0)*state_fab(i + 2, j + 2, k, Idx::u) + (1.0/18.0)*state_fab(i + 2, j - 1, k, Idx::u) - 1.0/144.0*state_fab(i + 2, j - 2, k, Idx::u) - 4.0/9.0*state_fab(i - 1, j + 1, k, Idx::u) + (1.0/18.0)*state_fab(i - 1, j + 2, k, Idx::u) + (4.0/9.0)*state_fab(i - 1, j - 1, k, Idx::u) - 1.0/18.0*state_fab(i - 1, j - 2, k, Idx::u) + (1.0/18.0)*state_fab(i - 2, j + 1, k, Idx::u) - 1.0/144.0*state_fab(i - 2, j + 2, k, Idx::u) - 1.0/18.0*state_fab(i - 2, j - 1, k, Idx::u) + (1.0/144.0)*state_fab(i - 2, j - 2, k, Idx::u))/(dx[0]*dx[1]);
        amrex::Real ddu11 = ((4.0/3.0)*state_fab(i, j + 1, k, Idx::u) - 1.0/12.0*state_fab(i, j + 2, k, Idx::u) + (4.0/3.0)*state_fab(i, j - 1, k, Idx::u) - 1.0/12.0*state_fab(i, j - 2, k, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[1], 2);
        amrex::Real ddu12 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::u) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::u) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::u) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::u) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::u) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::u) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::u) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::u) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::u) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::u) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::u) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::u) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::u) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::u) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::u) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::u))/(dx[1]*dx[2]);
        amrex::Real ddu20 = ((4.0/9.0)*state_fab(i + 1, j, k + 1, Idx::u) - 1.0/18.0*state_fab(i + 1, j, k + 2, Idx::u) - 4.0/9.0*state_fab(i + 1, j, k - 1, Idx::u) + (1.0/18.0)*state_fab(i + 1, j, k - 2, Idx::u) - 1.0/18.0*state_fab(i + 2, j, k + 1, Idx::u) + (1.0/144.0)*state_fab(i + 2, j, k + 2, Idx::u) + (1.0/18.0)*state_fab(i + 2, j, k - 1, Idx::u) - 1.0/144.0*state_fab(i + 2, j, k - 2, Idx::u) - 4.0/9.0*state_fab(i - 1, j, k + 1, Idx::u) + (1.0/18.0)*state_fab(i - 1, j, k + 2, Idx::u) + (4.0/9.0)*state_fab(i - 1, j, k - 1, Idx::u) - 1.0/18.0*state_fab(i - 1, j, k - 2, Idx::u) + (1.0/18.0)*state_fab(i - 2, j, k + 1, Idx::u) - 1.0/144.0*state_fab(i - 2, j, k + 2, Idx::u) - 1.0/18.0*state_fab(i - 2, j, k - 1, Idx::u) + (1.0/144.0)*state_fab(i - 2, j, k - 2, Idx::u))/(dx[0]*dx[2]);
        amrex::Real ddu21 = ((4.0/9.0)*state_fab(i, j + 1, k + 1, Idx::u) - 1.0/18.0*state_fab(i, j + 1, k + 2, Idx::u) - 4.0/9.0*state_fab(i, j + 1, k - 1, Idx::u) + (1.0/18.0)*state_fab(i, j + 1, k - 2, Idx::u) - 1.0/18.0*state_fab(i, j + 2, k + 1, Idx::u) + (1.0/144.0)*state_fab(i, j + 2, k + 2, Idx::u) + (1.0/18.0)*state_fab(i, j + 2, k - 1, Idx::u) - 1.0/144.0*state_fab(i, j + 2, k - 2, Idx::u) - 4.0/9.0*state_fab(i, j - 1, k + 1, Idx::u) + (1.0/18.0)*state_fab(i, j - 1, k + 2, Idx::u) + (4.0/9.0)*state_fab(i, j - 1, k - 1, Idx::u) - 1.0/18.0*state_fab(i, j - 1, k - 2, Idx::u) + (1.0/18.0)*state_fab(i, j - 2, k + 1, Idx::u) - 1.0/144.0*state_fab(i, j - 2, k + 2, Idx::u) - 1.0/18.0*state_fab(i, j - 2, k - 1, Idx::u) + (1.0/144.0)*state_fab(i, j - 2, k - 2, Idx::u))/(dx[1]*dx[2]);
        amrex::Real ddu22 = ((4.0/3.0)*state_fab(i, j, k + 1, Idx::u) - 1.0/12.0*state_fab(i, j, k + 2, Idx::u) + (4.0/3.0)*state_fab(i, j, k - 1, Idx::u) - 1.0/12.0*state_fab(i, j, k - 2, Idx::u) - 5.0/2.0*state_fab(i, j, k, Idx::u))/std::pow(dx[2], 2);
        amrex::Real KOu = 0.0015625000000000001*(15*state_fab(i, j, k + 1, Idx::u) - 6*state_fab(i, j, k + 2, Idx::u) + state_fab(i, j, k + 3, Idx::u) + 15*state_fab(i, j, k - 1, Idx::u) - 6*state_fab(i, j, k - 2, Idx::u) + state_fab(i, j, k - 3, Idx::u) - 20*state_fab(i, j, k, Idx::u))/dx[2] + 0.0015625000000000001*(15*state_fab(i, j + 1, k, Idx::u) - 6*state_fab(i, j + 2, k, Idx::u) + state_fab(i, j + 3, k, Idx::u) + 15*state_fab(i, j - 1, k, Idx::u) - 6*state_fab(i, j - 2, k, Idx::u) + state_fab(i, j - 3, k, Idx::u) - 20*state_fab(i, j, k, Idx::u))/dx[1] + 0.0015625000000000001*(15*state_fab(i + 1, j, k, Idx::u) - 6*state_fab(i + 2, j, k, Idx::u) + state_fab(i + 3, j, k, Idx::u) + 15*state_fab(i - 1, j, k, Idx::u) - 6*state_fab(i - 2, j, k, Idx::u) + state_fab(i - 3, j, k, Idx::u) - 20*state_fab(i, j, k, Idx::u))/dx[0];

        amrex::Real r = std::sqrt(std::pow(x0, 2) + std::pow(x1, 2) + std::pow(x2, 2));

        amrex::Real c10 = 4.891; //1.168642873;
        amrex::Real c11 = 0.0;
        amrex::Real c12 = 0.0;
        amrex::Real c20 = -4.891; //-1.168642873;
        amrex::Real c21 = 0.0;
        amrex::Real c22 = 0.0;
        amrex::Real M1 = 0.4856; //0.45300000000000001;

        amrex::Real M2 = 0.4856; //0.45300000000000001;

        amrex::Real P10 = 0.0;
        amrex::Real P11 = 0.0969; //0.33319174979999999;
        amrex::Real P12 = 0.0;
        amrex::Real P20 = 0.0;
        amrex::Real P21 = -0.0969; //-0.33319174979999999;
        amrex::Real P22 = 0.0;
        amrex::Real S10 = 0;
        amrex::Real S11 = 0;
        amrex::Real S12 = 0;
        amrex::Real S20 = 0;
        amrex::Real S21 = 0;
        amrex::Real S22 = 0;
        amrex::Real xc10 = -c10 + x0;
        amrex::Real xc11 = -c11 + x1;
        amrex::Real xc12 = -c12 + x2;
        amrex::Real xc20 = -c20 + x0;
        amrex::Real xc21 = -c21 + x1;
        amrex::Real xc22 = -c22 + x2;
        amrex::Real rc1 = std::max(0.0100000000000000, std::sqrt(std::pow(xc10, 2) + std::pow(xc11, 2) + std::pow(xc12, 2)));

        amrex::Real rc2 = std::max(0.0100000000000000, std::sqrt(std::pow(xc20, 2) + std::pow(xc21, 2) + std::pow(xc22, 2)));

        amrex::Real n10 = xc10/rc1;
        amrex::Real n11 = xc11/rc1;
        amrex::Real n12 = xc12/rc1;
        amrex::Real n20 = xc20/rc2;
        amrex::Real n21 = xc21/rc2;
        amrex::Real n22 = xc22/rc2;
        amrex::Real Abar1UU00 = -3.0/2.0*P10*n10*(1 - std::pow(n10, 2))/std::pow(rc1, 2) + 3*P10*n10/std::pow(rc1, 2) - 3.0/2.0*P11*n11*(1 - std::pow(n10, 2))/std::pow(rc1, 2) - 3.0/2.0*P12*n12*(1 - std::pow(n10, 2))/std::pow(rc1, 2) + 6*S11*n10*n12/std::pow(rc1, 3) - 6*S12*n10*n11/std::pow(rc1, 3);
        amrex::Real Abar1UU01 = (3.0/2.0)*P10*std::pow(n10, 2)*n11/std::pow(rc1, 2) + (3.0/2.0)*P11*n10*std::pow(n11, 2)/std::pow(rc1, 2) + (3.0/2.0)*P12*n10*n11*n12/std::pow(rc1, 2) - 3*S10*n10*n12/std::pow(rc1, 3) + 3*S11*n11*n12/std::pow(rc1, 3) + 3*S12*std::pow(n10, 2)/std::pow(rc1, 3) - 3*S12*std::pow(n11, 2)/std::pow(rc1, 3) + (3.0/2.0)*(P10*n11 + P11*n10)/std::pow(rc1, 2);
        amrex::Real Abar1UU02 = (3.0/2.0)*P10*std::pow(n10, 2)*n12/std::pow(rc1, 2) + (3.0/2.0)*P11*n10*n11*n12/std::pow(rc1, 2) + (3.0/2.0)*P12*n10*std::pow(n12, 2)/std::pow(rc1, 2) + 3*S10*n10*n11/std::pow(rc1, 3) - 3*S11*std::pow(n10, 2)/std::pow(rc1, 3) + 3*S11*std::pow(n12, 2)/std::pow(rc1, 3) - 3*S12*n11*n12/std::pow(rc1, 3) + (3.0/2.0)*(P10*n12 + P12*n10)/std::pow(rc1, 2);
        amrex::Real Abar1UU10 = (3.0/2.0)*P10*std::pow(n10, 2)*n11/std::pow(rc1, 2) + (3.0/2.0)*P11*n10*std::pow(n11, 2)/std::pow(rc1, 2) + (3.0/2.0)*P12*n10*n11*n12/std::pow(rc1, 2) - 3*S10*n10*n12/std::pow(rc1, 3) + 3*S11*n11*n12/std::pow(rc1, 3) + 3*S12*std::pow(n10, 2)/std::pow(rc1, 3) - 3*S12*std::pow(n11, 2)/std::pow(rc1, 3) + (3.0/2.0)*(P10*n11 + P11*n10)/std::pow(rc1, 2);
        amrex::Real Abar1UU11 = -3.0/2.0*P10*n10*(1 - std::pow(n11, 2))/std::pow(rc1, 2) - 3.0/2.0*P11*n11*(1 - std::pow(n11, 2))/std::pow(rc1, 2) + 3*P11*n11/std::pow(rc1, 2) - 3.0/2.0*P12*n12*(1 - std::pow(n11, 2))/std::pow(rc1, 2) - 6*S10*n11*n12/std::pow(rc1, 3) + 6*S12*n10*n11/std::pow(rc1, 3);
        amrex::Real Abar1UU12 = (3.0/2.0)*P10*n10*n11*n12/std::pow(rc1, 2) + (3.0/2.0)*P11*std::pow(n11, 2)*n12/std::pow(rc1, 2) + (3.0/2.0)*P12*n11*std::pow(n12, 2)/std::pow(rc1, 2) + 3*S10*std::pow(n11, 2)/std::pow(rc1, 3) - 3*S10*std::pow(n12, 2)/std::pow(rc1, 3) - 3*S11*n10*n11/std::pow(rc1, 3) + 3*S12*n10*n12/std::pow(rc1, 3) + (3.0/2.0)*(P11*n12 + P12*n11)/std::pow(rc1, 2);
        amrex::Real Abar1UU20 = (3.0/2.0)*P10*std::pow(n10, 2)*n12/std::pow(rc1, 2) + (3.0/2.0)*P11*n10*n11*n12/std::pow(rc1, 2) + (3.0/2.0)*P12*n10*std::pow(n12, 2)/std::pow(rc1, 2) + 3*S10*n10*n11/std::pow(rc1, 3) - 3*S11*std::pow(n10, 2)/std::pow(rc1, 3) + 3*S11*std::pow(n12, 2)/std::pow(rc1, 3) - 3*S12*n11*n12/std::pow(rc1, 3) + (3.0/2.0)*(P10*n12 + P12*n10)/std::pow(rc1, 2);
        amrex::Real Abar1UU21 = (3.0/2.0)*P10*n10*n11*n12/std::pow(rc1, 2) + (3.0/2.0)*P11*std::pow(n11, 2)*n12/std::pow(rc1, 2) + (3.0/2.0)*P12*n11*std::pow(n12, 2)/std::pow(rc1, 2) + 3*S10*std::pow(n11, 2)/std::pow(rc1, 3) - 3*S10*std::pow(n12, 2)/std::pow(rc1, 3) - 3*S11*n10*n11/std::pow(rc1, 3) + 3*S12*n10*n12/std::pow(rc1, 3) + (3.0/2.0)*(P11*n12 + P12*n11)/std::pow(rc1, 2);
        amrex::Real Abar1UU22 = -3.0/2.0*P10*n10*(1 - std::pow(n12, 2))/std::pow(rc1, 2) - 3.0/2.0*P11*n11*(1 - std::pow(n12, 2))/std::pow(rc1, 2) - 3.0/2.0*P12*n12*(1 - std::pow(n12, 2))/std::pow(rc1, 2) + 3*P12*n12/std::pow(rc1, 2) + 6*S10*n11*n12/std::pow(rc1, 3) - 6*S11*n10*n12/std::pow(rc1, 3);
        amrex::Real Abar2UU00 = -3.0/2.0*P20*n20*(1 - std::pow(n20, 2))/std::pow(rc2, 2) + 3*P20*n20/std::pow(rc2, 2) - 3.0/2.0*P21*n21*(1 - std::pow(n20, 2))/std::pow(rc2, 2) - 3.0/2.0*P22*n22*(1 - std::pow(n20, 2))/std::pow(rc2, 2) + 6*S21*n20*n22/std::pow(rc2, 3) - 6*S22*n20*n21/std::pow(rc2, 3);
        amrex::Real Abar2UU01 = (3.0/2.0)*P20*std::pow(n20, 2)*n21/std::pow(rc2, 2) + (3.0/2.0)*P21*n20*std::pow(n21, 2)/std::pow(rc2, 2) + (3.0/2.0)*P22*n20*n21*n22/std::pow(rc2, 2) - 3*S20*n20*n22/std::pow(rc2, 3) + 3*S21*n21*n22/std::pow(rc2, 3) + 3*S22*std::pow(n20, 2)/std::pow(rc2, 3) - 3*S22*std::pow(n21, 2)/std::pow(rc2, 3) + (3.0/2.0)*(P20*n21 + P21*n20)/std::pow(rc2, 2);
        amrex::Real Abar2UU02 = (3.0/2.0)*P20*std::pow(n20, 2)*n22/std::pow(rc2, 2) + (3.0/2.0)*P21*n20*n21*n22/std::pow(rc2, 2) + (3.0/2.0)*P22*n20*std::pow(n22, 2)/std::pow(rc2, 2) + 3*S20*n20*n21/std::pow(rc2, 3) - 3*S21*std::pow(n20, 2)/std::pow(rc2, 3) + 3*S21*std::pow(n22, 2)/std::pow(rc2, 3) - 3*S22*n21*n22/std::pow(rc2, 3) + (3.0/2.0)*(P20*n22 + P22*n20)/std::pow(rc2, 2);
        amrex::Real Abar2UU10 = (3.0/2.0)*P20*std::pow(n20, 2)*n21/std::pow(rc2, 2) + (3.0/2.0)*P21*n20*std::pow(n21, 2)/std::pow(rc2, 2) + (3.0/2.0)*P22*n20*n21*n22/std::pow(rc2, 2) - 3*S20*n20*n22/std::pow(rc2, 3) + 3*S21*n21*n22/std::pow(rc2, 3) + 3*S22*std::pow(n20, 2)/std::pow(rc2, 3) - 3*S22*std::pow(n21, 2)/std::pow(rc2, 3) + (3.0/2.0)*(P20*n21 + P21*n20)/std::pow(rc2, 2);
        amrex::Real Abar2UU11 = -3.0/2.0*P20*n20*(1 - std::pow(n21, 2))/std::pow(rc2, 2) - 3.0/2.0*P21*n21*(1 - std::pow(n21, 2))/std::pow(rc2, 2) + 3*P21*n21/std::pow(rc2, 2) - 3.0/2.0*P22*n22*(1 - std::pow(n21, 2))/std::pow(rc2, 2) - 6*S20*n21*n22/std::pow(rc2, 3) + 6*S22*n20*n21/std::pow(rc2, 3);
        amrex::Real Abar2UU12 = (3.0/2.0)*P20*n20*n21*n22/std::pow(rc2, 2) + (3.0/2.0)*P21*std::pow(n21, 2)*n22/std::pow(rc2, 2) + (3.0/2.0)*P22*n21*std::pow(n22, 2)/std::pow(rc2, 2) + 3*S20*std::pow(n21, 2)/std::pow(rc2, 3) - 3*S20*std::pow(n22, 2)/std::pow(rc2, 3) - 3*S21*n20*n21/std::pow(rc2, 3) + 3*S22*n20*n22/std::pow(rc2, 3) + (3.0/2.0)*(P21*n22 + P22*n21)/std::pow(rc2, 2);
        amrex::Real Abar2UU20 = (3.0/2.0)*P20*std::pow(n20, 2)*n22/std::pow(rc2, 2) + (3.0/2.0)*P21*n20*n21*n22/std::pow(rc2, 2) + (3.0/2.0)*P22*n20*std::pow(n22, 2)/std::pow(rc2, 2) + 3*S20*n20*n21/std::pow(rc2, 3) - 3*S21*std::pow(n20, 2)/std::pow(rc2, 3) + 3*S21*std::pow(n22, 2)/std::pow(rc2, 3) - 3*S22*n21*n22/std::pow(rc2, 3) + (3.0/2.0)*(P20*n22 + P22*n20)/std::pow(rc2, 2);
        amrex::Real Abar2UU21 = (3.0/2.0)*P20*n20*n21*n22/std::pow(rc2, 2) + (3.0/2.0)*P21*std::pow(n21, 2)*n22/std::pow(rc2, 2) + (3.0/2.0)*P22*n21*std::pow(n22, 2)/std::pow(rc2, 2) + 3*S20*std::pow(n21, 2)/std::pow(rc2, 3) - 3*S20*std::pow(n22, 2)/std::pow(rc2, 3) - 3*S21*n20*n21/std::pow(rc2, 3) + 3*S22*n20*n22/std::pow(rc2, 3) + (3.0/2.0)*(P21*n22 + P22*n21)/std::pow(rc2, 2);
        amrex::Real Abar2UU22 = -3.0/2.0*P20*n20*(1 - std::pow(n22, 2))/std::pow(rc2, 2) - 3.0/2.0*P21*n21*(1 - std::pow(n22, 2))/std::pow(rc2, 2) - 3.0/2.0*P22*n22*(1 - std::pow(n22, 2))/std::pow(rc2, 2) + 3*P22*n22/std::pow(rc2, 2) + 6*S20*n21*n22/std::pow(rc2, 3) - 6*S21*n20*n22/std::pow(rc2, 3);
        amrex::Real AbarUU00 = Abar1UU00 + Abar2UU00;
        amrex::Real AbarUU01 = Abar1UU01 + Abar2UU01;
        amrex::Real AbarUU02 = Abar1UU02 + Abar2UU02;
        amrex::Real AbarUU10 = Abar1UU10 + Abar2UU10;
        amrex::Real AbarUU11 = Abar1UU11 + Abar2UU11;
        amrex::Real AbarUU12 = Abar1UU12 + Abar2UU12;
        amrex::Real AbarUU20 = Abar1UU20 + Abar2UU20;
        amrex::Real AbarUU21 = Abar1UU21 + Abar2UU21;
        amrex::Real AbarUU22 = Abar1UU22 + Abar2UU22;
        amrex::Real alpha = 1.0/((1.0/2.0)*M1/rc1 + (1.0/2.0)*M2/rc2);

        amrex::Real beta = 0.125*std::pow(AbarUU00, 2)*std::pow(alpha, 7) + 0.125*std::pow(AbarUU01, 2)*std::pow(alpha, 7) + 0.125*std::pow(AbarUU02, 2)*std::pow(alpha, 7) + 0.125*std::pow(AbarUU10, 2)*std::pow(alpha, 7) + 0.125*std::pow(AbarUU11, 2)*std::pow(alpha, 7) + 0.125*std::pow(AbarUU12, 2)*std::pow(alpha, 7) + 0.125*std::pow(AbarUU20, 2)*std::pow(alpha, 7) + 0.125*std::pow(AbarUU21, 2)*std::pow(alpha, 7) + 0.125*std::pow(AbarUU22, 2)*std::pow(alpha, 7);

        rhs_fab(i, j, k, Idx::u) = KOu + beta/std::pow(alpha*(u + 1) + 1, 7) + ddu00 + ddu11 + ddu22;

}
#endif
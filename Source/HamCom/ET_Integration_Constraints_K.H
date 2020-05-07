#ifndef ET_INTEGRATION_CONSTRAINTS_K_H 
#define ET_INTEGRATION_CONSTRAINTS_K_H 

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_diagnostics(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& diag, 
        amrex::Array4<amrex::Real const> const& state_fab, 
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo(); 

        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(0) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(0) + domain_xlo[2];
        
        amrex::Real M1 = 1.0;
        amrex::Real M2 = 1.0;
        
        amrex::Real r = std::sqrt(std::pow(x0, 2) + std::pow(x1, 2) + std::pow(x2, 2));
        amrex::Real c10 = 3.0;
        amrex::Real c11 = 0.0;
        amrex::Real c12 = 0.0;
        amrex::Real c20 = -3.0;
        amrex::Real c21 = 0.0;
        amrex::Real c22 = 0.0;
        amrex::Real P10 = 0.0;
        amrex::Real P11 = 0.20000000000000001;
        amrex::Real P12 = 0.0;
        amrex::Real P20 = 0.0;
        amrex::Real P21 = -0.20000000000000001;
        amrex::Real P22 = 0.0;
        amrex::Real xc10 = -c10 + x0;
        amrex::Real xc11 = -c11 + x1;
        amrex::Real xc12 = -c12 + x2;
        amrex::Real xc20 = -c20 + x0;
        amrex::Real xc21 = -c21 + x1;
        amrex::Real xc22 = -c22 + x2;
        amrex::Real rc1 = std::sqrt(std::pow(xc10, 2) + std::pow(xc11, 2) + std::pow(xc12, 2));
        amrex::Real rc2 = std::sqrt(std::pow(xc20, 2) + std::pow(xc21, 2) + std::pow(xc22, 2));
        
        amrex::Real eps1;
        amrex::Real eps2;
        
        if(rc1<=0.5)
          eps1 = 0.01;
        else
          eps1 = 0;

        if(rc2<=0.5)
          eps2 = 0.01;
        else
          eps2 = 0;
        
        amrex::Real n10 = xc10/rc1;
        amrex::Real n11 = xc11/rc1;
        amrex::Real n12 = xc12/rc1;
        amrex::Real n20 = xc20/rc2;
        amrex::Real n21 = xc21/rc2;
        amrex::Real n22 = xc22/rc2;
        amrex::Real Atilde1UU00 = -3.0/2.0*P10*n10*(1 - std::pow(n10, 2))/std::pow(eps1 + rc1, 2) + 3*P10*n10/std::pow(eps1 + rc1, 2) - 3.0/2.0*P11*n11*(1 - std::pow(n10, 2))/std::pow(eps1 + rc1, 2) - 3.0/2.0*P12*n12*(1 - std::pow(n10, 2))/std::pow(eps1 + rc1, 2);
        amrex::Real Atilde1UU01 = (3.0/2.0)*P10*std::pow(n10, 2)*n11/std::pow(eps1 + rc1, 2) + (3.0/2.0)*P11*n10*std::pow(n11, 2)/std::pow(eps1 + rc1, 2) + (3.0/2.0)*P12*n10*n11*n12/std::pow(eps1 + rc1, 2) + (3.0/2.0)*(P10*n11 + P11*n10)/std::pow(eps1 + rc1, 2);
        amrex::Real Atilde1UU02 = (3.0/2.0)*P10*std::pow(n10, 2)*n12/std::pow(eps1 + rc1, 2) + (3.0/2.0)*P11*n10*n11*n12/std::pow(eps1 + rc1, 2) + (3.0/2.0)*P12*n10*std::pow(n12, 2)/std::pow(eps1 + rc1, 2) + (3.0/2.0)*(P10*n12 + P12*n10)/std::pow(eps1 + rc1, 2);
        amrex::Real Atilde1UU10 = (3.0/2.0)*P10*std::pow(n10, 2)*n11/std::pow(eps1 + rc1, 2) + (3.0/2.0)*P11*n10*std::pow(n11, 2)/std::pow(eps1 + rc1, 2) + (3.0/2.0)*P12*n10*n11*n12/std::pow(eps1 + rc1, 2) + (3.0/2.0)*(P10*n11 + P11*n10)/std::pow(eps1 + rc1, 2);
        amrex::Real Atilde1UU11 = -3.0/2.0*P10*n10*(1 - std::pow(n11, 2))/std::pow(eps1 + rc1, 2) - 3.0/2.0*P11*n11*(1 - std::pow(n11, 2))/std::pow(eps1 + rc1, 2) + 3*P11*n11/std::pow(eps1 + rc1, 2) - 3.0/2.0*P12*n12*(1 - std::pow(n11, 2))/std::pow(eps1 + rc1, 2);
        amrex::Real Atilde1UU12 = (3.0/2.0)*P10*n10*n11*n12/std::pow(eps1 + rc1, 2) + (3.0/2.0)*P11*std::pow(n11, 2)*n12/std::pow(eps1 + rc1, 2) + (3.0/2.0)*P12*n11*std::pow(n12, 2)/std::pow(eps1 + rc1, 2) + (3.0/2.0)*(P11*n12 + P12*n11)/std::pow(eps1 + rc1, 2);
        amrex::Real Atilde1UU20 = (3.0/2.0)*P10*std::pow(n10, 2)*n12/std::pow(eps1 + rc1, 2) + (3.0/2.0)*P11*n10*n11*n12/std::pow(eps1 + rc1, 2) + (3.0/2.0)*P12*n10*std::pow(n12, 2)/std::pow(eps1 + rc1, 2) + (3.0/2.0)*(P10*n12 + P12*n10)/std::pow(eps1 + rc1, 2);
        amrex::Real Atilde1UU21 = (3.0/2.0)*P10*n10*n11*n12/std::pow(eps1 + rc1, 2) + (3.0/2.0)*P11*std::pow(n11, 2)*n12/std::pow(eps1 + rc1, 2) + (3.0/2.0)*P12*n11*std::pow(n12, 2)/std::pow(eps1 + rc1, 2) + (3.0/2.0)*(P11*n12 + P12*n11)/std::pow(eps1 + rc1, 2);
        amrex::Real Atilde1UU22 = -3.0/2.0*P10*n10*(1 - std::pow(n12, 2))/std::pow(eps1 + rc1, 2) - 3.0/2.0*P11*n11*(1 - std::pow(n12, 2))/std::pow(eps1 + rc1, 2) - 3.0/2.0*P12*n12*(1 - std::pow(n12, 2))/std::pow(eps1 + rc1, 2) + 3*P12*n12/std::pow(eps1 + rc1, 2);
        amrex::Real Atilde2UU00 = -3.0/2.0*P20*n20*(1 - std::pow(n20, 2))/std::pow(eps2 + rc2, 2) + 3*P20*n20/std::pow(eps2 + rc2, 2) - 3.0/2.0*P21*n21*(1 - std::pow(n20, 2))/std::pow(eps2 + rc2, 2) - 3.0/2.0*P22*n22*(1 - std::pow(n20, 2))/std::pow(eps2 + rc2, 2);
        amrex::Real Atilde2UU01 = (3.0/2.0)*P20*std::pow(n20, 2)*n21/std::pow(eps2 + rc2, 2) + (3.0/2.0)*P21*n20*std::pow(n21, 2)/std::pow(eps2 + rc2, 2) + (3.0/2.0)*P22*n20*n21*n22/std::pow(eps2 + rc2, 2) + (3.0/2.0)*(P20*n21 + P21*n20)/std::pow(eps2 + rc2, 2);
        amrex::Real Atilde2UU02 = (3.0/2.0)*P20*std::pow(n20, 2)*n22/std::pow(eps2 + rc2, 2) + (3.0/2.0)*P21*n20*n21*n22/std::pow(eps2 + rc2, 2) + (3.0/2.0)*P22*n20*std::pow(n22, 2)/std::pow(eps2 + rc2, 2) + (3.0/2.0)*(P20*n22 + P22*n20)/std::pow(eps2 + rc2, 2);
        amrex::Real Atilde2UU10 = (3.0/2.0)*P20*std::pow(n20, 2)*n21/std::pow(eps2 + rc2, 2) + (3.0/2.0)*P21*n20*std::pow(n21, 2)/std::pow(eps2 + rc2, 2) + (3.0/2.0)*P22*n20*n21*n22/std::pow(eps2 + rc2, 2) + (3.0/2.0)*(P20*n21 + P21*n20)/std::pow(eps2 + rc2, 2);
        amrex::Real Atilde2UU11 = -3.0/2.0*P20*n20*(1 - std::pow(n21, 2))/std::pow(eps2 + rc2, 2) - 3.0/2.0*P21*n21*(1 - std::pow(n21, 2))/std::pow(eps2 + rc2, 2) + 3*P21*n21/std::pow(eps2 + rc2, 2) - 3.0/2.0*P22*n22*(1 - std::pow(n21, 2))/std::pow(eps2 + rc2, 2);
        amrex::Real Atilde2UU12 = (3.0/2.0)*P20*n20*n21*n22/std::pow(eps2 + rc2, 2) + (3.0/2.0)*P21*std::pow(n21, 2)*n22/std::pow(eps2 + rc2, 2) + (3.0/2.0)*P22*n21*std::pow(n22, 2)/std::pow(eps2 + rc2, 2) + (3.0/2.0)*(P21*n22 + P22*n21)/std::pow(eps2 + rc2, 2);
        amrex::Real Atilde2UU20 = (3.0/2.0)*P20*std::pow(n20, 2)*n22/std::pow(eps2 + rc2, 2) + (3.0/2.0)*P21*n20*n21*n22/std::pow(eps2 + rc2, 2) + (3.0/2.0)*P22*n20*std::pow(n22, 2)/std::pow(eps2 + rc2, 2) + (3.0/2.0)*(P20*n22 + P22*n20)/std::pow(eps2 + rc2, 2);
        amrex::Real Atilde2UU21 = (3.0/2.0)*P20*n20*n21*n22/std::pow(eps2 + rc2, 2) + (3.0/2.0)*P21*std::pow(n21, 2)*n22/std::pow(eps2 + rc2, 2) + (3.0/2.0)*P22*n21*std::pow(n22, 2)/std::pow(eps2 + rc2, 2) + (3.0/2.0)*(P21*n22 + P22*n21)/std::pow(eps2 + rc2, 2);
        amrex::Real Atilde2UU22 = -3.0/2.0*P20*n20*(1 - std::pow(n22, 2))/std::pow(eps2 + rc2, 2) - 3.0/2.0*P21*n21*(1 - std::pow(n22, 2))/std::pow(eps2 + rc2, 2) - 3.0/2.0*P22*n22*(1 - std::pow(n22, 2))/std::pow(eps2 + rc2, 2) + 3*P22*n22/std::pow(eps2 + rc2, 2);
        amrex::Real AtildeUU00 = Atilde1UU00 + Atilde2UU00;
        amrex::Real AtildeUU01 = Atilde1UU01 + Atilde2UU01;
        amrex::Real AtildeUU02 = Atilde1UU02 + Atilde2UU02;
        amrex::Real AtildeUU10 = Atilde1UU10 + Atilde2UU10;
        amrex::Real AtildeUU11 = Atilde1UU11 + Atilde2UU11;
        amrex::Real AtildeUU12 = Atilde1UU12 + Atilde2UU12;
        amrex::Real AtildeUU20 = Atilde1UU20 + Atilde2UU20;
        amrex::Real AtildeUU21 = Atilde1UU21 + Atilde2UU21;
        amrex::Real AtildeUU22 = Atilde1UU22 + Atilde2UU22;
        
        diag(i, j, k, Diag::phi) = std::log(M1/rc1 + M2/rc2 + 1 + state_fab(i, j, k, Idx::u));

        diag(i, j, k, Diag::gamtildeLL00) = 1;

        diag(i, j, k, Diag::gamtildeLL01) = 0;

        diag(i, j, k, Diag::gamtildeLL02) = 0;

        diag(i, j, k, Diag::gamtildeLL11) = 1;

        diag(i, j, k, Diag::gamtildeLL12) = 0;

        diag(i, j, k, Diag::gamtildeLL22) = 1;

        diag(i, j, k, Diag::Khat) = 0;

        diag(i, j, k, Diag::AtildeLL00) = AtildeUU00;

        diag(i, j, k, Diag::AtildeLL01) = AtildeUU01;

        diag(i, j, k, Diag::AtildeLL02) = AtildeUU02;

        diag(i, j, k, Diag::AtildeLL11) = AtildeUU11;

        diag(i, j, k, Diag::AtildeLL12) = AtildeUU12;

        diag(i, j, k, Diag::AtildeLL22) = AtildeUU22;

        diag(i, j, k, Diag::theta) = 0;

        diag(i, j, k, Diag::GamtildeU0) = 0;

        diag(i, j, k, Diag::GamtildeU1) = 0;

        diag(i, j, k, Diag::GamtildeU2) = 0;

        diag(i, j, k, Diag::alpha) = std::pow(M1/rc1 + M2/rc2 + 1 + state_fab(i, j, k, Idx::u), -2);

        diag(i, j, k, Diag::betaU0) = 0;

        diag(i, j, k, Diag::betaU1) = 0;

        diag(i, j, k, Diag::betaU2) = 0;

        diag(i, j, k, Diag::BU0) = 0;

        diag(i, j, k, Diag::BU1) = 0;

        diag(i, j, k, Diag::BU2) = 0;

}
#endif
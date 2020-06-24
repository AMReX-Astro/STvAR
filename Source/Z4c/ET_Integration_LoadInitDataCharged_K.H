#ifndef ET_INTEGRATION_LOADINITDATA_K_H 
#define ET_INTEGRATION_LOADINITDATA_K_H 

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup.H> 

AMREX_GPU_DEVICE 
inline 
void 
initialize_from_data(int i, int j, int k, 
                     amrex::Array4<amrex::Real> const& state, 
                     amrex::Array4<amrex::Real const> const& initial_data, 
                     const amrex::GeometryData& geom) noexcept 
{

        const auto domain_xlo = geom.ProbLo(); 

        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(0) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(0) + domain_xlo[2];
        amrex::Real u = initial_data(i, j, k, InitIdx::u);
        amrex::Real psi = initial_data(i, j, k, InitIdx::psi);
        amrex::Real kappa = initial_data(i, j, k, InitIdx::kappa);
        amrex::Real phi = initial_data(i, j, k, InitIdx::phi);
        
        amrex::Real r = std::sqrt(std::pow(x0, 2) + std::pow(x1, 2) + std::pow(x2, 2));
        amrex::Real c10 = 1.168642873;
        amrex::Real c11 = 0.0;
        amrex::Real c12 = 0.0;
        amrex::Real c20 = -1.168642873;
        amrex::Real c21 = 0.0;
        amrex::Real c22 = 0.0;
        amrex::Real M1 = 0.45300000000000001;
        amrex::Real M2 = 0.45300000000000001;
        amrex::Real P10 = 0.0;
        amrex::Real P11 = 0.33319174979999999;
        amrex::Real P12 = 0.0;
        amrex::Real P20 = 0.0;
        amrex::Real P21 = -0.33319174979999999;
        amrex::Real P22 = 0.0;
        amrex::Real S10 = 0;
        amrex::Real S11 = 0;
        amrex::Real S12 = 0.5;
        amrex::Real S20 = 0;
        amrex::Real S21 = 0;
        amrex::Real S22 = 0.5;
        amrex::Real Q1 = 0.3;
        amrex::Real Q2 = -0.3;
        amrex::Real xc10 = -c10 + x0;
        amrex::Real xc11 = -c11 + x1;
        amrex::Real xc12 = -c12 + x2;
        amrex::Real xc20 = -c20 + x0;
        amrex::Real xc21 = -c21 + x1;
        amrex::Real xc22 = -c22 + x2;
        amrex::Real rc1 = std::sqrt(std::pow(xc10, 2) + std::pow(xc11, 2) + std::pow(xc12, 2));
        amrex::Real rc2 = std::sqrt(std::pow(xc20, 2) + std::pow(xc21, 2) + std::pow(xc22, 2));
        amrex::Real n10 = xc10/rc1;
        amrex::Real n11 = xc11/rc1;
        amrex::Real n12 = xc12/rc1;
        amrex::Real n20 = xc20/rc2;
        amrex::Real n21 = xc21/rc2;
        amrex::Real n22 = xc22/rc2;
        amrex::Real Ebar1U0 = Q1*n10/std::pow(rc1, 2);
        amrex::Real Ebar1U1 = Q1*n11/std::pow(rc1, 2);
        amrex::Real Ebar1U2 = Q1*n12/std::pow(rc1, 2);
        amrex::Real Ebar2U0 = Q2*n20/std::pow(rc2, 2);
        amrex::Real Ebar2U1 = Q2*n21/std::pow(rc2, 2);
        amrex::Real Ebar2U2 = Q2*n22/std::pow(rc2, 2);
        amrex::Real EbarU0 = Ebar1U0 + Ebar2U0;
        amrex::Real EbarU1 = Ebar1U1 + Ebar2U1;
        amrex::Real EbarU2 = Ebar1U2 + Ebar2U2;
        amrex::Real Enbar = 0 + (1.0/8.0)*std::pow(EbarU0, 2)/M_PI + (1.0/8.0)*std::pow(EbarU1, 2)/M_PI + (1.0/8.0)*std::pow(EbarU2, 2)/M_PI;
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
        amrex::Real alphafunc = 1.0/((1.0/2.0)*M1/rc1 + (1.0/2.0)*M2/rc2);
        //amrex::Real psi = u + 1 + 1.0/alphafunc;
        /*amrex::Real r = std::sqrt(std::pow(x0, 2) + std::pow(x1, 2) + std::pow(x2, 2));
        amrex::Real c10 = 1.168642873;
        amrex::Real c11 = 0.0;
        amrex::Real c12 = 0.0;
        amrex::Real c20 = -1.168642873;
        amrex::Real c21 = 0.0;
        amrex::Real c22 = 0.0;
        amrex::Real M1 = 0.45300000000000001;
        amrex::Real M2 = 0.45300000000000001;
        amrex::Real P10 = 0.0;
        amrex::Real P11 = 0.33319174979999999;
        amrex::Real P12 = 0.0;
        amrex::Real P20 = 0.0;
        amrex::Real P21 = -0.33319174979999999;
        amrex::Real P22 = 0.0;
        amrex::Real xc10 = -c10 + x0;
        amrex::Real xc11 = -c11 + x1;
        amrex::Real xc12 = -c12 + x2;
        amrex::Real xc20 = -c20 + x0;
        amrex::Real xc21 = -c21 + x1;
        amrex::Real xc22 = -c22 + x2;
        amrex::Real rc1 = std::sqrt(std::pow(xc10, 2) + std::pow(xc11, 2) + std::pow(xc12, 2));
        amrex::Real rc2 = std::sqrt(std::pow(xc20, 2) + std::pow(xc21, 2) + std::pow(xc22, 2));
        amrex::Real n10 = xc10/rc1;
        amrex::Real n11 = xc11/rc1;
        amrex::Real n12 = xc12/rc1;
        amrex::Real n20 = xc20/rc2;
        amrex::Real n21 = xc21/rc2;
        amrex::Real n22 = xc22/rc2;
        amrex::Real Abar1UU00 = -3.0/2.0*P10*n10*(1 - std::pow(n10, 2))/std::pow(rc1, 2) + 3*P10*n10/std::pow(rc1, 2) - 3.0/2.0*P11*n11*(1 - std::pow(n10, 2))/std::pow(rc1, 2) - 3.0/2.0*P12*n12*(1 - std::pow(n10, 2))/std::pow(rc1, 2);
        amrex::Real Abar1UU01 = (3.0/2.0)*P10*std::pow(n10, 2)*n11/std::pow(rc1, 2) + (3.0/2.0)*P11*n10*std::pow(n11, 2)/std::pow(rc1, 2) + (3.0/2.0)*P12*n10*n11*n12/std::pow(rc1, 2) + (3.0/2.0)*(P10*n11 + P11*n10)/std::pow(rc1, 2);
        amrex::Real Abar1UU02 = (3.0/2.0)*P10*std::pow(n10, 2)*n12/std::pow(rc1, 2) + (3.0/2.0)*P11*n10*n11*n12/std::pow(rc1, 2) + (3.0/2.0)*P12*n10*std::pow(n12, 2)/std::pow(rc1, 2) + (3.0/2.0)*(P10*n12 + P12*n10)/std::pow(rc1, 2);
        amrex::Real Abar1UU10 = (3.0/2.0)*P10*std::pow(n10, 2)*n11/std::pow(rc1, 2) + (3.0/2.0)*P11*n10*std::pow(n11, 2)/std::pow(rc1, 2) + (3.0/2.0)*P12*n10*n11*n12/std::pow(rc1, 2) + (3.0/2.0)*(P10*n11 + P11*n10)/std::pow(rc1, 2);
        amrex::Real Abar1UU11 = -3.0/2.0*P10*n10*(1 - std::pow(n11, 2))/std::pow(rc1, 2) - 3.0/2.0*P11*n11*(1 - std::pow(n11, 2))/std::pow(rc1, 2) + 3*P11*n11/std::pow(rc1, 2) - 3.0/2.0*P12*n12*(1 - std::pow(n11, 2))/std::pow(rc1, 2);
        amrex::Real Abar1UU12 = (3.0/2.0)*P10*n10*n11*n12/std::pow(rc1, 2) + (3.0/2.0)*P11*std::pow(n11, 2)*n12/std::pow(rc1, 2) + (3.0/2.0)*P12*n11*std::pow(n12, 2)/std::pow(rc1, 2) + (3.0/2.0)*(P11*n12 + P12*n11)/std::pow(rc1, 2);
        amrex::Real Abar1UU20 = (3.0/2.0)*P10*std::pow(n10, 2)*n12/std::pow(rc1, 2) + (3.0/2.0)*P11*n10*n11*n12/std::pow(rc1, 2) + (3.0/2.0)*P12*n10*std::pow(n12, 2)/std::pow(rc1, 2) + (3.0/2.0)*(P10*n12 + P12*n10)/std::pow(rc1, 2);
        amrex::Real Abar1UU21 = (3.0/2.0)*P10*n10*n11*n12/std::pow(rc1, 2) + (3.0/2.0)*P11*std::pow(n11, 2)*n12/std::pow(rc1, 2) + (3.0/2.0)*P12*n11*std::pow(n12, 2)/std::pow(rc1, 2) + (3.0/2.0)*(P11*n12 + P12*n11)/std::pow(rc1, 2);
        amrex::Real Abar1UU22 = -3.0/2.0*P10*n10*(1 - std::pow(n12, 2))/std::pow(rc1, 2) - 3.0/2.0*P11*n11*(1 - std::pow(n12, 2))/std::pow(rc1, 2) - 3.0/2.0*P12*n12*(1 - std::pow(n12, 2))/std::pow(rc1, 2) + 3*P12*n12/std::pow(rc1, 2);
        amrex::Real Abar2UU00 = -3.0/2.0*P20*n20*(1 - std::pow(n20, 2))/std::pow(rc2, 2) + 3*P20*n20/std::pow(rc2, 2) - 3.0/2.0*P21*n21*(1 - std::pow(n20, 2))/std::pow(rc2, 2) - 3.0/2.0*P22*n22*(1 - std::pow(n20, 2))/std::pow(rc2, 2);
        amrex::Real Abar2UU01 = (3.0/2.0)*P20*std::pow(n20, 2)*n21/std::pow(rc2, 2) + (3.0/2.0)*P21*n20*std::pow(n21, 2)/std::pow(rc2, 2) + (3.0/2.0)*P22*n20*n21*n22/std::pow(rc2, 2) + (3.0/2.0)*(P20*n21 + P21*n20)/std::pow(rc2, 2);
        amrex::Real Abar2UU02 = (3.0/2.0)*P20*std::pow(n20, 2)*n22/std::pow(rc2, 2) + (3.0/2.0)*P21*n20*n21*n22/std::pow(rc2, 2) + (3.0/2.0)*P22*n20*std::pow(n22, 2)/std::pow(rc2, 2) + (3.0/2.0)*(P20*n22 + P22*n20)/std::pow(rc2, 2);
        amrex::Real Abar2UU10 = (3.0/2.0)*P20*std::pow(n20, 2)*n21/std::pow(rc2, 2) + (3.0/2.0)*P21*n20*std::pow(n21, 2)/std::pow(rc2, 2) + (3.0/2.0)*P22*n20*n21*n22/std::pow(rc2, 2) + (3.0/2.0)*(P20*n21 + P21*n20)/std::pow(rc2, 2);
        amrex::Real Abar2UU11 = -3.0/2.0*P20*n20*(1 - std::pow(n21, 2))/std::pow(rc2, 2) - 3.0/2.0*P21*n21*(1 - std::pow(n21, 2))/std::pow(rc2, 2) + 3*P21*n21/std::pow(rc2, 2) - 3.0/2.0*P22*n22*(1 - std::pow(n21, 2))/std::pow(rc2, 2);
        amrex::Real Abar2UU12 = (3.0/2.0)*P20*n20*n21*n22/std::pow(rc2, 2) + (3.0/2.0)*P21*std::pow(n21, 2)*n22/std::pow(rc2, 2) + (3.0/2.0)*P22*n21*std::pow(n22, 2)/std::pow(rc2, 2) + (3.0/2.0)*(P21*n22 + P22*n21)/std::pow(rc2, 2);
        amrex::Real Abar2UU20 = (3.0/2.0)*P20*std::pow(n20, 2)*n22/std::pow(rc2, 2) + (3.0/2.0)*P21*n20*n21*n22/std::pow(rc2, 2) + (3.0/2.0)*P22*n20*std::pow(n22, 2)/std::pow(rc2, 2) + (3.0/2.0)*(P20*n22 + P22*n20)/std::pow(rc2, 2);
        amrex::Real Abar2UU21 = (3.0/2.0)*P20*n20*n21*n22/std::pow(rc2, 2) + (3.0/2.0)*P21*std::pow(n21, 2)*n22/std::pow(rc2, 2) + (3.0/2.0)*P22*n21*std::pow(n22, 2)/std::pow(rc2, 2) + (3.0/2.0)*(P21*n22 + P22*n21)/std::pow(rc2, 2);
        amrex::Real Abar2UU22 = -3.0/2.0*P20*n20*(1 - std::pow(n22, 2))/std::pow(rc2, 2) - 3.0/2.0*P21*n21*(1 - std::pow(n22, 2))/std::pow(rc2, 2) - 3.0/2.0*P22*n22*(1 - std::pow(n22, 2))/std::pow(rc2, 2) + 3*P22*n22/std::pow(rc2, 2);
        amrex::Real AbarUU00 = Abar1UU00 + Abar2UU00;
        amrex::Real AbarUU01 = Abar1UU01 + Abar2UU01;
        amrex::Real AbarUU02 = Abar1UU02 + Abar2UU02;
        amrex::Real AbarUU10 = Abar1UU10 + Abar2UU10;
        amrex::Real AbarUU11 = Abar1UU11 + Abar2UU11;
        amrex::Real AbarUU12 = Abar1UU12 + Abar2UU12;
        amrex::Real AbarUU20 = Abar1UU20 + Abar2UU20;
        amrex::Real AbarUU21 = Abar1UU21 + Abar2UU21;
        amrex::Real AbarUU22 = Abar1UU22 + Abar2UU22;
        amrex::Real alphafunc = 1.0/((1.0/2.0)*M1/rc1 + (1.0/2.0)*M2/rc2);
        amrex::Real psi = u + 1 + 1.0/alphafunc;*/
        state(i, j, k, Idx::chi) = std::pow(psi, -4);

        state(i, j, k, Idx::gamtildeLL00) = 1;

        state(i, j, k, Idx::gamtildeLL01) = 0;

        state(i, j, k, Idx::gamtildeLL02) = 0;

        state(i, j, k, Idx::gamtildeLL11) = 1;

        state(i, j, k, Idx::gamtildeLL12) = 0;

        state(i, j, k, Idx::gamtildeLL22) = 1;

        state(i, j, k, Idx::Khat) = 0;

        state(i, j, k, Idx::AtildeLL00) = AbarUU00/std::pow(psi, 6);

        state(i, j, k, Idx::AtildeLL01) = AbarUU01/std::pow(psi, 6);

        state(i, j, k, Idx::AtildeLL02) = AbarUU02/std::pow(psi, 6);

        state(i, j, k, Idx::AtildeLL11) = AbarUU11/std::pow(psi, 6);

        state(i, j, k, Idx::AtildeLL12) = AbarUU12/std::pow(psi, 6);

        state(i, j, k, Idx::AtildeLL22) = AbarUU22/std::pow(psi, 6);

        state(i, j, k, Idx::theta) = 0;

        state(i, j, k, Idx::GamtildeU0) = 0;

        state(i, j, k, Idx::GamtildeU1) = 0;

        state(i, j, k, Idx::GamtildeU2) = 0;
        
        state(i, j, k, Idx::EU0) = EbarU0/std::pow(psi,6);

        state(i, j, k, Idx::EU1) = EbarU1/std::pow(psi,6);

        state(i, j, k, Idx::EU2) = EbarU2/std::pow(psi,6);
        
        state(i, j, k, Idx::BU0) = 0;

        state(i, j, k, Idx::BU1) = 0;

        state(i, j, k, Idx::BU2) = 0;

        state(i, j, k, Idx::alpha) = std::pow(psi, -2);

        state(i, j, k, Idx::betaU0) = 0;

        state(i, j, k, Idx::betaU1) = 0;

        state(i, j, k, Idx::betaU2) = 0;

}
#endif
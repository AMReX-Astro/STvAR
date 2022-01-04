#ifndef ET_INTEGRATION_INIT_K_H
#define ET_INTEGRATION_INIT_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_init(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab, 
        amrex::Real time, const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo(); 
        
        amrex::Real x_U_0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x_U_1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x_U_2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        /*Robust Stability*/
        /*
        amrex::Real epsilon = std::pow(10,-8)/16; //std::pow(10,-1); 
        
        state_fab(i, j, k, Idx::chi) = 1 + 2*epsilon*(amrex::Random()-0.5);

        state_fab(i, j, k, Idx::gamtilde_LL_00) = 1 + 2*epsilon*(amrex::Random()-0.5);
        state_fab(i, j, k, Idx::gamtilde_LL_01) = 0 + 2*epsilon*(amrex::Random()-0.5);
        state_fab(i, j, k, Idx::gamtilde_LL_02) = 0 + 2*epsilon*(amrex::Random()-0.5);
        state_fab(i, j, k, Idx::gamtilde_LL_11) = 1 + 2*epsilon*(amrex::Random()-0.5);
        state_fab(i, j, k, Idx::gamtilde_LL_12) = 0 + 2*epsilon*(amrex::Random()-0.5);
        state_fab(i, j, k, Idx::gamtilde_LL_22) = 1 + 2*epsilon*(amrex::Random()-0.5);

        state_fab(i, j, k, Idx::Khat) = 0 + 2*epsilon*(amrex::Random()-0.5);

        state_fab(i, j, k, Idx::Atilde_LL_00) = 0 + 2*epsilon*(amrex::Random()-0.5);
        state_fab(i, j, k, Idx::Atilde_LL_01) = 0 + 2*epsilon*(amrex::Random()-0.5);
        state_fab(i, j, k, Idx::Atilde_LL_02) = 0 + 2*epsilon*(amrex::Random()-0.5);
        state_fab(i, j, k, Idx::Atilde_LL_11) = 0 + 2*epsilon*(amrex::Random()-0.5);
        state_fab(i, j, k, Idx::Atilde_LL_12) = 0 + 2*epsilon*(amrex::Random()-0.5);
        state_fab(i, j, k, Idx::Atilde_LL_22) = 0 + 2*epsilon*(amrex::Random()-0.5);

        state_fab(i, j, k, Idx::theta) = 0 + 2*epsilon*(amrex::Random()-0.5);

        state_fab(i, j, k, Idx::Gamtilde_U_0) = 0 + 2*epsilon*(amrex::Random()-0.5);
        state_fab(i, j, k, Idx::Gamtilde_U_1) = 0 + 2*epsilon*(amrex::Random()-0.5);
        state_fab(i, j, k, Idx::Gamtilde_U_2) = 0 + 2*epsilon*(amrex::Random()-0.5);

        state_fab(i, j, k, Idx::alpha) = 1 + 2*epsilon*(amrex::Random()-0.5);

        state_fab(i, j, k, Idx::beta_U_0) = 0 + 2*epsilon*(amrex::Random()-0.5);
        state_fab(i, j, k, Idx::beta_U_1) = 0 + 2*epsilon*(amrex::Random()-0.5);
        state_fab(i, j, k, Idx::beta_U_2) = 0 + 2*epsilon*(amrex::Random()-0.5);
        
        */
        /*Linear Wave Test*/
        /*
        amrex::Real A = std::pow(10,-8);
        amrex::Real d = 1.28;
        amrex::Real b = A*std::sin(2*M_PI*x_U_0/d);
        amrex::Real dtb = -(2*M_PI*A/d)*std::cos(2*M_PI*x_U_0/d);
        
        amrex::Real chi = std::pow((1-b*b),-1.0/3.0);
        
        amrex::Real Ksclr = -b*dtb/(1-b*b);
        
        amrex::Real gamma_LL_00 = 1;
        amrex::Real gamma_LL_11 = 1+b;
        amrex::Real gamma_LL_22 = 1-b;
        
        amrex::Real K_LL_00 = 0;
        amrex::Real K_LL_11 = 0.5*dtb;
        amrex::Real K_LL_22 = -0.5*dtb;
         
        
        state_fab(i, j, k, Idx::chi) = chi;
        state_fab(i, j, k, Idx::phi) = -(1.0/4.0)*std::log(chi);

        state_fab(i, j, k, Idx::gamtilde_LL_00) = chi;
        state_fab(i, j, k, Idx::gamtilde_LL_01) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_02) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_11) = chi*(1+b);
        state_fab(i, j, k, Idx::gamtilde_LL_12) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_22) = chi*(1-b);

        state_fab(i, j, k, Idx::Khat) = Ksclr;

        state_fab(i, j, k, Idx::Atilde_LL_00) = chi*(K_LL_00 - (1.0/3.0)*gamma_LL_00*Ksclr);
        state_fab(i, j, k, Idx::Atilde_LL_01) = 0;
        state_fab(i, j, k, Idx::Atilde_LL_02) = 0;
        state_fab(i, j, k, Idx::Atilde_LL_11) = chi*(K_LL_11 - (1.0/3.0)*gamma_LL_11*Ksclr);
        state_fab(i, j, k, Idx::Atilde_LL_12) = 0;
        state_fab(i, j, k, Idx::Atilde_LL_22) = chi*(K_LL_22 - (1.0/3.0)*gamma_LL_22*Ksclr);

        state_fab(i, j, k, Idx::theta) = 0;

        state_fab(i, j, k, Idx::Gamtilde_U_0) = 0;
        state_fab(i, j, k, Idx::Gamtilde_U_1) = 0;
        state_fab(i, j, k, Idx::Gamtilde_U_2) = 0;

        state_fab(i, j, k, Idx::alpha) = 1;

        state_fab(i, j, k, Idx::beta_U_0) = 0;
        state_fab(i, j, k, Idx::beta_U_1) = 0;
        state_fab(i, j, k, Idx::beta_U_2) = 0;
        
        */
        /*Unshifted Gauge wave Test*/
        
        amrex::Real A = 0.01;
        amrex::Real d = 1.28;
        amrex::Real b = A*std::sin(2*M_PI*x_U_0/d);
        amrex::Real dtb = -(2*M_PI*A/d)*std::cos(2*M_PI*x_U_0/d);
        amrex::Real dxb = (2*M_PI*A/d)*std::cos(2*M_PI*x_U_0/d);
        
        amrex::Real chi = 1/std::pow(1-b,1.0/3.0);
        
        amrex::Real Ksclr = dtb/(2*std::pow(1-b,3.0/2.0));
        
        amrex::Real gamma_LL_00 = 1-b;
        amrex::Real gamma_LL_11 = 1;
        amrex::Real gamma_LL_22 = 1;
        
        amrex::Real K_LL_00 = dtb/(2*std::pow(1-b,1.0/2.0));
        amrex::Real K_LL_11 = 0;
        amrex::Real K_LL_22 = 0;
         
        
        state_fab(i, j, k, Idx::chi) = chi;

        state_fab(i, j, k, Idx::gamtilde_LL_00) = chi*(1-b);
        state_fab(i, j, k, Idx::gamtilde_LL_01) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_02) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_11) = chi;
        state_fab(i, j, k, Idx::gamtilde_LL_12) = 0;
        state_fab(i, j, k, Idx::gamtilde_LL_22) = chi;

        state_fab(i, j, k, Idx::Khat) = Ksclr;

        state_fab(i, j, k, Idx::Atilde_LL_00) = chi*(K_LL_00 - (1.0/3.0)*gamma_LL_00*Ksclr);
        state_fab(i, j, k, Idx::Atilde_LL_01) = 0;
        state_fab(i, j, k, Idx::Atilde_LL_02) = 0;
        state_fab(i, j, k, Idx::Atilde_LL_11) = chi*(K_LL_11 - (1.0/3.0)*gamma_LL_11*Ksclr);
        state_fab(i, j, k, Idx::Atilde_LL_12) = 0;
        state_fab(i, j, k, Idx::Atilde_LL_22) = chi*(K_LL_22 - (1.0/3.0)*gamma_LL_22*Ksclr);

        state_fab(i, j, k, Idx::theta) = 0;

        state_fab(i, j, k, Idx::Gamtilde_U_0) = -2*dxb/(3*std::pow(1-b,5.0/3.0));
        state_fab(i, j, k, Idx::Gamtilde_U_1) = 0;
        state_fab(i, j, k, Idx::Gamtilde_U_2) = 0;

        state_fab(i, j, k, Idx::alpha) = std::pow(1-b,1.0/2.0);

        state_fab(i, j, k, Idx::beta_U_0) = 0;
        state_fab(i, j, k, Idx::beta_U_1) = 0;
        state_fab(i, j, k, Idx::beta_U_2) = 0;
        
}
#endif
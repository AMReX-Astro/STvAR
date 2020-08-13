#ifndef ET_INTEGRATION_AMR_CELL_TAG_K_H
#define ET_INTEGRATION_AMR_CELL_TAG_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
bool
state_is_tagged(int i, int j, int k, 
        amrex::Array4<amrex::Real const> const& state_fab,
        amrex::Real error_threshold,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{

        const auto domain_xlo = geom.ProbLo(); 

        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        amrex::Real u = state_fab(i, j, k, Idx::u);

        amrex::Real r = std::sqrt(std::pow(x0, 2) + std::pow(x1, 2) + std::pow(x2, 2));

        amrex::Real c10 = 1.168642873;
        amrex::Real c11 = 0.0;
        amrex::Real c12 = 0.0;
        amrex::Real c20 = -1.168642873;
        amrex::Real c21 = 0.0;
        amrex::Real c22 = 0.0;
        amrex::Real M1 = 0.755;

        amrex::Real M2 = 0.151;

        amrex::Real P10 = 0.0;
        amrex::Real P11 = 0.23319174979999999;
        amrex::Real P12 = 0.0;
        amrex::Real P20 = 0.0;
        amrex::Real P21 = -0.23319174979999999;
        amrex::Real P22 = 0.0;
        amrex::Real S10 = 0;
        amrex::Real S11 = 0;
        amrex::Real S12 = 0.10000000000000001;
        amrex::Real S20 = 0;
        amrex::Real S21 = 0.050000000000000003;
        amrex::Real S22 = 0.02;
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
        amrex::Real alpha = 1.0/((1.0/2.0)*M1/rc1 + (1.0/2.0)*M2/rc2);

        amrex::Real psi = u + 1 + 1.0/alpha;

        amrex::Real alphaTAG = std::pow(psi, -2);

        return alphaTAG < error_threshold;
}
#endif
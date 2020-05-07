#ifndef ET_INTEGRATION_INIT_K_H 
#define ET_INTEGRATION_INIT_K_H 

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_init(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab, 
        amrex::Real time, const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo(); 

        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0];
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1];
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2];
        
        amrex::Real xp = x0-3.0;
        amrex::Real xm = x0+3.0;
        amrex::Real r = std::sqrt(x1*x1+x2*x2);

        state_fab(i, j, k, Idx::u) = (1.0/(std::pow(5*xp,2)+50)+1.0/(std::pow(5*xm,2)+50)+0.002)*(1.0/(std::pow(5*r,2)+50)+0.001)/(1.0/(50)+0.001);

}
#endif
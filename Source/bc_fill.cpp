#include <bc_fill.H>
#include "ET_Integration.H"

using namespace amrex;

void AmrCoreFillCpu (Box const& bx, Array4<Real> const& data,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp)
{
    // do something for external Dirichlet (BCType::ext_dir)
    const auto& domain = geom.Domain();
    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    ParallelFor(bx, numcomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int mcomp) {
        
        const auto domain_xlo = geom.ProbLo();
        const auto domain_xhi = geom.ProbHi();
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        amrex::Real x0blo = (0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1blo = (0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2blo = (0.5)*geom.CellSize(2) + domain_xlo[2];
        
        amrex::Real x0blop1 = (1.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1blop1 = (1.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2blop1 = (1.5)*geom.CellSize(2) + domain_xlo[2];
        
        amrex::Real x0bhi = (-0.5)*geom.CellSize(0) + domain_xhi[0]; 
        amrex::Real x1bhi = (-0.5)*geom.CellSize(1) + domain_xhi[1]; 
        amrex::Real x2bhi = (-0.5)*geom.CellSize(2) + domain_xhi[2];
        
        amrex::Real x0bhim1 = (-1.5)*geom.CellSize(0) + domain_xhi[0]; 
        amrex::Real x1bhim1 = (-1.5)*geom.CellSize(1) + domain_xhi[1]; 
        amrex::Real x2bhim1 = (-1.5)*geom.CellSize(2) + domain_xhi[2];
        
        int n = dcomp + mcomp;
        if ((i < dom_lo.x && bcr[n].lo(0) == BCType::ext_dir) ||
            (i > dom_hi.x && bcr[n].hi(0) == BCType::ext_dir) ||

            (j < dom_lo.y && bcr[n].lo(1) == BCType::ext_dir) ||
            (j > dom_hi.y && bcr[n].hi(1) == BCType::ext_dir) ||

            (k < dom_lo.z && bcr[n].lo(2) == BCType::ext_dir) ||
            (k > dom_hi.z && bcr[n].hi(2) == BCType::ext_dir))
        {
            
            
        }
        
        
    });
}
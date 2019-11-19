#include <AMReX_BCUtil.H>
#include <AMReX_PhysBCFunct.H>

namespace amrex
{

namespace {

void dummy_cpu_fill_extdir (Box const& bx, Array4<Real> const& dest,
                            const int dcomp, const int numcomp,
                            GeometryData const& geom, const Real time,
                            const BCRec* bcr, const int bcomp,
                            const int orig_comp)
{
    // do something for external Dirichlet (BCType::ext_dir) if needed
    const auto& domain = geom.Domain();
    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    ParallelFor(bx, numcomp, [=](int i, int j, int k, int mcomp) {
        int n = dcomp + mcomp;

        if ((i < dom_lo.x && bcr[n].lo(0) == BCType::ext_dir) ||
            (i > dom_hi.x && bcr[n].hi(0) == BCType::ext_dir) ||

            (j < dom_lo.y && bcr[n].lo(1) == BCType::ext_dir) ||
            (j > dom_hi.y && bcr[n].hi(1) == BCType::ext_dir) ||

            (k < dom_lo.z && bcr[n].lo(2) == BCType::ext_dir) ||
            (k > dom_hi.z && bcr[n].hi(2) == BCType::ext_dir))
        {
                dest(i, j, k, n) = 0.0;
        }
    });
}

struct dummy_gpu_fill_extdir
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
        {
            // do something for external Dirichlet (BCType::ext_dir) if there are
            // NOT IMPLEMENTED FOR NOW
        }
};

}

void FillDomainBoundary (MultiFab& phi, const Geometry& geom, const Vector<BCRec>& bc)
{
    if (geom.isAllPeriodic()) return;
    if (phi.nGrow() == 0) return;

    AMREX_ALWAYS_ASSERT(phi.ixType().cellCentered());

    if (Gpu::inLaunchRegion())
    {
        GpuBndryFuncFab<dummy_gpu_fill_extdir> gpu_bndry_func(dummy_gpu_fill_extdir{});
        PhysBCFunct<GpuBndryFuncFab<dummy_gpu_fill_extdir> > physbcf
            (geom, bc, gpu_bndry_func);
        physbcf.FillBoundary(phi, 0, phi.nComp(), 0.0, 0);
    }
    else
    {
        CpuBndryFuncFab cpu_bndry_func(dummy_cpu_fill_extdir);;
        PhysBCFunct<CpuBndryFuncFab> physbcf(geom, bc, cpu_bndry_func);
        physbcf.FillBoundary(phi, 0, phi.nComp(), 0.0, 0);
    }
}

}

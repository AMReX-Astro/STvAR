#include <bc_fill.H>

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
        int n = dcomp + mcomp;

        if ((i < dom_lo.x && bcr[n].lo(0) == BCType::ext_dir) ||
            (i > dom_hi.x && bcr[n].hi(0) == BCType::ext_dir) ||

            (j < dom_lo.y && bcr[n].lo(1) == BCType::ext_dir) ||
            (j > dom_hi.y && bcr[n].hi(1) == BCType::ext_dir) ||

            (k < dom_lo.z && bcr[n].lo(2) == BCType::ext_dir) ||
            (k > dom_hi.z && bcr[n].hi(2) == BCType::ext_dir))
        {
                data(i, j, k, n) = 0.0;
        }
    });
}
#ifndef BCFILL_H
#define BCFILL_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

using namespace amrex;

void AmrCoreFillCpu (Box const& bx, Array4<Real> const& data,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp);

struct AmrCoreFillGpu
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& data,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
        {
            // do something for external Dirichlet (BCType::ext_dir)
            const auto& domain = geom.Domain();
            const auto& dom_lo = amrex::lbound(domain);
            const auto& dom_hi = amrex::ubound(domain);

            const auto cell = iv.dim3();
            const int i = cell.x;
            const int j = cell.y;
            const int k = cell.z;

            for (int n = dcomp; n < dcomp + numcomp; ++n) {
                if ((i < dom_lo.x && bcr[n].lo(0) == BCType::ext_dir) ||
                    (i > dom_hi.x && bcr[n].hi(0) == BCType::ext_dir) ||

                    (j < dom_lo.y && bcr[n].lo(1) == BCType::ext_dir) ||
                    (j > dom_hi.y && bcr[n].hi(1) == BCType::ext_dir) ||

                    (k < dom_lo.z && bcr[n].lo(2) == BCType::ext_dir) ||
                    (k > dom_hi.z && bcr[n].hi(2) == BCType::ext_dir))
                {
                        data(i, j, k, n) = 0.0;
                }
            }
        }
};

#endif

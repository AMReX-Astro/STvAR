#include "ET_Integration.H"
#include "ET_Integration_K.H"

using namespace amrex;

void fill_phi_rhs (MultiFab& phi_rhs_mf, MultiFab& phi_old_mf, const amrex::Geometry& geom)
{
  const auto dx = geom.CellSizeArray();

  for ( MFIter mfi(phi_rhs_mf); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.validbox();
    const auto ncomp = phi_old_mf.nComp();

    const auto& phi_rhs_fab = phi_rhs_mf.array(mfi);
    const auto& phi_old_fab = phi_old_mf.array(mfi);

    // For each grid, loop over all the valid points
    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      phi_rhs(i, j, k, n, phi_rhs_fab, phi_old_fab, dx);
    });
  }
}
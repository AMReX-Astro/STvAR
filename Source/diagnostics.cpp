#include "ET_Integration.H"

using namespace amrex;

void fill_state_diagnostics (MultiFab& diag_mf, const MultiFab& state_mf, const amrex::Geometry& geom)
{
  const auto dx = geom.CellSizeArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( MFIter mfi(diag_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& diag_fab = diag_mf.array(mfi);
    const auto& state_fab = state_mf.array(mfi);

    // For each grid, loop over all the valid points
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      state_diagnostics(i, j, k, diag_fab, state_fab, dx, geom.data());
    });
  }
}

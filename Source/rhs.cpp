#include "ET_Integration.H"

using namespace amrex;

void fill_state_rhs (MultiFab& rhs_mf, const MultiFab& state_mf, const amrex::Geometry& geom)
{
  const auto dx = geom.CellSizeArray();

  for ( MFIter mfi(rhs_mf); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.validbox();
    const auto ncomp = state_mf.nComp();

    const auto& rhs_fab = rhs_mf.array(mfi);
    const auto& state_fab = state_mf.array(mfi);

    // For each grid, loop over all the valid points
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      state_rhs(i, j, k, rhs_fab, state_fab, dx);
    });
  }
}

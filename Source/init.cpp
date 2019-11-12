#include "ET_Integration.H"

using namespace amrex;

void init(MultiFab& state_mf, Real time, const Geometry& geom)
{
    const auto dx = geom.CellSizeArray();
    int ncomp = state_mf.nComp();

    // Loop over grids 
    for ( MFIter mfi(state_mf); mfi.isValid(); ++mfi )
    {

      const Box& bx = mfi.validbox();

      const auto& state_fab = state_mf.array(mfi);

      // For each grid, loop over all the valid points
      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        state_init(i, j, k, state_fab, time, dx);
      });
    }

    // Fill ghost cells for each grid from valid regions of another grid
    state_mf.FillBoundary(geom.periodicity());

    Print() << "Domain initialized\n";
}

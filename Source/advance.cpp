#include "ET_Integration.H"

using namespace amrex;

void advance(MultiFab& state_new_mf, MultiFab& state_old_mf, Real time, Real dt, const Geometry& geom)
{
    int ncomp = state_new_mf.nComp();

    // Fill ghost cells for each grid from valid regions of another grid
    state_old_mf.FillBoundary(geom.periodicity());

    // Create a MultiFab containing the time integration RHS
    MultiFab rhs_mf(state_new_mf.boxArray(), state_new_mf.DistributionMap(), ncomp, 0);
    fill_state_rhs(rhs_mf, state_old_mf, geom);

    // Loop over grids to do a forward euler integration in time
    for ( MFIter mfi(state_new_mf); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.validbox();

      const auto& state_new_fab = state_new_mf.array(mfi);
      const auto& state_old_fab = state_old_mf.array(mfi);
      const auto&     rhs_fab =     rhs_mf.array(mfi);

      // For each grid, loop over all the valid points
      AMREX_FOR_4D(bx, ncomp, i, j, k, n,
      {
         // Right now rhs_fab(i,j,k,n) = 1 so this adds dt to every value in every time step
         state_new_fab(i,j,k,n) = state_old_fab(i,j,k,n) + dt * rhs_fab(i,j,k,n);
      });
    }

    // Fill ghost cells for each grid from valid regions of another grid
    state_new_mf.FillBoundary(geom.periodicity());
}

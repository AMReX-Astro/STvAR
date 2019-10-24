
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

using namespace amrex;

void advance_phi (MultiFab& phi_new_mf, MultiFab& phi_old_mf, Real time, Real dt, const Real* dx)
{
    int ncomp = phi_new_mf.nComp();

    // Fill ghost cells for each grid from valid regions of another grid
    phi_old_mf.FillBoundary();

    // HACK HACK HACK -- you will replace this
    MultiFab rhs_mf(phi_new_mf.boxArray(), phi_new_mf.DistributionMap(), ncomp, 0);
    rhs_mf.setVal(1.0);

    // Loop over grids 
    for ( MFIter mfi(phi_new_mf); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.validbox();

      const auto& phi_new_fab = phi_new_mf.array(mfi);
      const auto& phi_old_fab = phi_old_mf.array(mfi);
      const auto&     rhs_fab =     rhs_mf.array(mfi);

      // For each grid, loop over all the valid points
      AMREX_FOR_4D(bx, ncomp, i, j, k, n,
      {
         // Right now rhs_fab(i,j,k,n) = 1 so this adds to every value in every time step
         phi_new_fab(i,j,k,n) = phi_old_fab(i,j,k,n) + rhs_fab(i,j,k,n);
      });
    }

    // Fill ghost cells for each grid from valid regions of another grid
    phi_new_mf.FillBoundary();
}

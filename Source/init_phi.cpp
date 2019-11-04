#include "ET_Integration.H"
#include "ET_Integration_Indexes.H"

using namespace amrex;

void init_phi (MultiFab& state_mf, Real time, const Geometry& geom)
{
    const auto dx = geom.CellSizeArray();
    int ncomp = state_mf.nComp();

    // Loop over grids 
    for ( MFIter mfi(state_mf); mfi.isValid(); ++mfi )
    {

      const Box& bx = mfi.validbox();

      const auto& state_fab = state_mf.array(mfi);

      // For each grid, loop over all the valid points
      AMREX_FOR_3D(bx, i, j, k,
      {
         // init phi
         Real x = (i+0.5) * dx[0];
         state_fab(i, j, k, IPHI) = sin(2.0*M_PI*(x - time));
         // init pi
         /* state_fab(i, j, k, IPI) = -2.0*M_PI*cos(2.0*M_PI*(x-time)); */
         state_fab(i, j, k, IPI) = +2.0*M_PI*cos(2.0*M_PI*(x-time));
      });
    }

    // Fill ghost cells for each grid from valid regions of another grid
    state_mf.FillBoundary();

    std::cout << "Phi and Pi have been initialized" << std::endl;
}

#include "ET_Integration.H"

using namespace amrex;

void init_phi (MultiFab& phi_mf, Real time, const Real* dx)
{
    int ncomp = phi_mf.nComp();

    // Loop over grids 
    for ( MFIter mfi(phi_mf); mfi.isValid(); ++mfi )
    {

      const Box& bx = mfi.validbox();

      const auto& phi_fab = phi_mf.array(mfi);

      // For each grid, loop over all the valid points
      AMREX_FOR_4D(bx, ncomp, i, j, k, n,
      {
         Real x = (i+0.5) * dx[0];
         phi_fab(i,j,k,n) = sin(2.0*M_PI*(x - time));
      });
    }

    // Fill ghost cells for each grid from valid regions of another grid
    phi_mf.FillBoundary();

    std::cout << "Phi has been initialized" << std::endl;
}

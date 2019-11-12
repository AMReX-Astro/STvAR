#include "ET_Integration.H"
#include "AMReX_TimeIntegrator.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();
    return 0;
}

void main_main ()
{
    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = amrex::second();

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size, plot_int, nsteps;
    Real cfl = 0.9;
    Real end_time = 1.0;
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // periodic in all direction by default

    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of
        //   a square (or cubic) domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size",max_grid_size);

        // Default plot_int to -1, allow us to set it to something else in the inputs file
        //  If plot_int < 0 then no plot files will be writtenq
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // Query domain periodicity
        pp.queryarr("is_periodic", is_periodic);

        // Read CFL number
        pp.query("cfl", cfl);

        // Default nsteps to 10, allow us to set it to something else in the inputs file
        nsteps = 10;
        pp.query("nsteps", nsteps);

        // Time at end of simulation
        pp.query("end_time", end_time);
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL(0.0,0.0,0.0)},
                         {AMREX_D_DECL(1.0,1.0,1.0)});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    const Real* dx = geom.CellSize();

    // Nghost = number of ghost cells for each array
    int Nghost = NUM_GHOST_CELLS;

    // Ncomp = number of components for each array
    int Ncomp  = Idx::NumScalars;

    // Initialize variable names
    Variable::Initialize();

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate two state multifabs; one will store the old state, the other the new.
    MultiFab state_old(ba, dm, Ncomp, Nghost);
    MultiFab state_new(ba, dm, Ncomp, Nghost);

    // time = starting time in the simulation
    Real time = 0.0;

    // Initialize state_new by calling a C++ initializing routine.
    init(state_new, time, geom);

    // Compute the time step
    Real dt = cfl*dx[0]*dx[0] / (2.0*AMREX_SPACEDIM);

    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
    if (plot_int > 0)
    {
        int n = 0;
        const std::string& pltfile = amrex::Concatenate("plt",n,7);
        WriteSingleLevelPlotfile(pltfile, state_new, Variable::names, geom, time, 0);
    }

    // Create integrator with the old state, new state, and new state time
    TimeIntegrator integrator(state_old, state_new, time);

    // Create a RHS source function we will integrate
    auto source_fun = [&](MultiFab& rhs, const MultiFab& state, const Real time){
        fill_state_rhs(rhs, state, geom);
    };

    // Create a function to call after updating a state
    auto post_update_fun = [&](MultiFab& S_data){
        // Fill ghost cells for S_data from interior & periodic BCs
        S_data.FillBoundary(geom.periodicity());
    };

    // Create a post-timestep function
    auto post_timestep_fun = [&](){
        // Tell the I/O Processor to write out which step we're doing
        const int n = integrator.get_step_number();
        amrex::Print() << "Advanced step " << n << "\n";

        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if (plot_int > 0 && n % plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",n,7);
            WriteSingleLevelPlotfile(pltfile, integrator.get_new_data(), Variable::names, geom, integrator.get_time(), n);
        }
    };

    integrator.set_rhs(source_fun);
    integrator.set_post_update(post_update_fun);
    integrator.set_post_timestep(post_timestep_fun);
    integrator.integrate(dt, end_time, nsteps);

    // Write a final plotfile
    {
        const std::string& pltfile = "plt_End_Simulation";
        WriteSingleLevelPlotfile(pltfile, integrator.get_new_data(), Variable::names, geom, integrator.get_time(), integrator.get_step_number());
    }

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real stop_time = amrex::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
}

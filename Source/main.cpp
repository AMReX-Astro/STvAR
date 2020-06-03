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
    int n_cell, max_grid_size, plot_int, diag_int, nsteps, elliptic, checkpoint_int, initialize_from_data, write_checkpoint;
    Real cfl = 0.9;
    Real end_time = 1.0;
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // periodic in all direction by default
    Vector<int> domain_lo_bc_types(AMREX_SPACEDIM, BCType::int_dir);
    Vector<int> domain_hi_bc_types(AMREX_SPACEDIM, BCType::int_dir);

    Vector<Real> domain_lo(AMREX_SPACEDIM,0.0);
    Vector<Real> domain_hi(AMREX_SPACEDIM,1.0);

    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of
        //   a square (or cubic) domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size",max_grid_size);

        // Default plot_int to -1, allow us to set it to something else in the inputs file
        //  If plot_int < 0 then no plot files will be written
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // Default diag_int to -1, allow us to set it to something else in the inputs file
        //  If diag_int < 0 then no diagnostic files will be written
        diag_int = -1;
        pp.query("diag_int",diag_int);

        // Query domain periodicity
        pp.queryarr("is_periodic", is_periodic);
        pp.queryarr("domain_lo_bc_types", domain_lo_bc_types);
        pp.queryarr("domain_hi_bc_types", domain_hi_bc_types);

        // Query domain physical lo, hi bounds
        pp.queryarr("domain_lo", domain_lo);
        pp.queryarr("domain_hi", domain_hi);

        // Read CFL number
        pp.query("cfl", cfl);

        // Default nsteps to 10, allow us to set it to something else in the inputs file
        nsteps = 10;
        pp.query("nsteps", nsteps);

        // Time at end of simulation
        pp.query("end_time", end_time);
        
        pp.query("elliptic", elliptic);
        
        pp.query("initialize_from_data", initialize_from_data);
        pp.query("write_checkpoint",write_checkpoint);
        
        checkpoint_int = -1;
        pp.query("checkpoint_int",checkpoint_int);
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_index_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_index_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_index_lo, dom_index_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

       // This defines the physical box in each direction.
        RealBox real_box(domain_lo.dataPtr(), domain_hi.dataPtr());

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    const Real* dx = geom.CellSize();

    // Nghost = number of ghost cells for each array
    int Nghost = NUM_GHOST_CELLS;

    // Ncomp = number of components for each array
    int Ncomp  = Idx::NumScalars;

    Vector<BCRec> state_bc(Ncomp);
    for (int n = 0; n < Ncomp; ++n)
    {
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
        {
            // is_periodic overrides inputs in domain_(lo/hi)_bc_type
            if (geom.isPeriodic(i))
            {
                state_bc[n].setLo(i, BCType::int_dir);
                state_bc[n].setHi(i, BCType::int_dir);
            }
            else
            {
                state_bc[n].setLo(i, domain_lo_bc_types[i]);
                state_bc[n].setHi(i, domain_hi_bc_types[i]);
            }
        }
    }

    // Initialize variable names
    Variable::Initialize();

    // Initialize diagnostics names
    Diagnostics::Initialize();

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate two state multifabs; one will store the old state, the other the new.
    MultiFab state_old(ba, dm, Ncomp, Nghost);
    MultiFab state_new(ba, dm, Ncomp, Nghost);
    MultiFab diagnostics(ba, dm, Diag::NumScalars, Nghost);

    // time = starting time in the simulation
    Real time = 0.0;
    
    
    
    //////////////////////////////////////
    
    
    
    // Initialize state_new by calling a C++ initializing routine.
    if (initialize_from_data)
    {
        VisMF::Read(state_new, amrex::MultiFabFileFullPrefix(0, "End_data_Final", "Level_", "Cell"));
    }
    else 
        init(state_new, time, geom);

    // Fill ghost cells for state_new from interior & periodic BCs
    state_new.FillBoundary(geom.periodicity());
    FillDomainBoundary(state_new, geom, state_bc);

    // Write diagnostics for step 0
    {
        // Fill the diagnostics data
        fill_state_diagnostics(diagnostics, state_new, geom);

        // Write a plotfile of the diagnostic data
        const std::string& pltfile = amrex::Concatenate("diag_plt",0,7);
        WriteSingleLevelPlotfile(pltfile, diagnostics, Diagnostics::names, geom, 0.0, 0);
    }

    // Compute the time step
    Real dt;
    
    if (elliptic)
    {
        dt = cfl*dx[0]*dx[0];
#if AMREX_SPACEDIM > 1
        dt = amrex::min(dt, cfl*dx[1]*dx[1]);
#endif
#if AMREX_SPACEDIM > 2
        dt = amrex::min(dt, cfl*dx[2]*dx[2]);
#endif
    }
    else
    {
        dt = cfl*dx[0];
#if AMREX_SPACEDIM > 1
        dt = amrex::min(dt, cfl*dx[1]);
#endif
#if AMREX_SPACEDIM > 2
        dt = amrex::min(dt, cfl*dx[2]);
#endif
    }

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
        // Call user function to rescale state
        post_update(S_data, geom);

        // Fill ghost cells for S_data from interior & periodic BCs
        S_data.FillBoundary(geom.periodicity());
        FillDomainBoundary(S_data, geom, state_bc);
    };

    // Create a post-timestep function
    auto post_timestep_fun = [&](){
        // Tell the I/O Processor to write out which step we're doing
        const int n = integrator.get_step_number();
        amrex::Print() << "Advanced step " << n << "\n";

        // Write a plotfile of the current data
        if (plot_int > 0 && n % plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",n,7);
            WriteSingleLevelPlotfile(pltfile, integrator.get_new_data(), Variable::names, geom, integrator.get_time(), n);
        }

        // Calculate diagnostics and write to plotfile
        if (diag_int > 0 && n % diag_int == 0)
        {
            // Fill the diagnostics data
            fill_state_diagnostics(diagnostics, integrator.get_new_data(), geom);

            // Write a plotfile of the diagnostic data
            const std::string& pltfile = amrex::Concatenate("diag_plt",n,7);
            WriteSingleLevelPlotfile(pltfile, diagnostics, Diagnostics::names, geom, integrator.get_time(), n);
        }
        if (write_checkpoint && checkpoint_int > 0 && n % checkpoint_int == 0)
        {
            const std::string& checkpointname = amrex::Concatenate("End_data",n,7);
   
            amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

            const int nlevels = 1;

            bool callBarrier = true;

            // ---- prebuild a hierarchy of directories
            // ---- dirName is built first.  if dirName exists, it is renamed.  then build
            // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
            // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
            // ---- after all directories are built
            // ---- ParallelDescriptor::IOProcessor() creates the directories
            amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, callBarrier);
        
            VisMF::Write(diagnostics, amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "Cell"));
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
        
        if (write_checkpoint)
        {
            const std::string& checkpointname = "End_data_Final";

            amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

            const int nlevels = 1;

            bool callBarrier = true;

            // ---- prebuild a hierarchy of directories
            // ---- dirName is built first.  if dirName exists, it is renamed.  then build
            // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
            // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
            // ---- after all directories are built
            // ---- ParallelDescriptor::IOProcessor() creates the directories
            amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, callBarrier);
        
            VisMF::Write(diagnostics, amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "Cell"));
        }
        
    }

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real stop_time = amrex::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
}

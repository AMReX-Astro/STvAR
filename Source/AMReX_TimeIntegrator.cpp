#include AMReX_TimeIntegrator.H

using namespace amrex;

TimeIntegrator::TimeIntegrator(const BoxArray& ba, const DistributionMapping& dm, const int Ncomp, const int Nghost)
{
    ParmParse pp("integration");

    // Read an integrator type, if not recognized, then read weights/nodes/butcher tableau

    // Read weights/nodes/butcher tableau
    pp.queryarr("weights", weights);
    pp.queryarr("nodes", nodes);

    Vector<Real> btable; // flattened into row major format
    pp.queryarr("tableau", btable);

    // Sanity check the inputs
    if (weights.size() != nodes.size())
    {
        Error("integration.weights should be the same length as integration.nodes");
    } else {
        number_nodes = weights.size();
        const int nTableau = (number_nodes * (number_nodes + 1)) / 2; // includes diagonal
        if (btable.size() != nTableau)
        {
            Error("integration.tableau incorrect length - should include the Butcher Tableau diagonal.");
        }
    }

    // Fill tableau from the flattened entries
    int k = 0;
    for (int i = 0; i < number_nodes; ++i)
    {
        Vector<Real> stage_row;
        for (int j = 0; j < i; ++j)
        {
            stage_row.push_back(btable[k]);
            ++k;
        }

        tableau.push_back(stage_row);
    }

    // Check that this is an explicit method
    for (const auto& astage : tableau)
    {
        if (astage[-1] != 0.0)
        {
            Error("TimeIntegrator currently only supports explicit Butcher tableaus.");
        }
    }

    // Create MultiFabs for stages
    for (int i = 0; i < number_nodes; ++i)
    {
        F_nodes.emplace_back(ba, dm, Ncomp, Nghost);
    }

    // Create MultiFabs for solution
    for (int i = 0; i < 3; ++i)
    {
        S_val.emplace_back(ba, dm, Ncomp, Nghost);
    }
}

amrex::MultiFab& TimeIntegrator::get_old()
{
    if (S_val.size() == 0)
    {
        Error("Time integrator has not been initialized.")
    }

    return S_val[StateTimes::Old];
}

amrex::MultiFab& TimeIntegrator::get_new()
{
    if (S_val.size() == 0)
    {
        Error("Time integrator has not been initialized.")
    }

    return S_val[StateTimes::New];
}

amrex::MultiFab& TimeIntegrator::get_tmp()
{
    if (S_val.size() == 0)
    {
        Error("Time integrator has not been initialized.")
    }

    return S_val[StateTimes::Tmp];
}

Real TimeIntegrator::advance(std::function<void(MultiFab&, const MultiFab&, const Real)> F, const Real time, const Real timestep)
{
    // References to our state solution data
    MultiFab& S_old = get_old();
    MultiFab& S_tmp = get_tmp();
    MultiFab& S_new = get_new();

    // Assume before advance() that S_new is valid data at the current time ("time" argument)
    // So we update S_old by copying the current state.
    MultiFab::Copy(S_old, S_new, 0, 0, S_old.nComp(), S_old.nGrow());

    // Fill the RHS F_nodes at each stage
    for (int i = 0; i < number_nodes; ++i)
    {
        // Get current stage time, t = t_old + h * Ci
        Real stage_time = time + timestep * nodes[i];

        // Fill S_tmp with the solution value for evaluating F at the current stage
        // Copy S_tmp = S_old
        MultiFab::Copy(S_tmp, S_old, 0, 0, S_old.nComp(), S_old.nGrow());
        if (i > 0) {
            // Saxpy across the tableau row:
            // S_tmp += h * Aij * Fj
            // We should fuse these kernels ...
            for (int j = 0; j < i; ++j)
            {
                MultiFab::Saxpy(S_tmp, timestep * tableau[i][j], F_nodes[j], 0, 0, S_tmp.nComp(), S_tmp.nGrow());
            }
        }

        // Fill F[i], the RHS at the current stage
        // F[i] = RHS(y, t) at y = S_tmp, t = stage_time
        F(F_nodes[i], S_tmp, stage_time);
    }

    // Fill new State. S_new = S_old already, so we add the stage contributions.
    // Saxpy S_new += h * Wi * Fi for integration weights Wi
    // We should fuse these kernels ...
    for (int i = 0; i < number_nodes; ++i)
    {
        MultiFab::Saxpy(S_new, timestep * weights[i], F_nodes[i], 0, 0, S_new.nComp(), S_new.nGrow());
    }

    // If we are working with an extended Butcher tableau, we can estimate the error in S_tmp here,
    // and then calculate an adaptive timestep.

    // Return timestep
    return timestep;
}
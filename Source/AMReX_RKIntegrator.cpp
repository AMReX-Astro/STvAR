#include "AMReX_RKIntegrator.H"

using namespace amrex;

RKIntegrator::RKIntegrator(amrex::MultiFab& S_old_external,
                           amrex::MultiFab& S_new_external,
                           amrex::Real initial_time) : IntegratorBase(S_old_external, S_new_external, initial_time)
{
    initialize_parameters();
    initialize_stages();

    // Create temporary State MultiFab
    S_tmp_ptr = std::make_unique<MultiFab>(S_old.boxArray(), S_old.DistributionMap(), S_old.nComp(), S_old.nGrow());
}

void RKIntegrator::initialize_parameters()
{
    ParmParse pp("integration.rk");

    // Read an integrator type, if not recognized, then read weights/nodes/butcher tableau
    tableau_type = 0;
    pp.query("type", tableau_type);

    // By default, define no extended weights and no adaptive timestepping
    extended_weights = {};
    use_adaptive_timestep = false;
    pp.query("use_adaptive_timestep", use_adaptive_timestep);

    if (tableau_type == ButcherTableauTypes::User)
    {
        // Read weights/nodes/butcher tableau"
        pp.queryarr("weights", weights);
        pp.queryarr("extended_weights", extended_weights);
        pp.queryarr("nodes", nodes);

        Vector<Real> btable; // flattened into row major format
        pp.queryarr("tableau", btable);

        // Sanity check the inputs
        if (weights.size() != nodes.size())
        {
            Error("integration.rk.weights should be the same length as integration.rk.nodes");
        } else {
            number_nodes = weights.size();
            const int nTableau = (number_nodes * (number_nodes + 1)) / 2; // includes diagonal
            if (btable.size() != nTableau)
            {
                Error("integration.rk.tableau incorrect length - should include the Butcher Tableau diagonal.");
            }
        }

        // Fill tableau from the flattened entries
        int k = 0;
        for (int i = 0; i < number_nodes; ++i)
        {
            Vector<Real> stage_row;
            for (int j = 0; j <= i; ++j)
            {
                stage_row.push_back(btable[k]);
                ++k;
            }

            tableau.push_back(stage_row);
        }

        // Check that this is an explicit method
        for (const auto& astage : tableau)
        {
            if (astage.back() != 0.0)
            {
                Error("RKIntegrator currently only supports explicit Butcher tableaus.");
            }
        }
    } else if (tableau_type > ButcherTableauTypes::User && tableau_type < ButcherTableauTypes::NumTypes)
    {
        initialize_preset_tableau();
    } else {
        Error("RKIntegrator received invalid input for integration.rk.type");
    }
}

void RKIntegrator::initialize_preset_tableau()
{
    switch (tableau_type)
    {
        case ButcherTableauTypes::ForwardEuler:
            nodes = {0.0};
            tableau = {{0.0}};
            weights = {1.0};
            break;
        case ButcherTableauTypes::Trapezoid:
            nodes = {0.0,
                     1.0};
            tableau = {{0.0},
                       {1.0, 0.0}};
            weights = {0.5, 0.5};
            break;
        case ButcherTableauTypes::SSPRK3:
            nodes = {0.0,
                     1.0,
                     0.5};
            tableau = {{0.0},
                       {1.0, 0.0},
                       {0.25, 0.25, 0.0}};
            weights = {1./6., 1./6., 2./3.};
            break;
        case ButcherTableauTypes::RK4:
            nodes = {0.0,
                     0.5,
                     0.5,
                     1.0};
            tableau = {{0.0},
                       {0.5, 0.0},
                       {0.0, 0.5, 0.0},
                       {0.0, 0.0, 1.0, 0.0}};
            weights = {1./6., 1./3., 1./3., 1./6.};
            break;
        default:
            Error("Invalid RK Integrator tableau type");
            break;
    }

    number_nodes = weights.size();
}

void RKIntegrator::initialize_stages()
{
    // Create MultiFabs for stages
    for (int i = 0; i < number_nodes; ++i)
    {
        F_nodes.emplace_back(S_old.boxArray(), S_old.DistributionMap(), S_old.nComp(), S_old.nGrow());
    }
}

Real RKIntegrator::advance(const Real timestep)
{
    // Get a reference to our temporary State workspace
    auto& S_tmp = *S_tmp_ptr;

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

            // Call the post-update hook for S_tmp
            post_update(S_tmp);
        }

        // Fill F[i], the RHS at the current stage
        // F[i] = RHS(y, t) at y = S_tmp, t = stage_time
        rhs(F_nodes[i], S_tmp, stage_time);
    }

    // Fill new State. S_new = S_old already, so we add the stage contributions.
    // Saxpy S_new += h * Wi * Fi for integration weights Wi
    // We should fuse these kernels ...
    for (int i = 0; i < number_nodes; ++i)
    {
        MultiFab::Saxpy(S_new, timestep * weights[i], F_nodes[i], 0, 0, S_new.nComp(), S_new.nGrow());
    }

    // Call the post-update hook for S_new
    post_update(S_new);

    // Update time
    time += timestep;

    // If we are working with an extended Butcher tableau, we can estimate the error in S_tmp here,
    // and then calculate an adaptive timestep.

    // Return timestep
    return timestep;
}

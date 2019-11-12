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
        case ButcherTableauTypes::Ralston4:
            nodes = {0.0,
                     0.4,
                     0.45573725,
                     1.0};
            tableau = {{0.0},
                       {0.4, 0.0},
                       {0.29697761, 0.15875964, 0.0},
                       {0.21810040, -3.05096516, 3.83286476, 0.0}};
            weights = {0.17476028, -0.55148066, 1.20553560, 0.17118478};
            break;
        case ButcherTableauTypes::HeunEuler21:
            nodes = {0.0,
                     1.0};
            tableau = {{0.0},
                       {1.0, 0.0}};
            weights = {0.5, 0.5};
            extended_weights = {1.0, 0.0};
        case ButcherTableauTypes::Fehlberg21:
            nodes = {0.0,
                     0.5,
                     1.0};
            tableau = {{0.0},
                       {0.5, 0.0},
                       {1./256., 255./256., 0.0}};
            weights = {1./256., 255./256., 0.0};
            extended_weights = {1./512., 255./256., 1./512.};
        case ButcherTableauTypes::BogackiShampine32:
            nodes = {0.0,
                     0.5,
                     0.75,
                     1.0};
            tableau = {{0.0},
                       {0.5, 0.0},
                       {0.0, 0.75, 0.0},
                       {2./9., 1./3., 4./9., 0.0}};
            weights = {2./9., 1./3., 4./9., 0.0};
            extended_weights = {7./24., 1./4., 1./3., 1./8.};
        case ButcherTableauTypes::Fehlberg54:
            nodes = {0.0,
                     0.25,
                     3./8.,
                     12./13.,
                     1.0,
                     0.5};
            tableau = {{0.0},
                       {1./4., 0.0},
                       {3./32., 9./32., 0.0},
                       {1932./2197., -7200./2197., 7296./2197., 0.0},
                       {439./216., -8.0, 3680./513., -845./4104., 0.0},
                       {-8./27., 2.0, -3544./2565., 1859./4104., -11./40., 0.0}};
            weights = {16./135., 0.0, 6656./12825., 28561./56430., -9./50., 2./55.};
            extended_weights = {25./216., 0.0, 1408./2565., 2197./4104., -1./5., 0.0};
        case ButcherTableauTypes::CashKarp54:
            nodes = {0.0,
                     1./5.,
                     3./10.,
                     3./5.,
                     1.0,
                     7./8.};
            tableau = {{0.0},
                       {1./5., 0.0},
                       {3./40., 9./40., 0.0},
                       {3./10., -9./10., 6./5., 0.0},
                       {-11./54., 5./2., -70./27., 35./27., 0.0},
                       {1631./55296., 175./512., 575./13824., 44275./110592., 253./4096., 0.0}};
            weights = {37./378., 0.0, 250./621., 125./594., 0.0, 512./1771.};
            extended_weights = {2825./27648., 0.0, 18575./48384., 13525./55296., 277./14336., 1./4.};
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

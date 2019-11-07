#include "AMReX_TimeIntegrator.H"

using namespace amrex;

TimeIntegrator::TimeIntegrator(amrex::MultiFab& S_old_external,
                               amrex::MultiFab& S_new_external,
                               amrex::Real initial_time)
{
    ParmParse pp("integration");

    int integrator_type = IntegratorTypes::ForwardEuler;
    pp.query("type", integrator_type);

    switch (integrator_type)
    {
        case IntegratorTypes::ForwardEuler:
            integrator_ptr = std::make_unique<FEIntegrator>(S_old_external, S_new_external, initial_time);
            break;
        case IntegratorTypes::ExplicitRungeKutta:
            integrator_ptr = std::make_unique<RKIntegrator>(S_old_external, S_new_external, initial_time);
            break;
        default:
            Error("integration.type did not name a valid integrator type.");
            break;
    }

    // Default nsteps to 10, allow us to set it to something else in the inputs file
    max_num_steps = 10;
    pp.query("nsteps", max_num_steps);

    // Stopping criteria
    end_time = 1.0;
    pp.query("end_time", end_time);

    // By default, do nothing post-timestep
    set_post_timestep([](){});

    // By default, set the RHS to 0.0
    set_rhs([](MultiFab& S_rhs, const MultiFab& S_data, const Real time){ S_rhs = 0.0; });

    // Set the initial time and step number
    time = initial_time;
    step_number = 0;
}

void TimeIntegrator::integrate(const amrex::Real start_timestep)
{
    Real timestep = start_timestep;
    bool stop_advance = false;
    for (step_number = 1; step_number <= max_num_steps && !stop_advance; ++step_number)
    {
        if (end_time - time < timestep) {
            timestep = end_time - time;
            stop_advance = true;
        }

        // Call the time integrator advance
        integrator_ptr->advance(timestep);

        // Update our time variable
        time = integrator_ptr->get_time();

        // Call the post-timestep hook
        post_timestep();
    }
}

void TimeIntegrator::set_post_timestep(std::function<void ()> F)
{
    post_timestep = F;
}

void TimeIntegrator::set_rhs(std::function<void (MultiFab& S_rhs, const MultiFab& S_data, const Real time)> F)
{
    integrator_ptr->set_rhs(F);
}

amrex::MultiFab& TimeIntegrator::get_new_data()
{
    return integrator_ptr->get_new_data();
}

amrex::MultiFab& TimeIntegrator::get_old_data()
{
    return integrator_ptr->get_old_data();
}

amrex::Real TimeIntegrator::get_time()
{
    return time;
}

amrex::Real TimeIntegrator::get_timestep()
{
    return timestep;
}

int TimeIntegrator::get_step_number()
{
    return step_number;
}

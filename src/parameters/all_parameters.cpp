#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include "parameters/all_parameters.h"

namespace PHiLiP {
namespace Parameters {

AllParameters::AllParameters ()
    : manufactured_convergence_study_param(ManufacturedConvergenceStudyParam())
    , ode_solver_param(ODESolverParam())
    , linear_solver_param(LinearSolverParam())
    , euler_param(EulerParam())
    , pcout(std::cout, dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
{ }
void AllParameters::declare_parameters (dealii::ParameterHandler &prm)
{
    const int mpi_rank = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
    dealii::ConditionalOStream pcout(std::cout, mpi_rank==0);
    pcout << "Declaring inputs." << std::endl;
    prm.declare_entry("dimension", "1",
                      dealii::Patterns::Integer(),
                      "Number of dimensions");

    prm.declare_entry("use_weak_form", "true",
                      dealii::Patterns::Bool(),
                      "Use weak form by default. If false, use strong form.");

    prm.declare_entry("use_collocated_nodes", "false",
                      dealii::Patterns::Bool(),
                      "Use Gauss-Legendre by default. Otherwise, use Gauss-Lobatto to collocate.");

    prm.declare_entry("use_split_form", "false",
                      dealii::Patterns::Bool(),
                      "Use original form by defualt. Otherwise, split the fluxes.");

    prm.declare_entry("use_periodic_bc", "false",
                      dealii::Patterns::Bool(),
                      "Use other boundary conditions by default. Otherwise use periodic (for 1d burgers only");

    prm.declare_entry("test_type", "run_control",
                      dealii::Patterns::Selection(
                      " run_control | "
                      " burgers_energy_stability | "
                      " euler_gaussian_bump | "
                      " euler_cylinder | "
                      " euler_vortex | "
                      " euler_entropy_waves | "
                      " numerical_flux_convervation | "
                      " jacobian_regression |"
                      " advection_periodicity |"
                      " euler_split_taylor_green |"
                      " collocation"),
                      "The type of test we want to solve. "
                      "Choices are (only run control has been coded up for now)" 
                      " <run_control | " 
                      "  burgers_energy_stability | "
                      "  euler_gaussian_bump | "
                      "  euler_cylinder | "
                      "  euler_vortex | "
                      "  euler_entropy_waves | "
                      "  numerical_flux_convervation | "
                      "  jacobian_regression |"
					  "  euler_split_taylor_green |"
					  "  advection_periodicity |"
					  "  euler_split_taylor_green |"
					  "  collocation >.");

    prm.declare_entry("pde_type", "advection",
                      dealii::Patterns::Selection(
                          " advection | "
                          " diffusion | "
                          " convection_diffusion | "
                          " advection_vector | "
                          " burgers_inviscid | "
                          " euler |"
                          " mhd"),
                      "The PDE we want to solve. "
                      "Choices are " 
                      " <advection | " 
                      "  diffusion | "
                      "  convection_diffusion | "
                      "  advection_vector | "
                      "  burgers_inviscid | "
                      "  euler | "
                      "  mhd>.");
    prm.declare_entry("conv_num_flux", "lax_friedrichs",
                      dealii::Patterns::Selection("lax_friedrichs | roe | split_form"),
                      "Convective numerical flux. "
                      "Choices are <lax_friedrichs | roe | split_form>.");

    prm.declare_entry("diss_num_flux", "symm_internal_penalty",
                      dealii::Patterns::Selection("symm_internal_penalty"),
                      "Dissipative numerical flux. "
                      "Choices are <symm_internal_penalty>.");

    Parameters::LinearSolverParam::declare_parameters (prm);
    Parameters::ManufacturedConvergenceStudyParam::declare_parameters (prm);
    Parameters::ODESolverParam::declare_parameters (prm);

    Parameters::EulerParam::declare_parameters (prm);

    pcout << "Done declaring inputs." << std::endl;
}

void AllParameters::parse_parameters (dealii::ParameterHandler &prm)
{
    pcout << "Parsing main input..." << std::endl;

    dimension                   = prm.get_integer("dimension");

    const std::string test_string = prm.get("test_type");
    if (test_string == "run_control") { test_type = run_control; }
    else if (test_string == "burgers_energy_stability") { test_type = burgers_energy_stability; }
    else if (test_string == "euler_gaussian_bump") { test_type = euler_gaussian_bump; }
    else if (test_string == "euler_cylinder") { test_type = euler_cylinder; }
    else if (test_string == "euler_vortex") { test_type = euler_vortex; }
    else if (test_string == "euler_entropy_waves") { test_type = euler_entropy_waves; }
    else if (test_string == "numerical_flux_convervation") { test_type = numerical_flux_convervation; }
    else if (test_string == "jacobian_regression") { test_type = jacobian_regression; }
    else if (test_string == "advection_periodicity") {test_type = advection_periodicity; }
    else if (test_string == "euler_split_taylor_green") {test_type = euler_split_taylor_green;}

    const std::string pde_string = prm.get("pde_type");
    if (pde_string == "advection") {
        pde_type = advection;
        nstate = 1;
    } else if (pde_string == "advection_vector") {
        pde_type = advection_vector;
        nstate = 2;
    } else if (pde_string == "diffusion") {
        pde_type = diffusion;
        nstate = 1;
    } else if (pde_string == "convection_diffusion") {
        pde_type = convection_diffusion;
        nstate = 1;
    } else if (pde_string == "burgers_inviscid") {
        pde_type = burgers_inviscid;
        nstate = dimension;
    } else if (pde_string == "euler") {
        pde_type = euler;
        nstate = dimension+2;
    }

    use_weak_form = prm.get_bool("use_weak_form");
    use_collocated_nodes = prm.get_bool("use_collocated_nodes");
    use_split_form = prm.get_bool("use_split_form");
    use_periodic_bc = prm.get_bool("use_periodic_bc");

    const std::string conv_num_flux_string = prm.get("conv_num_flux");
    if (conv_num_flux_string == "lax_friedrichs") conv_num_flux_type = lax_friedrichs;
    if (conv_num_flux_string == "split_form") conv_num_flux_type = split_form;
    if (conv_num_flux_string == "roe") conv_num_flux_type = roe;

    const std::string diss_num_flux_string = prm.get("diss_num_flux");
    if (diss_num_flux_string == "symm_internal_penalty") diss_num_flux_type = symm_internal_penalty;


    pcout << "Parsing linear solver subsection..." << std::endl;
    linear_solver_param.parse_parameters (prm);

    pcout << "Parsing ODE solver subsection..." << std::endl;
    ode_solver_param.parse_parameters (prm);

    pcout << "Parsing manufactured convergence study subsection..." << std::endl;
    manufactured_convergence_study_param.parse_parameters (prm);

    pcout << "Parsing euler subsection..." << std::endl;
    euler_param.parse_parameters (prm);

    pcout << "Done parsing." << std::endl;
}

} // Parameters namespace
} // PHiLiP namespace

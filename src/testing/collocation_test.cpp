#include "collocation_test.h"

#define TOLERANCE 1E-12

namespace PHiLiP {
namespace Tests {

template <int dim, int nstate>
CollocationTest<dim, nstate>::CollocationTest(const Parameters::AllParameters *const parameters_input)
:
TestsBase::TestsBase(parameters_input)
, mpi_communicator(MPI_COMM_WORLD)
, pcout(std::cout, dealii::Utilities::MPI::this_mpi_process(mpi_communicator)==0)
{}

template <int dim, int nstate>
int CollocationTest<dim, nstate>::run_test() const
{
    //dealii::Triangulation<dim> grid;
#if PHILIP_DIM==1 // dealii::parallel::distributed::Triangulation<dim> does not work for 1D
        dealii::Triangulation<dim> grid;
#else
    dealii::parallel::distributed::Triangulation<dim> grid(this->mpi_communicator);
#endif
    double left = 0.0;
    double right = 2.0;
    const bool colorize = true;
    int n_refinements = 4;
    unsigned int poly_degree = 3;
    dealii::GridGenerator::hyper_cube(grid, left, right, colorize);

    grid.refine_global(n_refinements);

    std::shared_ptr < PHiLiP::DGBase<dim, double> > dg = PHiLiP::DGFactory<dim,double>::create_discontinuous_galerkin(all_parameters, poly_degree, &grid);
    //dg->allocate_system ();
    dg->evaluate_mass_matrices(false);

    pcout << "EVALUATING COLLOCATION" << std::endl;

    for (unsigned int i = 0; i < dg->global_mass_matrix.m(); ++i)
    {
        for (unsigned int j = 0; j < dg->global_mass_matrix.n   (); ++j)
        {
            if (i != j)
            {
                if (dg->global_mass_matrix.el(i,j) >= TOLERANCE)
                {
                    pcout << "matrix element (" << i <<", " << j <<") is not zero. It is " << dg->global_mass_matrix.el(i,j) << std::endl;
                    return 1;
                    break;
                }
            }
        }
    }
    return 0;
}

//#if PHILIP_DIM==3
    template class CollocationTest <PHILIP_DIM,1>;
    template class CollocationTest <PHILIP_DIM,2>;
    template class CollocationTest <PHILIP_DIM,3>;
    template class CollocationTest <PHILIP_DIM,4>;
    template class CollocationTest <PHILIP_DIM,5>;
}
}




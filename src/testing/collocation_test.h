
#ifndef _COLLOCATION_TEST_H_
#define _COLLOCATION_TEST_H_

#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/fe/mapping_q.h>
#include "tests.h"


#include "parameters/all_parameters.h"
#include "parameters/parameters.h"
#include "numerical_flux/numerical_flux.h"
#include "physics/physics_factory.h"
#include "physics/physics.h"
#include "dg/dg.h"
#include "ode_solver/ode_solver.h"

#include<fenv.h>

namespace PHiLiP {
namespace Tests {

template <int dim, int nstate>
class CollocationTest : public TestsBase
{
public:
    CollocationTest() = delete;
    CollocationTest(const Parameters::AllParameters *const parameters_input);
    int run_test() const override;

private:
    const MPI_Comm mpi_communicator;
    dealii::ConditionalOStream pcout;
};

}
}



#endif

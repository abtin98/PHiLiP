#include "split_form.h"

namespace PHiLiP {
namespace splitform {

template <int dim, typename real>
SplitForm<dim, real>::SplitForm(Parameters::AllParameters &parameters, int state_numbers) :
prm(parameters),
nstate (state_numbers)
{

}

} //splitform namespace
} //PHiLiP namespace

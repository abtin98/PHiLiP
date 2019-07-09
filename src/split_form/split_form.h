#ifndef SPLITFORM
#define SPLITFORM

#include <iostream>
#include <vector>
#include "all_parameters.h"

namespace PHiLiP {
namespace splitform {
template <int dim, typename real>
class SplitForm
{
	std::vector<double> alpha;
	std::vector<std::string> f;
	std::vector<std::string> g;
	Parameters::AllParameters prm;
	int nstate;

	SplitForm() = delete;
	SplitForm(Parameters::AllParameters &parameters, int state_numbers);
};
} //splitform namespace
} //PHiLiP namespace
#endif

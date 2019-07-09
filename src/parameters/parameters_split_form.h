#ifndef __PARAMETERS_SPLIT_FORM_H__
#define __PARAMETERS_SPLIT_FORM_H__

#include <deal.II/base/parameter_handler.h>
#include "parameters/parameters.h"

namespace PHiLiP {
namespace Parameters {

class SplitFormParam
{
	SplitFormParam ();


	/// Declares the possible variables and sets the defaults.
	static void declare_parameters (dealii::ParameterHandler &prm);
	/// Parses input file and sets the variables.
	void parse_parameters (dealii::ParameterHandler &prm);
};

} //Parameters namespace
} //PHiLiP namespace

#endif

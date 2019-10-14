#ifndef __MHD__
#define __MHD__

#include <deal.II/base/tensor.h>
#include "physics.h"

namespace PHiLiP {
namespace Physics {


template <int dim, int nstate, typename real>
class MHD : public PhysicsBase <dim, nstate, real>
{
public:
    /// Constructor
    MHD (const double gamma_gas)
    : gam(gamma_gas)
    , gamm1(gam-1.0)
    {
        static_assert(nstate==8, "Physics::MHD() should be created with nstate=8");

    };
    /// Destructor
    ~MHD ()
    {};

    /// Constant heat capacity ratio of fluid
    const double gam;
    /// Gamma-1.0 used often
    const double gamm1;

    std::array<real,nstate> manufactured_solution (const dealii::Point<dim,double> &pos) const;

    /// Dissipative flux: 0
    std::array<dealii::Tensor<1,dim,real>,nstate> dissipative_flux (
        const std::array<real,nstate> &conservative_soln,
        const std::array<dealii::Tensor<1,dim,real>,nstate> &solution_gradient) const;

    /// Source term is zero or depends on manufactured solution
    std::array<real,nstate> source_term (
        const dealii::Point<dim,double> &pos,
        const std::array<real,nstate> &conservative_soln) const;

    /// Given conservative variables [density, [momentum], total energy],
    /// returns primitive variables [density, [velocities], pressure].
    ///
    /// Opposite of convert_primitive_to_conservative
    std::array<real,nstate> convert_conservative_to_primitive ( const std::array<real,nstate> &conservative_soln ) const;

    /// Given primitive variables [density, [velocities], pressure],
    /// returns conservative variables [density, [momentum], total energy].
    ///
    /// Opposite of convert_primitive_to_conservative
    std::array<real,nstate> convert_primitive_to_conservative ( const std::array<real,nstate> &primitive_soln ) const;

    /// Evaluate pressure from conservative variables
    real compute_pressure ( const std::array<real,nstate> &conservative_soln ) const;

    dealii::Tensor<1,dim,real> compute_velocities ( const std::array<real,nstate> &conservative_soln ) const;

    real compute_velocity_squared (const dealii::Tensor<1,dim,real> &velocities) const;

    dealii::Tensor<1,dim,real> extract_velocities_from_primitive ( const std::array<real,nstate> &primitive_soln ) const;

    real compute_total_energy (const std::array<real,nstate> &primitive_soln) const;


    /// Evaluate Magnetic Energy
    real compute_magnetic_energy (const std::array<real,nstate> &conservative_or_primitive_soln) const;



    std::array<real,nstate> convective_eigenvalues (
        const std::array<real,nstate> &conservative_soln,
        const dealii::Tensor<1,dim,real> &normal) const;

    real max_convective_eigenvalue (const std::array<real,nstate> &conservative_soln) const;

    std::array<dealii::Tensor<1,dim,real>,nstate> convective_flux (const std::array<real,nstate> &conservative_soln) const;

    std::array<dealii::Tensor<1,dim,real>,nstate> convective_numerical_split_flux(const std::array<real,nstate> &soln_const,
                                       const std::array<real,nstate> &soln_loop) const;
protected:


};

} // Physics namespace
} // PHiLiP namespace

#endif


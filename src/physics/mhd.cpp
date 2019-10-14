#include <cmath>
#include <vector>

#include <Sacado.hpp>
#include <deal.II/differentiation/ad/sacado_math.h>
#include <deal.II/differentiation/ad/sacado_number_types.h>
#include <deal.II/differentiation/ad/sacado_product_types.h>

#include "physics.h"
#include "mhd.h"


namespace PHiLiP {
namespace Physics {

template <int dim, int nstate, typename real>
std::array<real,nstate> MHD<dim,nstate,real>
::source_term (
    const dealii::Point<dim,double> &/*pos*/,
    const std::array<real,nstate> &/*conservative_soln*/) const
{
    std::array<real,nstate> source_term;
    for (int s=0; s<nstate; s++) {
        source_term[s] = 0;
    }

    return source_term;
}

//incomplete
template <int dim, int nstate, typename real>
inline std::array<real,nstate> MHD<dim,nstate,real>
::convert_conservative_to_primitive ( const std::array<real,nstate> &conservative_soln ) const
{
    std::array<real, nstate> primitive_soln;

    real density = conservative_soln[0];
    dealii::Tensor<1,dim,real> vel = compute_velocities (conservative_soln);
    real pressure = 0;//compute_pressure (conservative_soln);


    primitive_soln[0] = density;
    for (unsigned int d = 0; d<dim; ++d) {
        primitive_soln[1+d] = vel[d];
        //magnetic field
        primitive_soln[5+d] = conservative_soln[5+d];
    }
    primitive_soln[4] = pressure;

    return primitive_soln;
}

//incomplete
template <int dim, int nstate, typename real>
inline std::array<real,nstate> MHD<dim,nstate,real>
::convert_primitive_to_conservative ( const std::array<real,nstate> &primitive_soln ) const
{
    const real density = primitive_soln[0];
    const dealii::Tensor<1,dim,real> velocities = extract_velocities_from_primitive(primitive_soln);

    std::array<real, nstate> conservative_soln;
    conservative_soln[0] = density;
    for (int d=0; d<dim; ++d) {
        conservative_soln[1+d] = density*velocities[d];

        //magnetic field
        conservative_soln[5+d] = primitive_soln[5+d];
    }
    conservative_soln[4] = compute_total_energy(primitive_soln);

    return conservative_soln;
}

template <int dim, int nstate, typename real>
inline dealii::Tensor<1,dim,real> MHD<dim,nstate,real>
::compute_velocities ( const std::array<real,nstate> &conservative_soln ) const
{
    const real density = conservative_soln[0];
    dealii::Tensor<1,dim,real> vel;
    for (int d=0; d<dim; ++d) { vel[d] = conservative_soln[1+d]/density; }
    return vel;
}

template <int dim, int nstate, typename real>
inline real MHD<dim,nstate,real>
::compute_velocity_squared ( const dealii::Tensor<1,dim,real> &velocities ) const
{
    real vel2 = 0.0;
    for (int d=0; d<dim; d++) { vel2 = vel2 + velocities[d]*velocities[d]; }
    return vel2;
}

template <int dim, int nstate, typename real>
inline dealii::Tensor<1,dim,real> MHD<dim,nstate,real>
::extract_velocities_from_primitive ( const std::array<real,nstate> &primitive_soln ) const
{
    dealii::Tensor<1,dim,real> velocities;
    for (int d=0; d<dim; d++) { velocities[d] = primitive_soln[1+d]; }
    return velocities;
}

template <int dim, int nstate, typename real>
inline real MHD<dim,nstate,real>
::compute_total_energy ( const std::array<real,nstate> &primitive_soln ) const
{
    const real density = primitive_soln[0];
    const real pressure = primitive_soln[nstate-1];
    const dealii::Tensor<1,dim,real> velocities = extract_velocities_from_primitive(primitive_soln);
    const real vel2 = compute_velocity_squared(velocities);
    const real magnetic_energy = compute_magnetic_energy(primitive_soln);
    const real tot_energy = pressure / gamm1 + 0.5*density*vel2 + magnetic_energy;
    return tot_energy;
}

template <int dim, int nstate, typename real>
inline real MHD<dim,nstate,real>
::compute_magnetic_energy (const std::array<real,nstate> &conservative_or_primitive_soln) const
{
    real magnetic_energy = 0;
    for (int i = 1; i <= 3; ++i)
        magnetic_energy += 1./2. * (conservative_or_primitive_soln[nstate - i] * conservative_or_primitive_soln[nstate - i] );
    return magnetic_energy;
}

template <int dim, int nstate, typename real>
std::array<dealii::Tensor<1,dim,real>,nstate> MHD<dim,nstate,real>
::dissipative_flux (
    const std::array<real,nstate> &/*conservative_soln*/,
    const std::array<dealii::Tensor<1,dim,real>,nstate> &/*solution_gradient*/) const
{
    std::array<dealii::Tensor<1,dim,real>,nstate> diss_flux;
    // No dissipation
    for (int i=0; i<nstate; i++) {
        diss_flux[i] = 0;
    }
    return diss_flux;
}

template <int dim, int nstate, typename real>
std::array<real,nstate> MHD<dim,nstate,real>
::convective_eigenvalues (
    const std::array<real,nstate> &conservative_soln,
    const dealii::Tensor<1,dim,real> &normal) const
{
    const dealii::Tensor<1,dim,real> vel = compute_velocities(conservative_soln);
    std::array<real,nstate> eig;
    real vel_dot_n = 0.0;
    for (int d=0;d<dim;++d) { vel_dot_n += vel[d]*normal[d]; };
    for (int i=0; i<nstate; i++) {
        eig[i] = vel_dot_n;
        //eig[i] = advection_speed*normal;

        //eig[i] = 1.0;
        //eig[i] = -1.0;
    }
    return eig;
}
template <int dim, int nstate, typename real>
real MHD<dim,nstate,real>
::max_convective_eigenvalue (const std::array<real,nstate> &conservative_soln) const
{
    (void) conservative_soln;

    const real max_eig = 0;//sqrt(vel2) + sound;
    //std::cout << "max eig calculated" << std::endl;

    return max_eig;
}

template <int dim, int nstate, typename real>
std::array<dealii::Tensor<1,dim,real>,nstate> MHD<dim,nstate,real>
::convective_flux (const std::array<real,nstate> &conservative_soln) const
{
    std::array<dealii::Tensor<1,dim,real>,nstate> conv_flux;
    const real density = conservative_soln[0];
    const real pressure = 0;//compute_pressure (conservative_soln);
    const dealii::Tensor<1,dim,real> vel = compute_velocities(conservative_soln);
    const real specific_total_energy = conservative_soln[nstate-1]/conservative_soln[0];
    const real specific_total_enthalpy = specific_total_energy + pressure/density;

    for (int flux_dim=0; flux_dim<dim; ++flux_dim) {
        // Density equation
        conv_flux[0][flux_dim] = conservative_soln[1+flux_dim];
        // Momentum equation
        for (int velocity_dim=0; velocity_dim<dim; ++velocity_dim){
            conv_flux[1+velocity_dim][flux_dim] = density*vel[flux_dim]*vel[velocity_dim];
        }
        conv_flux[1+flux_dim][flux_dim] += pressure; // Add diagonal of pressure
        // Energy equation
        conv_flux[nstate-1][flux_dim] = density*vel[flux_dim]*specific_total_enthalpy;
    }
    return conv_flux;
}

template <int dim, int nstate, typename real>
std::array<dealii::Tensor<1,dim,real>,nstate> MHD<dim, nstate, real>
::convective_numerical_split_flux(const std::array<real,nstate> &soln_const,
                                  const std::array<real,nstate> &soln_loop) const
{
    (void) soln_loop;

    return convective_flux(soln_const);
}

// Instantiate explicitly
template class MHD < PHILIP_DIM, 8, double >;
template class MHD < PHILIP_DIM, 8, Sacado::Fad::DFad<double>  >;

} // Physics namespace
} // PHiLiP namespace

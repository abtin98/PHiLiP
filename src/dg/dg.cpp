#include<fstream>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

//#include <deal.II/lac/solver_control.h>
//#include <deal.II/lac/trilinos_precondition.h>
//#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/loop.h>

//#include <deal.II/fe/mapping_q1.h> // Might need mapping_q
#include <deal.II/fe/mapping_q.h> // Might need mapping_q
#include <deal.II/fe/mapping_manifold.h> 
#include <deal.II/fe/mapping_fe_field.h> 

// Finally, we take our exact solution from the library as well as volume_quadrature
// and additional tools.
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_dof_data.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools.templates.h>


#include "dg.h"
#include "post_processor/physics_post_processor.h"

//template class dealii::MappingFEField<PHILIP_DIM,PHILIP_DIM,dealii::LinearAlgebra::distributed::Vector<double>, dealii::hp::DoFHandler<PHILIP_DIM> >;
namespace PHiLiP {
#if PHILIP_DIM==1 // dealii::parallel::distributed::Triangulation<dim> does not work for 1D
    template <int dim> using Triangulation = dealii::Triangulation<dim>;
#else
    template <int dim> using Triangulation = dealii::parallel::distributed::Triangulation<dim>;
#endif

// DGFactory ***********************************************************************
template <int dim, typename real>
std::shared_ptr< DGBase<dim,real> >
DGFactory<dim,real>
::create_discontinuous_galerkin(
    const Parameters::AllParameters *const parameters_input,
    const unsigned int degree,
    Triangulation *const triangulation_input)
{
    using PDE_enum = Parameters::AllParameters::PartialDifferentialEquation;

    PDE_enum pde_type = parameters_input->pde_type;
    if (parameters_input->use_weak_form) {
        if (pde_type == PDE_enum::advection) {
            return std::make_shared< DGWeak<dim,1,real> >(parameters_input, degree, triangulation_input);
        } else if (pde_type == PDE_enum::advection_vector) {
            return std::make_shared< DGWeak<dim,2,real> >(parameters_input, degree, triangulation_input);
        } else if (pde_type == PDE_enum::diffusion) {
            return std::make_shared< DGWeak<dim,1,real> >(parameters_input, degree, triangulation_input);
        } else if (pde_type == PDE_enum::convection_diffusion) {
            return std::make_shared< DGWeak<dim,1,real> >(parameters_input, degree, triangulation_input);
        } else if (pde_type == PDE_enum::burgers_inviscid) {
            return std::make_shared< DGWeak<dim,dim,real> >(parameters_input, degree, triangulation_input);
        } else if (pde_type == PDE_enum::euler) {
            return std::make_shared< DGWeak<dim,dim+2,real> >(parameters_input, degree, triangulation_input);
        }
    } else {
        if (pde_type == PDE_enum::advection) {
            return std::make_shared< DGStrong<dim,1,real> >(parameters_input, degree, triangulation_input);
        } else if (pde_type == PDE_enum::advection_vector) {
            return std::make_shared< DGStrong<dim,2,real> >(parameters_input, degree, triangulation_input);
        } else if (pde_type == PDE_enum::diffusion) {
            return std::make_shared< DGStrong<dim,1,real> >(parameters_input, degree, triangulation_input);
        } else if (pde_type == PDE_enum::convection_diffusion) {
            return std::make_shared< DGStrong<dim,1,real> >(parameters_input, degree, triangulation_input);
        } else if (pde_type == PDE_enum::burgers_inviscid) {
            return std::make_shared< DGStrong<dim,dim,real> >(parameters_input, degree, triangulation_input);
        } else if (pde_type == PDE_enum::euler) {
            return std::make_shared< DGStrong<dim,dim+2,real> >(parameters_input, degree, triangulation_input);
        }
    }
    std::cout << "Can't create DGBase in create_discontinuous_galerkin(). Invalid PDE type: " << pde_type << std::endl;
    return nullptr;
}

// DGBase ***************************************************************************
template <int dim, typename real>
DGBase<dim,real>::DGBase(
    const int nstate_input,
    const Parameters::AllParameters *const parameters_input,
    const unsigned int max_degree_input,
    Triangulation *const triangulation_input)
    : DGBase<dim,real>(nstate_input, parameters_input, max_degree_input, triangulation_input, this->create_collection_tuple(max_degree_input, nstate_input, parameters_input))
{ }

template <int dim, typename real>
DGBase<dim,real>::DGBase( 
    const int nstate_input,
    const Parameters::AllParameters *const parameters_input,
    const unsigned int max_degree_input,
    Triangulation *const triangulation_input,
    const MassiveCollectionTuple collection_tuple)
    : all_parameters(parameters_input)
    , nstate(nstate_input)
    , max_degree(max_degree_input)
    , triangulation(triangulation_input)
    , fe_collection(std::get<0>(collection_tuple))
    , volume_quadrature_collection(std::get<1>(collection_tuple))
    , face_quadrature_collection(std::get<2>(collection_tuple))
    , oned_quadrature_collection(std::get<3>(collection_tuple))
    , fe_collection_lagrange(std::get<4>(collection_tuple))
    , dof_handler(*triangulation)
    , high_order_grid(all_parameters, max_degree_input+1, triangulation)
    , mpi_communicator(MPI_COMM_WORLD)
    , pcout(std::cout, dealii::Utilities::MPI::this_mpi_process(mpi_communicator)==0)
{ 

    dof_handler.initialize(*triangulation, fe_collection);

    set_all_cells_fe_degree (fe_collection.size()-1); // Always sets it to the maximum degree

}

template <int dim, typename real> 
std::tuple<
        //dealii::hp::MappingCollection<dim>, // Mapping
        dealii::hp::FECollection<dim>, // Solution FE
        dealii::hp::QCollection<dim>,  // Volume quadrature
        dealii::hp::QCollection<dim-1>, // Face quadrature
        dealii::hp::QCollection<1>, // 1D quadrature for strong form
        dealii::hp::FECollection<dim> >   // Lagrange polynomials for strong form
DGBase<dim,real>::create_collection_tuple(const unsigned int max_degree, const int nstate, const Parameters::AllParameters *const parameters_input) const
{
    //dealii::hp::MappingCollection<dim> mapping_coll;
    dealii::hp::FECollection<dim>      fe_coll;
    dealii::hp::QCollection<dim>       volume_quad_coll;
    dealii::hp::QCollection<dim-1>     face_quad_coll;
    dealii::hp::QCollection<1>         oned_quad_coll;

    dealii::hp::FECollection<dim>      fe_coll_lagr;

    // for p=0, we use a p=1 FE for collocation, since there's no p=0 quadrature for Gauss Lobatto
    if (parameters_input->use_collocated_nodes==true)
    {
        int degree = 1;
		//const dealii::MappingQ<dim,dim> mapping(degree, true);
		//const dealii::MappingQ<dim,dim> mapping(degree+1, true);
		//const dealii::MappingManifold<dim,dim> mapping;
		//mapping_coll.push_back(mapping);

		const dealii::FE_DGQ<dim> fe_dg(degree);
		const dealii::FESystem<dim,dim> fe_system(fe_dg, nstate);
		fe_coll.push_back (fe_system);

		//

		dealii::Quadrature<1>     oned_quad(degree+1);
		dealii::Quadrature<dim>   volume_quad(degree+1);
		dealii::Quadrature<dim-1> face_quad(degree+1); //removed const

		if (parameters_input->use_collocated_nodes)
		{
			dealii::QGaussLobatto<1> oned_quad_Gauss_Lobatto (degree+1);
			dealii::QGaussLobatto<dim> vol_quad_Gauss_Lobatto (degree+1);
			oned_quad = oned_quad_Gauss_Lobatto;
			volume_quad = vol_quad_Gauss_Lobatto;

			if(dim == 1)
			{
				dealii::QGauss<dim-1> face_quad_Gauss_Legendre (degree+1);
				face_quad = face_quad_Gauss_Legendre;
			}
			else
			{
				dealii::QGaussLobatto<dim-1> face_quad_Gauss_Lobatto (degree+1);
				face_quad = face_quad_Gauss_Lobatto;
			}
		}
		else
		{
			dealii::QGauss<1> oned_quad_Gauss_Legendre (degree+1);
			dealii::QGauss<dim> vol_quad_Gauss_Legendre (degree+1);
			dealii::QGauss<dim-1> face_quad_Gauss_Legendre (degree+1);
			oned_quad = oned_quad_Gauss_Legendre;
			volume_quad = vol_quad_Gauss_Legendre;
			face_quad = face_quad_Gauss_Legendre;
		}
		volume_quad_coll.push_back (volume_quad);
		face_quad_coll.push_back (face_quad);
		oned_quad_coll.push_back (oned_quad);

		dealii::FE_DGQArbitraryNodes<dim,dim> lagrange_poly(oned_quad);
		fe_coll_lagr.push_back (lagrange_poly);
    }

    int minimum_degree = (parameters_input->use_collocated_nodes==true) ?  1 :  0;
    for (unsigned int degree=minimum_degree; degree<=max_degree; ++degree) {
        //const dealii::MappingQ<dim,dim> mapping(degree, true);
        //const dealii::MappingQ<dim,dim> mapping(degree+1, true);
        //const dealii::MappingManifold<dim,dim> mapping;
        //mapping_coll.push_back(mapping);

        // Solution FECollection
        const dealii::FE_DGQ<dim> fe_dg(degree);
        const dealii::FESystem<dim,dim> fe_system(fe_dg, nstate);
        fe_coll.push_back (fe_system);

        dealii::Quadrature<1>     oned_quad(degree+1);
        dealii::Quadrature<dim>   volume_quad(degree+1);
        dealii::Quadrature<dim-1> face_quad(degree+1); //removed const


        if (parameters_input->use_collocated_nodes) {
            dealii::QGaussLobatto<1> oned_quad_Gauss_Lobatto (degree+1);
            dealii::QGaussLobatto<dim> vol_quad_Gauss_Lobatto (degree+1);
            oned_quad = oned_quad_Gauss_Lobatto;
            volume_quad = vol_quad_Gauss_Lobatto;

            if(dim == 1)
            {
                dealii::QGauss<dim-1> face_quad_Gauss_Legendre (degree+1);
                face_quad = face_quad_Gauss_Legendre;
            }
            else
            {
                dealii::QGaussLobatto<dim-1> face_quad_Gauss_Lobatto (degree+1);
                face_quad = face_quad_Gauss_Lobatto;
            }
        } else {
            dealii::QGauss<1> oned_quad_Gauss_Legendre (degree+1);
            dealii::QGauss<dim> vol_quad_Gauss_Legendre (degree+1);
            dealii::QGauss<dim-1> face_quad_Gauss_Legendre (degree+1);
            oned_quad = oned_quad_Gauss_Legendre;
            volume_quad = vol_quad_Gauss_Legendre;
            face_quad = face_quad_Gauss_Legendre;
        }

        volume_quad_coll.push_back (volume_quad);
        face_quad_coll.push_back (face_quad);
        oned_quad_coll.push_back (oned_quad);

        dealii::FE_DGQArbitraryNodes<dim,dim> lagrange_poly(oned_quad);
        fe_coll_lagr.push_back (lagrange_poly);
    }
    return std::make_tuple(fe_coll, volume_quad_coll, face_quad_coll, oned_quad_coll, fe_coll_lagr);
}


template <int dim, typename real>
void DGBase<dim,real>::set_all_cells_fe_degree ( const unsigned int degree )
{
    triangulation->prepare_coarsening_and_refinement();
    for (auto cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
    {
        if (cell->is_locally_owned()) cell->set_future_fe_index (degree);
    }

    triangulation->execute_coarsening_and_refinement();
}



// Destructor
template <int dim, typename real>
DGBase<dim,real>::~DGBase () 
{ 
    dof_handler.clear ();
}

// template <int dim, typename real>
// void DGBase<dim,real>::set_triangulation(Triangulation *triangulation_input)
// { 
//     // dof_handler.clear();
// 
//     // triangulation = triangulation_input;
// 
//     // dof_handler.initialize(*triangulation, fe_collection);
// 
//     //set_all_cells_fe_degree (fe_collection.size()-1); // Always sets it to the maximum degree
// 
//     //set_all_cells_fe_degree ( max_degree-min_degree);
// }

template <int dim, typename real>
template<typename DoFCellAccessorType>
real DGBase<dim,real>::evaluate_penalty_scaling (
    const DoFCellAccessorType &cell,
    const int iface,
    const dealii::hp::FECollection<dim> fe_collection) const
{

    const unsigned int fe_index = cell->active_fe_index();
    const unsigned int degree = fe_collection[fe_index].tensor_degree();
    const unsigned int degsq = (degree == 0) ? 1 : degree * (degree+1);

    const unsigned int normal_direction = dealii::GeometryInfo<dim>::unit_normal_direction[iface];
    const real vol_div_facearea = cell->extent_in_direction(normal_direction);

    const real penalty = degsq / vol_div_facearea;

    return penalty;
}

template <int dim, typename real>
template<typename DoFCellAccessorType1, typename DoFCellAccessorType2>
bool DGBase<dim,real>::current_cell_should_do_the_work (const DoFCellAccessorType1 &current_cell, const DoFCellAccessorType2 &neighbor_cell) const
{
    if (neighbor_cell->is_ghost()) {
    // In the case the neighbor is a ghost cell, we let the processor with the lower rank do the work on that face
    // We cannot use the cell->index() because the index is relative to the distributed triangulation
    // Therefore, the cell index of a ghost cell might be different to the physical cell index even if they refer to the same cell
        return (current_cell->subdomain_id() < neighbor_cell->subdomain_id());
    } else {
    // Locally owned neighbor cell
        Assert(neighbor_cell->is_locally_owned(), dealii::ExcMessage("If not ghost, neighbor should be locally owned.")); 

        if (current_cell->index() < neighbor_cell->index()) {
        // Cell with lower index does work
            return true;
        } else if (neighbor_cell->index() == current_cell->index()) {
        // If both cells have same index
        // See https://www.dealii.org/developer/doxygen/deal.II/classTriaAccessorBase.html#a695efcbe84fefef3e4c93ee7bdb446ad
        // then cell at the lower level does the work
            return (current_cell->level() < neighbor_cell->level());
        }
        return false;
    }
    Assert(0==1, dealii::ExcMessage("Should not have reached here. Somehow another possible case has not been considered when two cells have the same coarseness."));
    return false;
}

template <int dim, typename real>
template<typename DoFCellAccessorType>
void DGBase<dim,real>::assemble_cell_residual (
    const DoFCellAccessorType &current_cell,
    const bool compute_dRdW,
    dealii::hp::FEValues<dim,dim>        &fe_values_collection_volume,
    dealii::hp::FEFaceValues<dim,dim>    &fe_values_collection_face_int,
    dealii::hp::FEFaceValues<dim,dim>    &fe_values_collection_face_ext,
    dealii::hp::FESubfaceValues<dim,dim> &fe_values_collection_subface,
    dealii::hp::FEValues<dim,dim>        &fe_values_collection_volume_lagrange,
    dealii::LinearAlgebra::distributed::Vector<double> &rhs)
{
    std::vector<dealii::types::global_dof_index> current_dofs_indices;
    std::vector<dealii::types::global_dof_index> neighbor_dofs_indices;

    // Current reference element related to this physical cell
    const int i_fele = current_cell->active_fe_index();
    const int i_quad = i_fele;
    const int i_mapp = 0;

    const dealii::FESystem<dim,dim> &current_fe_ref = fe_collection[i_fele];
    const unsigned int n_dofs_curr_cell = current_fe_ref.n_dofs_per_cell();

    // Local vector contribution from each cell
    dealii::Vector<real> current_cell_rhs (n_dofs_curr_cell); // Defaults to 0.0 initialization

    // Obtain the mapping from local dof indices to global dof indices
    current_dofs_indices.resize(n_dofs_curr_cell);
    current_cell->get_dof_indices (current_dofs_indices);

    fe_values_collection_volume.reinit (current_cell, i_quad, i_mapp, i_fele);
    const dealii::FEValues<dim,dim> &fe_values_volume = fe_values_collection_volume.get_present_fe_values();

    dealii::TriaIterator<dealii::CellAccessor<dim,dim>> cell_iterator = static_cast<dealii::TriaIterator<dealii::CellAccessor<dim,dim>> > (current_cell);
    //if (!(all_parameters->use_weak_form)) fe_values_collection_volume_lagrange.reinit (current_cell, i_quad, i_mapp, i_fele);
    fe_values_collection_volume_lagrange.reinit (cell_iterator, i_quad, i_mapp, i_fele);
    const dealii::FEValues<dim,dim> &fe_values_lagrange = fe_values_collection_volume_lagrange.get_present_fe_values();

    if ( compute_dRdW ) {
        assemble_volume_terms_implicit (fe_values_volume, current_dofs_indices, current_cell_rhs, fe_values_lagrange);
    } else {
        assemble_volume_terms_explicit (fe_values_volume, current_dofs_indices, current_cell_rhs, fe_values_lagrange);
    }

    for (unsigned int iface=0; iface < dealii::GeometryInfo<dim>::faces_per_cell; ++iface) {

        auto current_face = current_cell->face(iface);
        auto neighbor_cell = current_cell->neighbor(iface);

        // Case 1: Face at boundary
        if (current_face->at_boundary() && !current_cell->has_periodic_neighbor(iface) ) {

            fe_values_collection_face_int.reinit(current_cell, iface, i_quad, i_mapp, i_fele);

            // Case 1.1: 1D Periodic boundary condition
            if(current_face->at_boundary() && all_parameters->use_periodic_bc == true && dim == 1) {

                const int neighbor_iface = (iface == 1) ? 0 : 1;

                int cell_index = current_cell->index();
                if (cell_index == 0 && iface == 0) {
                // First cell of the domain, neighbor is the last.
                    neighbor_cell = dof_handler.begin_active();
                    for (unsigned int i = 0 ; i < triangulation->n_active_cells() - 1; ++i) {
                        ++neighbor_cell;
                    }
                } else if (cell_index == (int) triangulation->n_active_cells() - 1 && iface == 1) {
                // Last cell of the domain, neighbor is the first.
                    neighbor_cell = dof_handler.begin_active();
                }

                const int i_fele_n = neighbor_cell->active_fe_index(), i_quad_n = i_fele_n, i_mapp_n = 0;
                const unsigned int n_dofs_neigh_cell = fe_collection[i_fele_n].n_dofs_per_cell();
                neighbor_dofs_indices.resize(n_dofs_neigh_cell);
                fe_values_collection_face_ext.reinit(neighbor_cell, neighbor_iface, i_quad_n, i_mapp_n, i_fele_n);

                const dealii::FEFaceValues<dim,dim> &fe_values_face_int = fe_values_collection_face_int.get_present_fe_values();
                const dealii::FEFaceValues<dim,dim> &fe_values_face_ext = fe_values_collection_face_ext.get_present_fe_values();

                dealii::Vector<real> neighbor_cell_rhs (n_dofs_neigh_cell); // Defaults to 0.0 initialization

                const real penalty1 = evaluate_penalty_scaling (current_cell, iface, fe_collection);
                const real penalty2 = evaluate_penalty_scaling (neighbor_cell, neighbor_iface, fe_collection);
                const real penalty = 0.5 * (penalty1 + penalty2);

                if ( compute_dRdW ) {
                    assemble_face_term_implicit (
                                                fe_values_face_int, fe_values_face_ext,
                                                penalty,
                                                current_dofs_indices, neighbor_dofs_indices,
                                                current_cell_rhs, neighbor_cell_rhs);
                } else {
                    assemble_face_term_explicit (
                                                fe_values_face_int, fe_values_face_ext,
                                                penalty,
                                                current_dofs_indices, neighbor_dofs_indices,
                                                current_cell_rhs, neighbor_cell_rhs);
                }

            } else {
            // Actual boundary term
                const dealii::FEFaceValues<dim,dim> &fe_values_face_int = fe_values_collection_face_int.get_present_fe_values();

                const real penalty = evaluate_penalty_scaling (current_cell, iface, fe_collection);

                const unsigned int boundary_id = current_face->boundary_id();
                // Need to somehow get boundary type from the mesh
                if ( compute_dRdW ) {
                    assemble_boundary_term_implicit (boundary_id, fe_values_face_int, penalty, current_dofs_indices, current_cell_rhs);
                } else {
                    assemble_boundary_term_explicit (boundary_id, fe_values_face_int, penalty, current_dofs_indices, current_cell_rhs);
                }
            }

        //CASE 1.5: periodic boundary conditions
        //note that periodicity is not adapted for hp adaptivity yet. this needs to be figured out in the future
        } else if (current_face->at_boundary() && current_cell->has_periodic_neighbor(iface)){

            neighbor_cell = current_cell->periodic_neighbor(iface);
            //std::cout << "cell " << current_cell->index() << " at boundary" <<std::endl;
            //std::cout << "periodic neighbour on face " << iface << " is " << neighbor_cell->index() << std::endl;


            if (!current_cell->periodic_neighbor_is_coarser(iface) && current_cell_should_do_the_work(current_cell, neighbor_cell)) {

                Assert (current_cell->periodic_neighbor(iface).state() == dealii::IteratorState::valid, dealii::ExcInternalError());

                const unsigned int n_dofs_neigh_cell = fe_collection[neighbor_cell->active_fe_index()].n_dofs_per_cell();
                dealii::Vector<real> neighbor_cell_rhs (n_dofs_neigh_cell); // Defaults to 0.0 initialization

                // Obtain the mapping from local dof indices to global dof indices for neighbor cell
                neighbor_dofs_indices.resize(n_dofs_neigh_cell);
                neighbor_cell->get_dof_indices (neighbor_dofs_indices);

                fe_values_collection_face_int.reinit (current_cell, iface, i_quad, i_mapp, i_fele);
                const dealii::FEFaceValues<dim,dim> &fe_values_face_int = fe_values_collection_face_int.get_present_fe_values();

                // Corresponding face of the neighbor.
                const unsigned int neighbor_iface = current_cell->periodic_neighbor_of_periodic_neighbor(iface);

                const int i_fele_n = neighbor_cell->active_fe_index(), i_quad_n = i_fele_n, i_mapp_n = 0;
                fe_values_collection_face_ext.reinit (neighbor_cell, neighbor_iface, i_quad_n, i_mapp_n, i_fele_n);
                const dealii::FEFaceValues<dim,dim> &fe_values_face_ext = fe_values_collection_face_ext.get_present_fe_values();

                const real penalty1 = evaluate_penalty_scaling (current_cell, iface, fe_collection);
                const real penalty2 = evaluate_penalty_scaling (neighbor_cell, neighbor_iface, fe_collection);
                const real penalty = 0.5 * (penalty1 + penalty2);

                if ( compute_dRdW ) {
                    assemble_face_term_implicit (
                            fe_values_face_int, fe_values_face_ext,
                            penalty,
                            current_dofs_indices, neighbor_dofs_indices,
                            current_cell_rhs, neighbor_cell_rhs);
                } else {
                    assemble_face_term_explicit (
                            fe_values_face_int, fe_values_face_ext,
                            penalty,
                            current_dofs_indices, neighbor_dofs_indices,
                            current_cell_rhs, neighbor_cell_rhs);
                }

                // Add local contribution from neighbor cell to global vector
                for (unsigned int i=0; i<n_dofs_neigh_cell; ++i) {
                    rhs[neighbor_dofs_indices[i]] += neighbor_cell_rhs[i];
                }
            } else {
                //do nothing
            }


        // Case 2:
        // Neighbour is finer occurs if the face has children
        // In this case, we loop over the current large face's subfaces and visit multiple neighbors
        } else if (current_face->has_children()) {

            Assert (current_cell->neighbor(iface).state() == dealii::IteratorState::valid, dealii::ExcInternalError());
            Assert (current_cell->neighbor(iface)->has_children(), dealii::ExcInternalError());

            // Obtain cell neighbour
            const unsigned int neighbor_iface = current_cell->neighbor_face_no(iface);

            for (unsigned int subface_no=0; subface_no < current_face->number_of_children(); ++subface_no) {

                // Get neighbor on ith subface
                auto neighbor_cell = current_cell->neighbor_child_on_subface (iface, subface_no);
                // Since the neighbor cell is finer than the current cell, it should not have more children
                Assert (!neighbor_cell->has_children(), dealii::ExcInternalError());
                Assert (neighbor_cell->neighbor(neighbor_iface) == current_cell, dealii::ExcInternalError());

                const int i_fele_n = neighbor_cell->active_fe_index(), i_quad_n = i_fele_n, i_mapp_n = 0;

                const unsigned int n_dofs_neigh_cell = fe_collection[i_fele_n].n_dofs_per_cell();
                dealii::Vector<real> neighbor_cell_rhs (n_dofs_neigh_cell); // Defaults to 0.0 initialization

                // Obtain the mapping from local dof indices to global dof indices for neighbor cell
                neighbor_dofs_indices.resize(n_dofs_neigh_cell);
                neighbor_cell->get_dof_indices (neighbor_dofs_indices);

                fe_values_collection_subface.reinit (current_cell, iface, subface_no, i_quad, i_mapp, i_fele);
                const dealii::FESubfaceValues<dim,dim> &fe_values_face_int = fe_values_collection_subface.get_present_fe_values();

                fe_values_collection_face_ext.reinit (neighbor_cell, neighbor_iface, i_quad_n, i_mapp_n, i_fele_n);
                const dealii::FEFaceValues<dim,dim> &fe_values_face_ext = fe_values_collection_face_ext.get_present_fe_values();

                const real penalty1 = evaluate_penalty_scaling (current_cell, iface, fe_collection);
                const real penalty2 = evaluate_penalty_scaling (neighbor_cell, neighbor_iface, fe_collection);
                const real penalty = 0.5 * (penalty1 + penalty2);

                if ( compute_dRdW ) {
                    assemble_face_term_implicit (
                            fe_values_face_int, fe_values_face_ext,
                            penalty,
                            current_dofs_indices, neighbor_dofs_indices,
                            current_cell_rhs, neighbor_cell_rhs);
                } else {
                    assemble_face_term_explicit (
                        fe_values_face_int, fe_values_face_ext,
                        penalty,
                        current_dofs_indices, neighbor_dofs_indices,
                        current_cell_rhs, neighbor_cell_rhs);
                }
                // Add local contribution from neighbor cell to global vector
                for (unsigned int i=0; i<n_dofs_neigh_cell; ++i) {
                    rhs[neighbor_dofs_indices[i]] += neighbor_cell_rhs[i];
                }
            }

        // Case 3:
        // Neighbor cell is NOT coarser
        // Therefore, they have the same coarseness, and we need to choose one of them to do the work
        } else if ( !(current_cell->neighbor_is_coarser(iface)) && current_cell_should_do_the_work(current_cell, neighbor_cell) ) {
            Assert (current_cell->neighbor(iface).state() == dealii::IteratorState::valid, dealii::ExcInternalError());

            auto neighbor_cell = current_cell->neighbor_or_periodic_neighbor(iface);
            // Corresponding face of the neighbor.
            // e.g. The 4th face of the current cell might correspond to the 3rd face of the neighbor
            const unsigned int neighbor_iface = current_cell->neighbor_of_neighbor(iface);

            // Get information about neighbor cell
            const unsigned int n_dofs_neigh_cell = fe_collection[neighbor_cell->active_fe_index()].n_dofs_per_cell();

            // Local rhs contribution from neighbor
            dealii::Vector<real> neighbor_cell_rhs (n_dofs_neigh_cell); // Defaults to 0.0 initialization

            // Obtain the mapping from local dof indices to global dof indices for neighbor cell
            neighbor_dofs_indices.resize(n_dofs_neigh_cell);
            neighbor_cell->get_dof_indices (neighbor_dofs_indices);

            fe_values_collection_face_int.reinit (current_cell, iface, i_quad, i_mapp, i_fele);
            const dealii::FEFaceValues<dim,dim> &fe_values_face_int = fe_values_collection_face_int.get_present_fe_values();

            const int i_fele_n = neighbor_cell->active_fe_index(), i_quad_n = i_fele_n, i_mapp_n = 0;
            fe_values_collection_face_ext.reinit (neighbor_cell, neighbor_iface, i_quad_n, i_mapp_n, i_fele_n);
            const dealii::FEFaceValues<dim,dim> &fe_values_face_ext = fe_values_collection_face_ext.get_present_fe_values();

            const real penalty1 = evaluate_penalty_scaling (current_cell, iface, fe_collection);
            const real penalty2 = evaluate_penalty_scaling (neighbor_cell, neighbor_iface, fe_collection);
            const real penalty = 0.5 * (penalty1 + penalty2);

            if ( compute_dRdW ) {
                assemble_face_term_implicit (
                        fe_values_face_int, fe_values_face_ext,
                        penalty,
                        current_dofs_indices, neighbor_dofs_indices,
                        current_cell_rhs, neighbor_cell_rhs);
            } else {
                assemble_face_term_explicit (
                        fe_values_face_int, fe_values_face_ext,
                        penalty,
                        current_dofs_indices, neighbor_dofs_indices,
                        current_cell_rhs, neighbor_cell_rhs);
            }

            // Add local contribution from neighbor cell to global vector
            for (unsigned int i=0; i<n_dofs_neigh_cell; ++i) {
                rhs[neighbor_dofs_indices[i]] += neighbor_cell_rhs[i];
            }
        } else {
            // Case 4: Neighbor is coarser
            // Do nothing.
            // The face contribution from the current cell will appear then the coarse neighbor checks for subfaces
        }

    } // end of face loop

    // Add local contribution from current cell to global vector
    for (unsigned int i=0; i<n_dofs_curr_cell; ++i) {
        rhs[current_dofs_indices[i]] += current_cell_rhs[i];
    }
}


template <int dim, typename real>
void DGBase<dim,real>::assemble_residual (const bool compute_dRdW)
{
    right_hand_side = 0;

    if (compute_dRdW) system_matrix = 0;

    //dealii::hp::MappingCollection<dim> mapping_collection(*(high_order_grid.mapping_fe_field));
    //const dealii::MappingManifold<dim,dim> mapping;
    //const dealii::MappingQ<dim,dim> mapping(max_degree+1);
    const auto mapping = (*(high_order_grid.mapping_fe_field));
    dealii::hp::MappingCollection<dim> mapping_collection(mapping);

    dealii::hp::FEValues<dim,dim>        fe_values_collection_volume (mapping_collection, fe_collection, volume_quadrature_collection, this->volume_update_flags); ///< FEValues of volume.
    dealii::hp::FEFaceValues<dim,dim>    fe_values_collection_face_int (mapping_collection, fe_collection, face_quadrature_collection, this->face_update_flags); ///< FEValues of interior face.
    dealii::hp::FEFaceValues<dim,dim>    fe_values_collection_face_ext (mapping_collection, fe_collection, face_quadrature_collection, this->neighbor_face_update_flags); ///< FEValues of exterior face.
    dealii::hp::FESubfaceValues<dim,dim> fe_values_collection_subface (mapping_collection, fe_collection, face_quadrature_collection, this->face_update_flags); ///< FEValues of subface.

    dealii::hp::FEValues<dim,dim>        fe_values_collection_volume_lagrange (mapping_collection, fe_collection_lagrange, volume_quadrature_collection, this->volume_update_flags);

    solution.update_ghost_values();

    for (auto current_cell = dof_handler.begin_active(); current_cell != dof_handler.end(); ++current_cell) {
        if (!current_cell->is_locally_owned()) continue;

        // Add right-hand side contributions this cell can compute
        assemble_cell_residual (current_cell, compute_dRdW,
            fe_values_collection_volume,
            fe_values_collection_face_int,
            fe_values_collection_face_ext,
            fe_values_collection_subface,
            fe_values_collection_volume_lagrange,
            right_hand_side);
    } // end of cell loop

    right_hand_side.compress(dealii::VectorOperation::add);
    if ( compute_dRdW ) system_matrix.compress(dealii::VectorOperation::add);

} // end of assemble_system_explicit ()


template <int dim, typename real>
double DGBase<dim,real>::get_residual_l2norm () const
{
    return right_hand_side.l2_norm();
}
template <int dim, typename real>
unsigned int DGBase<dim,real>::n_dofs () const
{
    return dof_handler.n_dofs();
}

template <int dim, typename real>
void DGBase<dim,real>::output_results_vtk (const unsigned int cycle)// const
{

    dealii::DataOut<dim, dealii::hp::DoFHandler<dim>> data_out;
    data_out.attach_dof_handler (dof_handler);

    //std::vector<std::string> solution_names;
    //for(int s=0;s<nstate;++s) {
    //    std::string varname = "u" + dealii::Utilities::int_to_string(s,1);
    //    solution_names.push_back(varname);
    //}
    //std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(nstate, dealii::DataComponentInterpretation::component_is_scalar);
    //data_out.add_data_vector (solution, solution_names, dealii::DataOut<dim>::type_dof_data, data_component_interpretation);

    const std::unique_ptr< dealii::DataPostprocessor<dim> > post_processor = Postprocess::PostprocessorFactory<dim>::create_Postprocessor(all_parameters);
    data_out.add_data_vector (solution, *post_processor);

    dealii::Vector<float> subdomain(triangulation->n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i) {
        subdomain(i) = triangulation->locally_owned_subdomain();
    }
    data_out.add_data_vector(subdomain, "subdomain", dealii::DataOut_DoFData<dealii::hp::DoFHandler<dim>,dim>::DataVectorType::type_cell_data);

    // Output the polynomial degree in each cell
    std::vector<unsigned int> active_fe_indices;
    dof_handler.get_active_fe_indices(active_fe_indices);
    dealii::Vector<double> active_fe_indices_dealiivector(active_fe_indices.begin(), active_fe_indices.end());
    dealii::Vector<double> cell_poly_degree = active_fe_indices_dealiivector;

//    int index = 0;
//    for (auto current_cell_poly = cell_poly_degree.begin(); current_cell_poly != cell_poly_degree.end(); ++current_cell_poly) {
//        current_cell_poly[index] = fe_collection[active_fe_indices_dealiivector[index]].tensor_degree();
//        index++;
//    }
//    //using DVTenum = dealii::DataOut_DoFData<dealii::hp::DoFHandler<dim>,dim>::DataVectorType;
//    data_out.add_data_vector (cell_poly_degree, "PolynomialDegree", dealii::DataOut_DoFData<dealii::hp::DoFHandler<dim>,dim>::DataVectorType::type_cell_data);
    data_out.add_data_vector (active_fe_indices_dealiivector, "PolynomialDegree", dealii::DataOut_DoFData<dealii::hp::DoFHandler<dim>,dim>::DataVectorType::type_cell_data);


    //assemble_residual (false);
    std::vector<std::string> residual_names;
    for(int s=0;s<nstate;++s) {
        std::string varname = "residual" + dealii::Utilities::int_to_string(s,1);
        residual_names.push_back(varname);
    }
    //std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(nstate, dealii::DataComponentInterpretation::component_is_scalar);
    //data_out.add_data_vector (right_hand_side, residual_names, dealii::DataOut<dim, dealii::hp::DoFHandler<dim>>::type_dof_data, data_component_interpretation);
    data_out.add_data_vector (right_hand_side, residual_names, dealii::DataOut_DoFData<dealii::hp::DoFHandler<dim>,dim>::DataVectorType::type_dof_data);


    const int iproc = dealii::Utilities::MPI::this_mpi_process(mpi_communicator);
    //data_out.build_patches (mapping_collection[mapping_collection.size()-1]);
    data_out.build_patches(*(high_order_grid.mapping_fe_field), max_degree, dealii::DataOut<dim, dealii::hp::DoFHandler<dim>>::CurvedCellRegion::no_curved_cells);
    //data_out.build_patches(*(high_order_grid.mapping_fe_field), fe_collection.size(), dealii::DataOut<dim>::CurvedCellRegion::curved_inner_cells);
    std::string filename = "solution-" + dealii::Utilities::int_to_string(dim, 1) +"D_maxpoly"+dealii::Utilities::int_to_string(max_degree, 2)+"-";
    filename += dealii::Utilities::int_to_string(cycle, 4) + ".";
    filename += dealii::Utilities::int_to_string(iproc, 4);
    filename += ".vtu";
    std::ofstream output(filename);
    data_out.write_vtu(output);

    if (iproc == 0) {
        std::vector<std::string> filenames;
        for (unsigned int iproc = 0; iproc < dealii::Utilities::MPI::n_mpi_processes(mpi_communicator); ++iproc) {
            std::string fn = "solution-" + dealii::Utilities::int_to_string(dim, 1) +"D_maxpoly"+dealii::Utilities::int_to_string(max_degree, 2)+"-";
            fn += dealii::Utilities::int_to_string(cycle, 4) + ".";
            fn += dealii::Utilities::int_to_string(iproc, 4);
            fn += ".vtu";
            filenames.push_back(fn);
        }
        std::string master_fn = "solution-" + dealii::Utilities::int_to_string(dim, 1) +"D_maxpoly"+dealii::Utilities::int_to_string(max_degree, 2)+"-";
        master_fn += dealii::Utilities::int_to_string(cycle, 4) + ".pvtu";
        std::ofstream master_output(master_fn);
        data_out.write_pvtu_record(master_output, filenames);
    }

}

template <int dim, typename real>
void DGBase<dim,real>::allocate_system ()
{
    pcout << "Allocating DG system and initializing FEValues" << std::endl;
    // This function allocates all the necessary memory to the 
    // system matrices and vectors.

    dof_handler.distribute_dofs(fe_collection);

    //dealii::MappingFEField<dim,dim,dealii::LinearAlgebra::distributed::Vector<double>, dealii::DoFHandler<dim>> mapping = high_order_grid.get_MappingFEField();
    //dealii::MappingFEField<dim,dim,dealii::LinearAlgebra::distributed::Vector<double>, dealii::DoFHandler<dim>> mapping = *(high_order_grid.mapping_fe_field);

    //int minimum_degree = (all_parameters->use_collocated_nodes==true) ?  1 :  0;
    //int current_fe_index = 0;
    //for (unsigned int degree=minimum_degree; degree<=max_degree; ++degree) {
    //    //mapping_collection.push_back(mapping);
    //    if(current_fe_index <= mapping_collection.size()) {
    //        mapping_collection.push_back(mapping);
    //    } else {
    //        mapping_collection[current_fe_index] = std::shared_ptr<const dealii::Mapping<dim, dim>>(mapping.clone());
    //    }
    //    current_fe_index++;
    //}

    // Solution and RHS
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    dealii::DoFTools::extract_locally_relevant_dofs(dof_handler, ghost_dofs);
    locally_relevant_dofs = ghost_dofs;
    ghost_dofs.subtract_set(locally_owned_dofs);
    //dealii::DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    solution.reinit(locally_owned_dofs, ghost_dofs, mpi_communicator);
    //right_hand_side.reinit(locally_owned_dofs, mpi_communicator);
    right_hand_side.reinit(locally_owned_dofs, ghost_dofs, mpi_communicator);

    // System matrix allocation
    dealii::DynamicSparsityPattern dsp(locally_relevant_dofs);
    dealii::DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
    dealii::SparsityTools::distribute_sparsity_pattern(dsp, dof_handler.compute_n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);

    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(locally_owned_dofs, sparsity_pattern, mpi_communicator);
}

template <int dim, typename real>
void DGBase<dim,real>::evaluate_mass_matrices (bool do_inverse_mass_matrix)
{
    // Mass matrix sparsity pattern
    //dealii::SparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs(), dof_handler.get_fe_collection().max_dofs_per_cell());
    //dealii::SparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs(), dof_handler.get_fe_collection().max_dofs_per_cell());
    //dealii::DynamicSparsityPattern dsp(dof_handler.n_locally_owned_dofs(), dof_handler.n_locally_owned_dofs(), dof_handler.get_fe_collection().max_dofs_per_cell());
    dealii::DynamicSparsityPattern dsp(dof_handler.n_dofs());
    std::vector<dealii::types::global_dof_index> dofs_indices;
    for (auto cell = dof_handler.begin_active(); cell!=dof_handler.end(); ++cell) {

        if (!cell->is_locally_owned()) continue;

        const unsigned int fe_index_curr_cell = cell->active_fe_index();

        // Current reference element related to this physical cell
        const dealii::FESystem<dim,dim> &current_fe_ref = fe_collection[fe_index_curr_cell];
        const unsigned int n_dofs_cell = current_fe_ref.n_dofs_per_cell();

        dofs_indices.resize(n_dofs_cell);
        cell->get_dof_indices (dofs_indices);
        for (unsigned int itest=0; itest<n_dofs_cell; ++itest) {
            for (unsigned int itrial=0; itrial<n_dofs_cell; ++itrial) {
                dsp.add(dofs_indices[itest], dofs_indices[itrial]);
            }
        }
    }
    dealii::SparsityTools::distribute_sparsity_pattern(dsp, dof_handler.compute_n_locally_owned_dofs_per_processor(), mpi_communicator, locally_owned_dofs);
    mass_sparsity_pattern.copy_from(dsp);
    if (do_inverse_mass_matrix == true) {
        global_inverse_mass_matrix.reinit(locally_owned_dofs, mass_sparsity_pattern);
    } else {
        global_mass_matrix.reinit(locally_owned_dofs, mass_sparsity_pattern);
    }

    //dealii::TrilinosWrappers::SparseMatrix 
    //    matrix_with_correct_size(locally_owned_dofs,
    //            mpi_communicator,
    //            dof_handler.get_fe_collection().max_dofs_per_cell());
    //pcout << "Before compress" << std::endl;
    //matrix_with_correct_size.compress(dealii::VectorOperation::unknown);
    //if (do_inverse_mass_matrix == true) {
    //    global_inverse_mass_matrix.reinit(matrix_with_correct_size);
    //} else {
    //    global_mass_matrix.reinit(matrix_with_correct_size);
    //}
    //pcout << "AFter reinit" << std::endl;

    //dealii::hp::MappingCollection<dim> mapping_collection(*(high_order_grid.mapping_fe_field));
    //const dealii::MappingManifold<dim,dim> mapping;
    //const dealii::MappingQ<dim,dim> mapping(max_degree+1);
    const auto mapping = (*(high_order_grid.mapping_fe_field));
    dealii::hp::MappingCollection<dim> mapping_collection(mapping);

    dealii::hp::FEValues<dim,dim> fe_values_collection_volume (mapping_collection, fe_collection, volume_quadrature_collection, this->volume_update_flags); ///< FEValues of volume.
    for (auto cell = dof_handler.begin_active(); cell!=dof_handler.end(); ++cell) {

        if (!cell->is_locally_owned()) continue;

        const unsigned int mapping_index = 0;
        const unsigned int fe_index_curr_cell = cell->active_fe_index();
        const unsigned int quad_index = fe_index_curr_cell;

        // Current reference element related to this physical cell
        const dealii::FESystem<dim,dim> &current_fe_ref = fe_collection[fe_index_curr_cell];
        const unsigned int n_dofs_cell = current_fe_ref.n_dofs_per_cell();
        const unsigned int n_quad_pts = volume_quadrature_collection[fe_index_curr_cell].size();

        dealii::FullMatrix<real> local_mass_matrix(n_dofs_cell);

        fe_values_collection_volume.reinit (cell, quad_index, mapping_index, fe_index_curr_cell);
        const dealii::FEValues<dim,dim> &fe_values_volume = fe_values_collection_volume.get_present_fe_values();

        for (unsigned int itest=0; itest<n_dofs_cell; ++itest) {
            const unsigned int istate_test = fe_values_volume.get_fe().system_to_component_index(itest).first;
            for (unsigned int itrial=itest; itrial<n_dofs_cell; ++itrial) {
                const unsigned int istate_trial = fe_values_volume.get_fe().system_to_component_index(itrial).first;
                real value = 0.0;
                for (unsigned int iquad=0; iquad<n_quad_pts; ++iquad) {
                    value +=
                        fe_values_volume.shape_value_component(itest,iquad,istate_test)
                        * fe_values_volume.shape_value_component(itrial,iquad,istate_trial)
                        * fe_values_volume.JxW(iquad);
                }
                local_mass_matrix[itrial][itest] = 0.0;
                local_mass_matrix[itest][itrial] = 0.0;
                if(istate_test==istate_trial) { 
                    local_mass_matrix[itrial][itest] = value;
                    local_mass_matrix[itest][itrial] = value;
                }
            }
        }

        dofs_indices.resize(n_dofs_cell);
        cell->get_dof_indices (dofs_indices);
        if (do_inverse_mass_matrix == true) {
            dealii::FullMatrix<real> local_inverse_mass_matrix(n_dofs_cell);
            local_inverse_mass_matrix.invert(local_mass_matrix);
            global_inverse_mass_matrix.set (dofs_indices, local_inverse_mass_matrix);
        } else {
            global_mass_matrix.set (dofs_indices, local_mass_matrix);
        }
    }

    if (do_inverse_mass_matrix == true) {
        global_inverse_mass_matrix.compress(dealii::VectorOperation::insert);
    } else {
        global_mass_matrix.compress(dealii::VectorOperation::insert);
    }

    return;
}
template<int dim, typename real>
void DGBase<dim,real>::add_mass_matrices(const real scale)
{
    system_matrix.add(scale, global_mass_matrix);
}

template <int dim, typename real>
std::vector<real> DGBase<dim,real>::evaluate_time_steps (const bool exact_time_stepping)
{
    // TO BE DONE
    std::vector<real> time_steps(10);
    if(exact_time_stepping) return time_steps;
    return time_steps;
}

template class DGBase <PHILIP_DIM, double>;
template class DGFactory <PHILIP_DIM, double>;

} // PHiLiP namespace

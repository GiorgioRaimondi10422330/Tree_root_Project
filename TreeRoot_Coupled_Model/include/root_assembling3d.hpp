/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   assembling3d.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   January 2016.
  @brief  Miscelleanous assembly routines for the 3D tissue problem.
 */
#ifndef M3D1D_ASSEMBLING_ROOT_3D_HPP_
#define M3D1D_ASSEMBLING_ROOT_3D_HPP_

#include <defines.hpp>
#include <utilities.hpp>

namespace getfem {

//! Build the mass and divergence matrices for the 3D Darcy's problem,
//! @f$ M = \int_{\Omega} \frac{1}{\kappa}\,\mathbf{u}\cdot\mathbf{v}~dx @f$ and
//! @f$ D = \int_{\Omega} div(\mathbf{u})\,p~dx @f$
/*!
	@param M     Darcy's mass matrix
	@param D     Darcy's divergence matrix
	@param mim   The integration method to be used
	@param mf_u  The finite element method for the velocity @f$ \mathbf{u} @f$
	@param mf_p  The finite element method for the pressure @f$ p @f$
	@param rg    The region where to integrate

	@ingroup asm
 */ 
template<typename MAT>
void 
asm_tissue_darcy_divergence
	(MAT & D,
	 const mesh_im & mim,
	 const mesh_fem & mf_u,
	 const mesh_fem & mf_p,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf_p.get_qdim() == 1, 
		"invalid data mesh fem for pressure (Qdim=1 required)");
	GMM_ASSERT1(mf_u.get_qdim() > 1, 
		"invalid data mesh fem for velocity (Qdim>1 required)");
	// Build the divergence matrix Dtt
	generic_assembly 
	assemP("M$1(#2,#1)+=comp(Base(#2).vGrad(#1))(:,:,i,i);");
	assemP.push_mi(mim);
	assemP.push_mf(mf_u);
	assemP.push_mf(mf_p);
	assemP.push_mat(D);
	assemP.assembly(rg);

}


template<typename VEC>
void 
asm_tissue_darcy_gravity
	(VEC & F,
	 const mesh_im & mim,
	 const mesh_fem & mf_u,
	 const scalar_type & C,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf_u.get_qdim() > 1, 
		"invalid data mesh fem for velocity (Qdim>1 required)");
	//Build the rhs gravity contribution
	generic_assembly 
	assemH( "V$1(#1)+=-comp(vBase(#1))(:,3);");
	assemH.push_mi(mim);
	assemH.push_mf(mf_u);
	assemH.push_vec(F);
	assemH.assembly(rg);

	gmm::scaled(F,C);
}




//! Aux function to extract the radius of the ith branch, R[i] 
template
<typename VEC>
scalar_type
simple_compute_radius
	(const mesh_im & mim,
	 const mesh_fem & mf_coef,
	 const VEC & R,
	 const size_type & rg
	 ) 
{
	vector_type dof_enum;
	size_type fine=0;
	for(getfem::mr_visitor mrv(mf_coef.linked_mesh().region(rg)); !mrv.finished(); ++mrv){
				for(auto b: mf_coef.ind_basic_dof_of_element(mrv.cv())){
					dof_enum.emplace_back(b);
					fine++;
				}
			}
	size_type first_=dof_enum[0];
	return R[first_];
}

}//End Namespace

#endif 

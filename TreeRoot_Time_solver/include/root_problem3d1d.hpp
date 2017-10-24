/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2017 Giorgio Raimondi
======================================================================*/
/*! 
  @file   c_problem3d1d.hpp
  @author Raimondi Giorgio <giorgio3.raimondi@mail.polimi.it>
  @date   May 2016.
  @brief  Declaration of the main class for the 3D/1D coupled curved problem.
 */

#ifndef C_M3D1D_PROBLEM3D1D_HPP_
#define C_M3D1D_PROBLEM3D1D_HPP_

//base class
#include <problem3d1d.hpp>
//root class headers
#include <root_descr3d1d.hpp>  
#include <root_descr3d1d.hpp>
#include <root_mesh1d.hpp>
#include <root_param3d1d.hpp>
#include <root_assembling3d.hpp>

namespace getfem {

//!	Main class defining the coupled 3D/1D fluid curved problem.
class root_problem3d1d: problem3d1d {

public:

	//! Initialize the problem
	/*!
		1. Read the .param filename from standard input
		2. Import problem descriptors (file paths, GetFEM types, ...)
		3. Import mesh for tissue (3D) and vessel network (1D)
		4. Set finite elements and integration methods
		5. Build problem parameters
		6. Build the list of tissue boundary data
		7. Build the list of vessel boundary (and junction) data
	 */
	void init (int argc, char *argv[]);
	//! Assemble the problem
	/*!
		1. Initialize problem matrices and vectors
		2. Build the monolithic matrix AM
		3. Build the monolithic rhs FM
	 */
	void assembly (void);
	//! Solve the problem
	/*!
		Solve the monolithic system AM*UM=FM (direct or iterative)
	 */
	bool solve (void);
	//! Export results into vtk files
	/*!
		Export solutions Ut, Pt, Uv, Pv from the monolithic array UM
	 */
	void export_vtk (const string & suff = "");

	//! Compute mean tissue pressure
	inline scalar_type mean_pt (void){ 
		return asm_mean(mf_Pt, mimt, 
			gmm::sub_vector(UM, gmm::sub_interval(dof.Ut(), dof.Pt()))); 
	}
	//! Compute mean vessel pressure
	inline scalar_type mean_pv (void){ 
		return asm_mean(mf_Pv, mimv, 
			gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()))); 
	}
 
	//! Compute mean vessel velocity
	inline scalar_type mean_uv (void){
		size_type dofi=0;
		size_type shift=dof.Ut()+dof.Pt();
		scalar_type mean=0.0;
		for(size_type branch=0; branch<mf_Uvi.size(); branch++){
			shift+=dofi;
			dofi=mf_Uvi[branch].nb_dof();
			mean+=asm_mean(mf_Uvi[branch],mimv,gmm::sub_vector(UM, gmm::sub_interval(shift, dofi)));
		}
		return mean/mf_Uvi.size();
	}

	//! Compute the dynamic pressure gain

	//! Compute total flow rate (network to tissue) - pressure
	inline scalar_type flow_rate (void) { return TFR; };

	void Uplot(void){ 
		size_type shift=dof.Ut()+dof.Pt();
		size_type Cb=0;
		cout<<"Uv root_problem3d1d:\n\n";
		for(size_type b=0;b<nb_branches;++b){
			cout<<" b"<<b<<"= [";
			if(b>0) shift+= Cb;
			Cb=mf_Uvi[b].nb_dof();
			for(size_type i=0; i<Cb; ++i)
				cout<<" "<<UM[shift+i];
			cout<<" ]\n"; 
		}
		cout<<"\n";
		shift+=Cb;
		cout<<" Pressure=[";
		for(size_type ss=shift; ss<shift+dof.Pv(); ss++)
			cout<<" "<<UM[ss];
		cout<<" ]\n\n";
	}

protected:
	//!	Algorithm description strings for curved model
	root_descr3d1d r_descr;
	//! Physical parameters (dimensionless) for curved model
	root_param3d1d r_param;
	//! Pressure in the poruos medium at a previous iteration
	vector_type Pt_old_ITER;
	//! Pressure in the porous medium at a previous iteration on data mesh_fem
	vector_type Pt_data;
	//! Non Linear matrix containing Darcy term
	sparse_matrix_type NLMtt;
	//! Monolitic matrix with both Linear and non Linear term
	sparse_matrix_type CM;
	//! Monolitic RHS with both Linear and non Linear term
	vector_type CFM;


	//! Timestep Variable
	//! Pressure of the porous medium at the previous time step
	vector_type Pt_Prev;
	//! Non Linear matrix containing Richards term
	sparse_matrix_type NLRtt;
	//! Non Linear Vector containing Richards term
	vector_type NLFtt;
	//! Time step
	scalar_type dT;
	//! Maximum Time 
	scalar_type maxT;
	//! Porosity
	vector_type Porosity;
	//! Pressure Head
	vector_type H;
	//! Conductivity
	vector_type Conductivity;

	// Aux methods for init
	//! Import algorithm specifications
	void import_data(void);
	//! Import mesh for tissue (3D) and vessel (1D)  
	void build_mesh(void); 
	//! Set finite elements methods and integration methods 
	void set_im_and_fem(void);
	//! Build problem parameters
	void build_param(void);
	//! Build the list of tissue boundary data 
	/*!	Face numbering:
		  0 : {x = 0 }  "back"
		  1 : {x = Lx}  "front"
		  2 : {y = 0 }  "left"
		  3 : {y = Ly}  "right"
		  4 : {z = 0 }  "bottom"
		  5 : {z = Lz}  "top"
	*/
	void build_tissue_boundary(void);
	//! Build the list of vessel boundary (and junctions) data 
	void build_vessel_boundary(void);
	// Aux methods for assembly
	//! Build the monolithic matrix AM by blocks
	void assembly_mat(bool Q0=false);
	//! Build the monolithic rhs FM by blocks
	void assembly_rhs(void);
	//! Assemble RHS source term for stand-alone tissue problem
	void assembly_tissue_test_rhs(void);
	//Auxiliary method for iterative solver
	//! Build non linear Matrix NL
	void assembly_nonlinear_mat(size_type iter=0, size_type tempo=0);
	//! Solve the monolitic sistem for a given iteration of the fixed point method
	bool solve_iter(size_type tempo=0);
	//! Solve prepare the monolitic system in order to solve a certain time step
	bool solve_time(size_type tempo=0);
	//! Set the initial condition
	bool set_initial_condition(void);

}; /* end of class problem3d1d */


} /* end of namespace */

#endif
 

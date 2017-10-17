/* -*- c++ -*- (enableMbars emacs c++ mode) */ 
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Curved Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2017 Giorgio Raimondi
======================================================================*/
/*! 
  @file   c_problem3d1d.cpp
  @author Giorgio Raimondi <giorgio3.raimondi@mail.polimi.it>
  @date   May 2017.
  @brief  Definition of the main class for the 3D/1D coupled curved problem.
 */    
#include <root_problem3d1d.hpp>
 
namespace getfem {


/////////// Initialize the problem ///////////////////////////////////// 
void 
root_problem3d1d::init(int argc, char *argv[])
{
	//1. Read the .param filename from standard input
	PARAM.read_command_line(argc, argv);
	//2. Import data (algorithm specifications, boundary conditions, ...)
	import_data();
	//3. Import mesh for tissue (3D) and vessel network (1D)
	build_mesh();
	//4. Set finite elements and integration methods
	set_im_and_fem();
	//5. Build problem parameters
	build_param();
	//6. Build the list of tissue boundary data
	build_tissue_boundary();
	//7. Build the list of tissue boundary (and junction) data
	build_vessel_boundary();
}

void
root_problem3d1d::import_data(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Importing descriptors for tissue and vessel problems ..." << endl;
	#endif
	descr.import(PARAM);
	if(PARAM.int_value("CURVE_PROBLEM"))
		r_descr.import(PARAM);
	#ifdef M3D1D_VERBOSE_
	cout << descr;
	if(PARAM.int_value("CURVE_PROBLEM"))
		cout << r_descr;
	#endif
}

void 
root_problem3d1d::build_mesh(void)
{
	bool test = 0;
	test = PARAM.int_value("TEST_GEOMETRY");
	if(test==0){
		#ifdef M3D1D_VERBOSE_
		cout << "Importing the 3D mesh for the tissue ...  "   << endl;
		#endif
		 import_msh_file(descr.MESH_FILET, mesht);
	}else{
		#ifdef M3D1D_VERBOSE_
		cout << "Building the regular 3D mesh for the tissue ...  "   << endl;
		#endif
		string st("GT='" + PARAM.string_value("GT_T") + "'; " +
					   "NSUBDIV=" + PARAM.string_value("NSUBDIV_T") + "; " +  
					   "ORG=" + PARAM.string_value("ORG_T") + "; " +  
					   "SIZES=" + PARAM.string_value("SIZES_T") + "; " +  
					   "NOISED=" + PARAM.string_value("NOISED_T")); 
		#ifdef M3D1D_VERBOSE_
		cout << "mesht description: " << st << endl;
		#endif
		regular_mesh(mesht, st); 
	}
 
	#ifdef M3D1D_VERBOSE_
	cout << "Importing the 1D mesh for the vessel ... "   << endl;
	#endif
	std::ifstream ifs(descr.MESH_FILEV);
	GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr.MESH_FILEV);
	bool Import=PARAM.int_value("IMPORT_CURVE");
	bool Curve=PARAM.int_value("CURVE_PROBLEM");

	if(Curve && !Import){
		import_pts_file(ifs, meshv, BCv, nb_vertices, descr.MESH_TYPEV,r_param);
	}
	else if(Import && !Curve){
		GMM_ASSERT1(0,"If you want to import the curvature, you need to enable CURVE_PROBLEM=1");
	}
	else if(Import && Curve){
		std::ifstream ifc(PARAM.string_value("CURVE_FILE","curvature file location"));
		GMM_ASSERT1(ifc.good(), "impossible to read from file " << PARAM.string_value("CURVE_FILE","curvature file location"));
		
		import_pts_file(ifs,ifc, meshv, BCv, nb_vertices, descr.MESH_TYPEV, r_param);

		ifc.close();
	} else{
		import_pts_file(ifs, meshv, BCv, nb_vertices, descr.MESH_TYPEV);
	}

	nb_branches = nb_vertices.size();
	ifs.close();
}

void
root_problem3d1d::set_im_and_fem(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs for tissue and vessel problems ..." << endl;
	#endif
	pintegration_method pim_t = int_method_descriptor(descr.IM_TYPET);
	pintegration_method pim_v = int_method_descriptor(descr.IM_TYPEV);
	mimt.set_integration_method(mesht.convex_index(), pim_t);
	mimv.set_integration_method(meshv.convex_index(), pim_v);
	
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEMs for tissue and vessel problems ..." << endl;
	#endif
	bgeot::pgeometric_trans pgt_t = bgeot::geometric_trans_descriptor(descr.MESH_TYPET);
	bgeot::pgeometric_trans pgt_v = bgeot::geometric_trans_descriptor(descr.MESH_TYPEV);
	pfem pf_Ut = fem_descriptor(descr.FEM_TYPET);
	pfem pf_Pt = fem_descriptor(descr.FEM_TYPET_P);
	pfem pf_Uv = fem_descriptor(descr.FEM_TYPEV); 
	pfem pf_Pv = fem_descriptor(descr.FEM_TYPEV_P);
	pfem pf_coeft = fem_descriptor(descr.FEM_TYPET_DATA);
	pfem pf_coefv = fem_descriptor(descr.FEM_TYPEV_DATA);
	DIMT = pgt_t->dim();	//DIMV = 1;  
	mf_Ut.set_qdim(bgeot::dim_type(DIMT)); 
	 
	mf_Ut.set_finite_element(mesht.convex_index(), pf_Ut);
	GMM_ASSERT1(mf_Ut.get_qdim() == mf_Ut.fem_of_element(0)->target_dim(), 
		"Intrinsic vectorial FEM used"); // RT0 IS INTRINSIC VECTORIAL!!!
	mf_Pt.set_finite_element(mesht.convex_index(),  pf_Pt);
	mf_coeft.set_finite_element(mesht.convex_index(), pf_coeft); 
	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for vessel branches ..." << endl;
	#endif
	mf_Uvi.reserve(nb_branches);
	mf_coefvi.reserve(nb_branches);
	for(size_type i=0; i<nb_branches; ++i){ 	
		mesh_fem mf_tmp(meshv);
		mf_tmp.set_finite_element(meshv.region(i).index(), pf_coefv);
		mf_coefvi.emplace_back(mf_tmp);
		mf_tmp.clear();
		
		mf_tmp.set_finite_element(meshv.region(i).index(), pf_Uv);
		mf_Uvi.emplace_back(mf_tmp);
		mf_tmp.clear();
	}
	mf_Pv.set_finite_element(meshv.convex_index(), pf_Pv);
	mf_coefv.set_finite_element(meshv.convex_index(), pf_coefv);
	
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEM dimensions for tissue and vessel problems ..." << endl;
	#endif
	dof.set(mf_Ut, mf_Pt, mf_Uvi, mf_Pv, mf_coeft, mf_coefv);
	#ifdef M3D1D_VERBOSE_
	cout << std::scientific << dof;
	#endif
}

void
root_problem3d1d::build_param(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Building parameters for tissue and vessel problems ..." << endl;
	#endif
	r_param.build(PARAM, mf_coeft, mf_coefv,mf_coefvi);
	#ifdef M3D1D_VERBOSE_
	cout << r_param ;
	#endif
}

void 
root_problem3d1d::build_tissue_boundary (void) 
{
 	#ifdef M3D1D_VERBOSE_
	cout << "Building tissue boundary ..." << endl;
	#endif
	BCt.clear();
	BCt.reserve(2*DIMT);
	// Parse BC data
	string label_in = PARAM.string_value("BClabel", "Array of tissue boundary labels");
	string value_in = PARAM.string_value("BCvalue", "Array of tissue boundary values");
	vector<string> labels = split(label_in, ' ');
	vector<string> values = split(value_in, ' ');
	GMM_ASSERT1(labels.size()==2*DIMT, "wrong number of BC labels");
	GMM_ASSERT1(values.size()==2*DIMT, "wrong number of BC values");
	for (unsigned f=0; f<2*DIMT; ++f) {
		BCt.emplace_back(labels[f], std::stof(values[f]), 0, f);
		#ifdef M3D1D_VERBOSE_
		cout << "  face " << f << " : " << BCt.back() << endl;
		#endif
	}
	// Build mesht regions
	mesh_region border_faces;
	outer_faces_of_mesh(mesht, border_faces);

	for (mr_visitor i(border_faces); !i.finished(); ++i) {

		assert(i.is_face());

		// Unit outward normal : used to identify faces
		//! \todo Use getfem 5.0's function select_faces_of_normal?
		base_node un = mesht.normal_of_face_of_convex(i.cv(), i.f());
		un /= gmm::vect_norm2(un);

		if (gmm::abs(un[0] + 1.0) < 1.0E-7)      // back
			mesht.region(0).add(i.cv(), i.f());
		else if (gmm::abs(un[0] - 1.0) < 1.0E-7) // front
			mesht.region(1).add(i.cv(), i.f());
		else if (gmm::abs(un[1] + 1.0) < 1.0E-7) // left
			mesht.region(2).add(i.cv(), i.f());
		else if (gmm::abs(un[1] - 1.0) < 1.0E-7) // right
			mesht.region(3).add(i.cv(), i.f());
		else if (gmm::abs(un[2] + 1.0) < 1.0E-7) // bottom
			mesht.region(4).add(i.cv(), i.f());
		else if (gmm::abs(un[2] - 1.0) < 1.0E-7) // top
			mesht.region(5).add(i.cv(), i.f());
		
	} /* end of border_faces loop */
	// Export an indicator function for BC regions
	if (PARAM.int_value("VTK_EXPORT")){
		vector_type ones(dof.coeft(), 1.0);
		vector_type indicator(dof.coeft());
		for (unsigned f=0; f<2*DIMT; ++f) {
			asm_source_term(indicator, mimt, mf_coeft, mf_coeft, 
				gmm::scaled(ones, BCt[f].value), mesht.region(BCt[f].rg));
		}
		vtk_export rgvtk(descr.OUTPUT+"mesht_boundary.vtk");
		rgvtk.exporting(mf_coeft);
		rgvtk.write_mesh();
		rgvtk.write_point_data(mf_coeft, indicator, "1t");
	}
} 

void //-------------------- DA MODIFICARE LE CONDIZIONI AL CONTORNO ----------------//
root_problem3d1d::build_vessel_boundary(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Building vessel boundary ..." << endl;
	#endif
  try {

	dal::bit_vector junctions; // global idx of junctions vertices in meshv
	dal::bit_vector extrema;   // global idx of extreme vertices in meshv

	Jv.clear();
	nb_extrema=0; 
	nb_junctions=0;
	
	size_type fer = nb_branches; // first empty region
	GMM_ASSERT1(meshv.has_region(fer)==0, 
		"Overload in meshv region assembling!");
	
	// List all the convexes
	dal::bit_vector nn = meshv.convex_index();
	bgeot::size_type cv;
	for (cv << nn; cv != bgeot::size_type(-1); cv << nn) {
		
		bgeot::pconvex_structure cvs = meshv.structure_of_convex(cv);
		if (cvs->nb_points()>2) 
			cerr << "Error: convex #" << cv << "has more than 2 vertices!" << endl;
		if (cvs->nb_faces()>2)  
			cerr << "Error: convex #" << cv << "has more than 2 faces!" << endl;

		// Build regions for BCs and junctions
		// Global idx of mesh vertices
		size_type i0 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(1)[0]];
		size_type i1 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(0)[0]];
		// Identify vertex type
		if (meshv.convex_to_point(i0).size()==1){ /* inflow extremum */
			// Update information
			extrema.add(i0);
			nb_extrema++;
			// Build a new region made by a single face
			GMM_ASSERT1(meshv.has_region(fer)==0, 
				"Overload in meshv region assembling!");
			meshv.region(fer).add(cv, 1);
			// Store the current index and then update it
			size_type bc = 0; 
			bool found = false;
			while (!found && (bc<BCv.size())) {
				found = (i0 == BCv[bc].idx);
				if (!found) bc++;
			}
			GMM_ASSERT1(found=true, "Miss a boundary node in BCv list!");
			BCv[bc].rg = fer; 
			fer++;
			// Store the containing branch index
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			BCv[bc].branches.emplace_back(branch); 
		}
		else if (meshv.convex_to_point(i0).size()==2){ /* trivial inflow junction */
			// DO NOTHING
		}
		else if (meshv.convex_to_point(i0).size()>=2){ /* non-trivial inflow junction */
			// Check if jucntion has been already stored, 
			// if not add to the junction list (J) and build a new region
			dal::bit_vector tmp; tmp.add(i0);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i0);
				nb_junctions++;
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshv.region(fer).add(cv, 1); // single-face region
				// Create a new junction node
				Jv.emplace_back("JUN", 0, i0, fer);
				fer++;
			}
			// Search for index of containing branch (\mathcal{P}^{in}_j)
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			// Add the inflow branch (to the right junction node)
			size_type jj = 0;
			bool found = false;
			while (!found && jj < nb_junctions){
				found = (i0 == Jv[jj].idx);
				if (!found) jj++;
			}
			cout << "Branch -" << branch << " added to junction " << jj << endl;
			Jv[jj].value += r_param.R(mimv, branch);
			Jv[jj].branches.emplace_back(-branch);
			GMM_ASSERT1(branch>0, 
				"Error in network labeling: -0 makes no sense");
		}
		
		if (meshv.convex_to_point(i1).size()==1){ 
			size_type bc = 0; 
			bool found = false;
			while (!found && (bc<BCv.size())) {
				found = (i1 == BCv[bc].idx);
				if (!found) bc++;
			}
			if (found){ /* outlow extremum */
				extrema.add(i1); 
				nb_extrema++; 
				// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 0);
				// Store the current index and then update it
				BCv[bc].value *= -1.0;
				BCv[bc].rg = fer; 
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv[bc].branches.emplace_back(branch); 
			}
			else { /* interior -> Mixed point */
				// "MIX" label via post-processing
				// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 0);
				BCv.emplace_back("MIX", 0.0, i1, fer);
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv.back().branches.emplace_back(branch); 
			}
		}
		else if (meshv.convex_to_point(i1).size()==2){ /* trivial outflow junction */

			// Search for index of first containing branch (\mathcal{P}^{out}_j)
			size_type firstbranch = 0; 
			bool contained = false;
			while (!contained && firstbranch<nb_branches ) {
				contained = meshv.region(firstbranch).is_in(cv);
				if (!contained) firstbranch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i1!");

			// Check if i1 is a trivial junction (or a INT point)
			size_type cv1 = meshv.convex_to_point(i1)[0];
			size_type cv2 = meshv.convex_to_point(i1)[1];
			bool is_junc = (meshv.region(firstbranch).is_in(cv1) < 1 ||
							meshv.region(firstbranch).is_in(cv2) < 1 );
							
			if (is_junc){
				cout << "Found a trivial junction at i1 = " << i1 << endl;
				// Check if jucntion has been already stored, 
				// if not add to the junction list (J) and build a new region
				dal::bit_vector tmp; tmp.add(i1);
				if(!junctions.contains(tmp)){
					// Store the junction vertex
					junctions.add(i1);
					nb_junctions++;
					GMM_ASSERT1(meshv.has_region(fer)==0, 
						"Overload in meshv region assembling!");
					// Build a new region with idx "first empty region"
					meshv.region(fer).add(cv, 0);
					// Create a new junction node
					Jv.emplace_back("JUN", 0, i1, fer);
					fer++;
					//} QUI C'ERA LA CHIUSURA DI JUNCTION.CONTAINS(TMP)
					// Search for index of second containing branch (\mathcal{P}^{out}_j)
					size_type secondbranch = 0; 
					size_type secondcv = (( cv1 == cv) ? cv2 : cv1);
					size_type firstcv  = (( cv1 != cv) ? cv2 : cv1);
					contained = false;
					while (!contained && secondbranch<nb_branches ) {
						if(secondbranch != firstbranch)
							contained = meshv.region(secondbranch).is_in(secondcv);
						if (!contained) secondbranch++;
					}
					GMM_ASSERT1(contained=true, "No branch region contains node i1!");
					// Add the two branches
					scalar_type in=0;
					if(meshv.ind_points_of_convex(firstcv)[0]==i1) in=-1;
					else if (meshv.ind_points_of_convex(firstcv)[1]==i1) in=+1;
					GMM_ASSERT1(in!=0, "There is something wrong in firstbranch convex index");
					Jv.back().branches.emplace_back(in*firstbranch);

					in=0;
					if(meshv.ind_points_of_convex(secondcv)[0]==i1) in=-1;
					else if (meshv.ind_points_of_convex(secondcv)[1]==i1) in=+1;
					GMM_ASSERT1(in!=0, "There is something wrong in secondbranch convex index");
					Jv.back().branches.emplace_back(in*secondbranch);
					Jv.back().value += r_param.R(mimv, firstbranch);
					Jv.back().value += r_param.R(mimv, secondbranch);
				}
			}
		}
		else if (meshv.convex_to_point(i1).size()>=2){ /* non-trivial outflow junction */

			// Search for index of containing branch (\mathcal{P}^{out}_j)
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");

			// Check if jucntion has been already stored, 
			// if not add to the junction list (J) and build a new region
			dal::bit_vector tmp; tmp.add(i1);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i1);
				nb_junctions++;
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshv.region(fer).add(cv, 0);
				// Create a new junction node
				Jv.emplace_back("JUN", 0, i1, fer);
				// Add the outflow branch
				Jv.back().branches.emplace_back(+branch);
				Jv.back().value += r_param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << i1 << endl;
				fer++;
			}
			else {
				// Add the outflow branch (to the right junction node)
				size_type jj = 0;
				bool found = false;
				while (!found && jj < nb_junctions){
					found = (i1 == Jv[jj].idx);
					if (!found) jj++;
				}
				Jv[jj].branches.emplace_back(+branch);
				Jv[jj].value += r_param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << jj << endl;
			}
		}
	} /* end of convexes loop */
	
	// Ckeck network assembly
	#ifdef M3D1D_VERBOSE_
	cout << "--- NETWORK ASSEMBLY ------------------ "   << endl;
	cout << "  Branches:   " << nb_branches << endl
		 << "  Vertices:   " << nn.size()+1 << endl;
	cout << "  Extrema:    " << extrema << endl;	  
	for (size_type i=0; i<BCv.size(); ++i)
		cout << "    -  label=" << BCv[i].label 
			 << ", value=" << BCv[i].value << ", ind=" << BCv[i].idx 
			 << ", rg=" << BCv[i].rg << ", branches=" << BCv[i].branches << endl; 
	cout << "  Junctions: " << junctions << endl;
	for (size_type i=0; i<Jv.size(); ++i)
		cout << "    -  label=" << Jv[i].label 
			 << ", value=" << Jv[i].value << ", ind=" << Jv[i].idx 
			 << ", rg=" << Jv[i].rg << ", branches=" << Jv[i].branches << endl; 
	cout << "---------------------------------------- "   << endl;
	#endif

  } 
  GMM_STANDARD_CATCH_ERROR; // catches standard errors
} /* end of build_vessel_boundary */


//////// Assemble the problem ////////////////////////////////////////// 
void
root_problem3d1d::assembly(void)
{
	//1 Build the monolithic matrix AM
	assembly_mat();
	//2 Build the monolithic rhs FM
	assembly_rhs();
}

void //--- Rimuoverei il termine di Assemblaggio Mtt ---//
root_problem3d1d::assembly_mat(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Allocating AM, UM, FM ..." << endl;
	#endif
	gmm::resize(AM, dof.tot(), dof.tot()); gmm::clear(AM);
	gmm::resize(UM, dof.tot()); gmm::clear(UM);
	gmm::resize(FM, dof.tot()); gmm::clear(FM);
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic matrix AM ..." << endl;
	#endif
	// Divergence matrix for the interstitial problem
	sparse_matrix_type Dtt(dof.Pt(), dof.Ut());
	// Junction compatibility matrix for the network problem
	sparse_matrix_type Jvv(dof.Pv(), dof.Uv());
	// Tissue-to-tissue exchange matrix
	sparse_matrix_type Btt(dof.Pt(), dof.Pt());
	// Vessel-to-tissue exchange matrix
	sparse_matrix_type Btv(dof.Pt(), dof.Pv());
	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Bvt(dof.Pv(), dof.Pt());
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Bvv(dof.Pv(), dof.Pv());
	// Aux tissue-to-vessel averaging matrix
	sparse_matrix_type Mbar(dof.Pv(), dof.Pt());
	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type Mlin(dof.Pv(), dof.Pt());



	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Mtt and Dtt ..." << endl;
	#endif
	asm_tissue_darcy_divergence(Dtt, mimt, mf_Ut, mf_Pt);
	// Copy -Dtt^T
	gmm::add(gmm::scaled(gmm::transposed(Dtt), -1.0),  
			  gmm::sub_matrix(AM, 
					gmm::sub_interval(0, dof.Ut()),
					gmm::sub_interval(dof.Ut(), dof.Pt()))); 
	// Copy Dtt
	gmm::add(Dtt,
			  gmm::sub_matrix(AM, 
					gmm::sub_interval(dof.Ut(), dof.Pt()),
					gmm::sub_interval(0, dof.Ut()))); 

	#ifdef M3D1D_VERBOSE_ 
	cout << "  Assembling the tangent versor ..." << endl;
	#endif
 
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Mvv and Dvv ..." << endl;
	#endif
	// Local matrices
	size_type shift = 0;
	for(size_type i=0; i<nb_branches; ++i){

		if(i>0) shift += mf_Uvi[i-1].nb_dof();
		scalar_type Ri = r_param.R(mimv, i);
		scalar_type kvi= r_param.kv(mimv, i);
		// Coefficient \pi^2*Ri'^4/\kappa_v *(1+Ci^2*Ri^2) //Adaptation to the curve model
		vector_type ci(mf_coefvi[i].nb_dof());
		for(size_type j=0; j<mf_coefvi[i].nb_dof(); ++j){
			ci[j]=pi*pi*Ri*Ri*Ri*Ri/kvi*(1+(r_param.Curv())[i][j]*(r_param.Curv())[i][j]*Ri*Ri);
		}
		// Allocate temp local matrices
		sparse_matrix_type Mvvi(mf_Uvi[i].nb_dof(), mf_Uvi[i].nb_dof());
		sparse_matrix_type Dvvi(dof.Pv(), mf_Uvi[i].nb_dof());
 
 		// Build Mvvi and Dvvi
		asm_network_poiseuille(Mvvi, Dvvi, 
			mimv, mf_Uvi[i], mf_Pv, mf_coefvi[i],
			ci, (r_param.lambdax())[i], (r_param.lambday())[i], (r_param.lambdaz())[i], meshv.region(i));

		gmm::scale(Dvvi, pi*Ri*Ri);

		// Copy Mvvi and Dvvi
		gmm::add(Mvvi, 
			gmm::sub_matrix(AM, 
				gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof()), 
				gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof()))); 

		gmm::add(gmm::scaled(gmm::transposed(Dvvi), -1.0),
			gmm::sub_matrix(AM, 
				gmm::sub_interval(dof.Ut()+dof.Pt()+shift,    mf_Uvi[i].nb_dof()),
				gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()))); 

		gmm::add(Dvvi, 
			gmm::sub_matrix(AM, 
				gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()),
				gmm::sub_interval(dof.Ut()+dof.Pt()+shift,     mf_Uvi[i].nb_dof()))); 
		gmm::clear(Mvvi); 
		gmm::clear(Dvvi);	
	} /* end of branches loop */
	 
    if (nb_junctions > 0){
		#ifdef M3D1D_VERBOSE_
		cout << "  Assembling Jvv" << " ..." << endl;
		#endif 
		asm_network_junctions(Jvv, mimv, mf_Uvi, mf_Pv, mf_coefv,  
			Jv, r_param.R());
		#ifdef M3D1D_VERBOSE_
		cout << "  Copying -Jvv^T" << " ..." << endl;
		#endif		
		gmm::add(gmm::scaled(gmm::transposed(Jvv), -1.0),
			gmm::sub_matrix(AM,
				 gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv()), 
				 gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())));
		#ifdef M3D1D_VERBOSE_
		cout << "  Copying Jvv" << " ..." << endl;
		#endif		
		gmm::add(Jvv,
			gmm::sub_matrix(AM,
				 gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()),
				 gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));
    }
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling aux exchange matrices Mbar and Mlin ..." << endl;
	#endif
	asm_exchange_aux_mat(Mbar, Mlin, 
			mimv, mf_Pt, mf_Pv, r_param.R(), descr.NInt);
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling exchange matrices ..." << endl;
	#endif
	bool NEWFORM = PARAM.int_value("NEW_FORMULATION");
	asm_exchange_mat(Btt, Btv, Bvt, Bvv,
			mimv, mf_Pv, mf_coefv, Mbar, Mlin, r_param.Q(), NEWFORM);

	// Copying Btt
	gmm::add(Btt, 
			  gmm::sub_matrix(AM, 
					gmm::sub_interval(dof.Ut(), dof.Pt()), 
					gmm::sub_interval(dof.Ut(), dof.Pt()))); 
	// Copying -Btv
	gmm::add(gmm::scaled(Btv, -1),
	 		  gmm::sub_matrix(AM, 
					gmm::sub_interval(dof.Ut(), dof.Pt()),
					gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()))); 
	// Copying -Bvt
	gmm::add(gmm::scaled(Bvt,-1), 
			  gmm::sub_matrix(AM, 
			  		gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv()	, dof.Pv()),
					gmm::sub_interval(dof.Ut(), dof.Pt()))); 
	// Copying Bvv
	gmm::add(Bvv, 
			  gmm::sub_matrix(AM, 
					gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()), 
					gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()))); 

	// De-allocate memory
	gmm::clear(Dtt); 
	gmm::clear(Mbar); gmm::clear(Mlin);
	gmm::clear(Btt);  gmm::clear(Btv);
	gmm::clear(Bvt);  gmm::clear(Bvv);
}

void //-----------Modifica temporanea rimosso Fg---------------
root_problem3d1d::assembly_rhs(void)//-Ho imposto il termine di gravità sul tessuto-//
{
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic rhs FM ... " << endl;
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Initializing RHS for FM ..." << endl;
	#endif
	// Right Hand Side for the interstitial problem 
	vector_type Ft(dof.Ut());
	// Right Hand Side for the network problem 
	vector_type Fv(dof.Uv());
	// Right Hand Sude for the interstitial gravity term
	vector_type Fg(dof.Ut());

	// Coefficients for tissue BCs
	scalar_type bcoef  = PARAM.real_value("BETA", "Coefficient for mixed BC");
	scalar_type p0coef = PARAM.real_value("P0"); // default: 0

	#ifdef M3D1D_VERBOSE_
	cout << "  Building tissue boundary term ..." << endl;
	#endif
	vector_type beta(dof.coeft(), 1.0/bcoef);
	vector_type P0(dof.coeft(), p0coef);
	vector_type P0_vel(dof.coefv(),p0coef);
	
	if (PARAM.int_value("TEST_RHS")) {
		#ifdef M3D1D_VERBOSE_
		cout << "  ... as the divergence of exact velocity ... " << endl;
		#endif
		assembly_tissue_test_rhs();
	}
	else {
		sparse_matrix_type Mtt(dof.Ut(), dof.Ut());
		asm_tissue_bc(Mtt, Ft, mimt, mf_Ut, mf_coeft, BCt, P0, beta);		gmm::add(Mtt, 
			gmm::sub_matrix(AM,
				gmm::sub_interval(0, dof.Ut()),
				gmm::sub_interval(0, dof.Ut())));
		gmm::add(Ft, gmm::sub_vector(FM, gmm::sub_interval(0, dof.Ut())));
		// De-allocate memory
		gmm::clear(Mtt); 
	}
	if(PARAM.int_value("Gravity")){
		asm_tissue_darcy_gravity(Fg,mimt,mf_Ut,r_param.rho()*r_param.g());
		gmm::add(Fg, gmm::sub_vector(FM, gmm::sub_interval(0, dof.Ut())));
	}

	#ifdef M3D1D_VERBOSE_
	cout << "  Building vessel boundary term ..." << endl;
	#endif
	sparse_matrix_type Mvv(dof.Uv(), dof.Uv());
	my_asm_network_bc(Mvv, Fv, 
			mimv, mf_Uvi, mf_coefv, BCv, P0_vel, r_param.R(), bcoef);
	gmm::add(Mvv, 
		gmm::sub_matrix(AM,
			gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv()),
			gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));
	gmm::add(Fv, gmm::sub_vector(FM, gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));
	// De-allocate memory
	gmm::clear(Ft); gmm::clear(Fv); gmm::clear(Mvv); gmm::clear(Fg);
}

void
root_problem3d1d::assembly_tissue_test_rhs(void)
{
	// Exact rhs (see defines.hpp)
	vector_type sol_Gt(dof.coeft());
	interpolation_function(mf_coeft, sol_Gt, sol_gt);
	#ifdef M3D1D_VERBOSE_
	cout << "    Assemble divergence source term ... " << endl;
	#endif
	vector_type Gt(dof.Pt());
	asm_source_term(Gt, mimt, mf_Pt, mf_coeft, sol_Gt);	
	gmm::add(Gt, gmm::sub_vector(FM, gmm::sub_interval(dof.Ut(), dof.Pt())));

  if (PARAM.int_value("VTK_EXPORT")) {

	#ifdef M3D1D_VERBOSE_
	cout << "    Compute theoretical expectations ... " << endl;
	#endif
	// FE spaces for exact velocity
	mesh_fem mf_data(mesht), mf_data_vec(mesht);
	mf_data_vec.set_qdim(DIMT);
	bgeot::pgeometric_trans pgt_t = bgeot::geometric_trans_descriptor(descr.MESH_TYPET);
	mf_data.set_classical_discontinuous_finite_element(1);
	mf_data_vec.set_classical_discontinuous_finite_element(1);

	GMM_ASSERT1(mf_data.nb_dof()*DIMT==mf_data_vec.nb_dof(), 
		"Wrong size of mf_data_vec"); 
	// Exact pressure (see defines.hpp)
	vector_type sol_Pt(dof.coeft());
	interpolation_function(mf_coeft, sol_Pt, sol_pt);
	// Exact x velocity (see defines.hpp)
	vector_type sol_Utx(mf_data.nb_dof());
	interpolation_function(mf_data, sol_Utx, sol_utx);
	// Exact y velocity (see defines.hpp)
	vector_type sol_Uty(mf_data.nb_dof());
	interpolation_function(mf_data, sol_Uty, sol_uty);
	// Exact z velocity (see defines.hpp)
	vector_type sol_Utz(mf_data.nb_dof());
	interpolation_function(mf_data, sol_Utz, sol_utz);
	// Exact velocity magnitude (see defines.hpp)
	vector_type sol_Utm(mf_data.nb_dof());
	interpolation_function(mf_data, sol_Utm, sol_utm);
	// Exact vectorial velocity (see defines.hpp)
	vector_type sol_Ut; sol_Ut.reserve(mf_data_vec.nb_dof());
	for( size_type i=0; i<mf_data_vec.nb_dof()/DIMT; ++i ){
		sol_Ut.emplace_back(sol_Utx[i]);
		sol_Ut.emplace_back(sol_Uty[i]);
		sol_Ut.emplace_back(sol_Utz[i]);
	}

	#ifdef M3D1D_VERBOSE_
	cout << "    Export theoretical expectations ... " << endl;
	#endif
	vtk_export vtk_sol_Gt(descr.OUTPUT+"sol_Gt.vtk");
	vtk_sol_Gt.exporting(mf_coeft);
	vtk_sol_Gt.write_mesh();
	vtk_sol_Gt.write_point_data(mf_coeft, sol_Gt, "sol_Gt");

	vtk_export vtk_sol_Pt(descr.OUTPUT+"sol_Pt.vtk");
	vtk_sol_Pt.exporting(mf_coeft);
	vtk_sol_Pt.write_mesh();
	vtk_sol_Pt.write_point_data(mf_coeft, sol_Pt, "sol_Pt");

	vtk_export vtk_sol_Utm(descr.OUTPUT+"sol_Utm.vtk");
	vtk_sol_Utm.exporting(mf_data);
	vtk_sol_Utm.write_mesh();
	vtk_sol_Utm.write_point_data(mf_data, sol_Utm, "sol_Utm");

	vtk_export vtk_sol_Ut(descr.OUTPUT+"sol_Ut.vtk");
	vtk_sol_Ut.exporting(mf_data_vec);
	vtk_sol_Ut.write_mesh();
	vtk_sol_Ut.write_point_data(mf_data_vec, sol_Ut, "sol_Ut");

  }
}

void 
root_problem3d1d::assembly_nonlinear_mat(size_type iter, size_type tempo)
{
	//Setting all constant
	scalar_type A=r_param.A_ret();
	scalar_type gam=r_param.gamma_ret();
	scalar_type al;
	scalar_type be;
	scalar_type Tr;
	scalar_type Ts;
	scalar_type kt=r_param.kt(0);
	scalar_type rho = r_param.rho();
	scalar_type g= r_param.g();
	vector_type Ct(dof.coeft());
	vector_type Bt(dof.coeft());

	bool TIME_STEP=PARAM.int_value("SOLVE_TIME_STEP");

	if(TIME_STEP){
		al=r_param.alfa_ret();
		be=r_param.beta_ret();
		Tr=r_param.Theta_r();
		Ts=r_param.Theta_s();
	}
	// Preparing vectors dimension
    gmm::resize(Pt_data,dof.coeft()); gmm::clear(Pt_data);
    //Assembling conductivity non linear term
    if(iter==0 && tempo==0){
		asm_tissue_non_linear_darcy	(NLMtt,mimt,mf_Ut,mf_coeft,Ct,kt,iter);
    }
    else{
	    // Interpolating Pt_old on mf_coeft
    	getfem::interpolation(mf_Pt,mf_coeft,Pt_old_ITER,Pt_data);

    	for(size_type DOF=0; DOF<dof.coeft(); DOF++){
    		Ct[DOF]=(A+pow(abs(Pt_data[DOF])/rho/g, gam)) / A ;
    	}
    	asm_tissue_non_linear_darcy(NLMtt,mimt,mf_Ut,mf_coeft,Ct,kt,iter);
    	cout<<"1/Kt    Valore Nuovo "<<Ct[0]/kt<<" Valore Vecchio "<<1.0/kt<<"\n";
    	if(TIME_STEP){
	    	for(size_type DOF=0;DOF<dof.coeft();DOF++){
	    		//Bt=- a*(Ts-Tr)*|psi|^(b-1)/(a+|psi|^b)^2*sign(psi)
	    		Bt[DOF]=(Ts-Tr)/(al+pow(abs(Pt_data[DOF])/rho/g, be))/(al+pow(abs(Pt_data[DOF])/rho/g, be));
	    		Bt[DOF]=-Bt[DOF]*al*pow(abs(Pt_data[DOF])/rho/g, be-1.0)*abs(Pt_data[DOF])/Pt_data[DOF];    		
	    	}
	    	asm_tissue_non_linear_richards(NLRtt,NLFtt,mimt, mf_Pt,mf_coeft,Pt_Prev,Bt,dT*rho*g);
	    	cout<<"(Ts-Tr) Valore Nuovo "<<Bt[0]<<" Valore Vecchio "<<Ts-Tr<<"\n";
	    }
	}

    gmm::clear(Pt_data);
    gmm::clear(Ct);
    gmm::clear(Bt);
}

////////// Solve the problem ///////////////////////////////////////////    
bool//
root_problem3d1d::solve_iter(size_type iter)
{
	#ifdef M3D1D_VERBOSE_
	cout << "  Solving the monolithic system ... " << endl;
	#endif
	
	// Using an alternative matrix to which add Non Linear Term at every iteration
	gmm::copy(AM,CM);
	gmm::copy(FM,CFM);
			//Adding non linear matrix term
	gmm::add( NLMtt, 
		gmm::sub_matrix(CM, 
			gmm::sub_interval(0, dof.Ut()), 
			gmm::sub_interval(0, dof.Ut()))); 
	if(PARAM.int_value("SOLVE_TIME_STEP")){
		gmm::add(NLRtt,
			gmm::sub_matrix(CM,
				gmm::sub_interval(dof.Ut(),dof.Pt()),
				gmm::sub_interval(dof.Ut(),dof.Pt())));
		gmm::add(NLFtt,
			gmm::sub_vector(CFM,
				gmm::sub_interval(dof.Ut(),dof.Pt())));
	}
	gmm::csc_matrix<scalar_type> A;
	gmm::clean(CM, 1E-12);
	gmm::copy(CM, A);
	gmm::clear(CM);
	gmm::clear(NLMtt);

	if ( descr.SOLVE_METHOD == "SuperLU" ) { // direct solver //
		#ifdef M3D1D_VERBOSE_
		//cout << "  Applying the SuperLU method ... " << endl;
		#endif
		scalar_type cond;
		gmm::SuperLU_solve(A, UM, CFM, cond);
		//cout << "  Condition number : " << cond << endl;
	}
	else { // Iterative solver //

		// Iterations
		gmm::iteration iter(descr.RES);  // iteration object with the max residu
		iter.set_noisy(1);               // output of iterations (2: sub-iteration)
		iter.set_maxiter(descr.MAXITER); // maximum number of iterations

		// Preconditioners
		//! \todo Add preconditioner choice to param file
		// See \link http://download.gna.org/getfem/html/homepage/gmm/iter.html
		gmm::identity_matrix PM; // no precond
		//gmm::diagonal_precond<sparse_matrix_type> PM(AM); // diagonal preocond
		//gmm::ilu_precond<sparse_matrix_type> PM(AM);
		// ...
		//gmm::clear(AM);
		// See <http://download.gna.org/getfem/doc/gmmuser.pdf>, pag 15
	
		if ( descr.SOLVE_METHOD == "CG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Conjugate Gradient method ... " << endl;
			#endif
			gmm::identity_matrix PS;  // optional scalar product
			gmm::cg(AM, UM, FM, PS, PM, iter);
		}
		else if ( descr.SOLVE_METHOD == "BiCGstab" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the BiConjugate Gradient Stabilized method ... " << endl;
			#endif
			gmm::bicgstab(AM, UM, FM, PM, iter);
		}
		else if ( descr.SOLVE_METHOD == "GMRES" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Generalized Minimum Residual method ... " << endl;
			#endif
			size_type restart = 50;
			gmm::gmres(A, UM, FM, PM, restart, iter);
		}
		else if ( descr.SOLVE_METHOD == "QMR" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Quasi-Minimal Residual method ... " << endl;
			#endif
			gmm::qmr(AM, UM, FM, PM, iter);
		}
		else if ( descr.SOLVE_METHOD == "LSCG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the unpreconditionned Least Square CG method ... " << endl;
			#endif
			gmm::least_squares_cg(AM, UM, FM, iter);
		}
		// Check
		if (iter.converged())
			cout << "  ... converged in " << iter.get_iteration() << " iterations." << endl;
		else if (iter.get_iteration() == descr.MAXITER)
			cerr << "  ... reached the maximum number of iterations!" << endl;
	}
	#ifdef M3D1D_VERBOSE_
	//cout << "Compute the total flow rate ... " << endl;
	#endif
	// Aux vector
	vector_type Uphi(dof.Pv()); 
	// Extracting matrices Bvt, Bvv
	sparse_matrix_type Bvt(dof.Pv(), dof.Pt());
	sparse_matrix_type Bvv(dof.Pv(), dof.Pv());
	gmm::copy(gmm::sub_matrix(A, 
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv()	, dof.Pv()),
			gmm::sub_interval(dof.Ut(), dof.Pt())),
				Bvt); 
	gmm::copy(gmm::sub_matrix(A, 
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()), 
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())),
				Bvv); 
	// Extracting solutions Pt, Pv 
	vector_type Pt(dof.Pt()); 
	vector_type Pv(dof.Pv()); 
	gmm::copy(gmm::sub_vector(UM, 
		gmm::sub_interval(dof.Ut(), dof.Pt())), Pt);
	gmm::copy(gmm::sub_vector(UM, 
		gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())), Pv);
	// Computing Bvv*Pv - Bvt*Pt
	gmm::mult(Bvt, Pt, Uphi);
	gmm::mult_add(Bvv, Pv, Uphi);
	TFR = std::accumulate(Uphi.begin(), Uphi.end(), 0.0);

	// De-allocate memory
	
	gmm::clear(Bvt); gmm::clear(Bvv);
	gmm::clear(Pt);  gmm::clear(Pv);  
	gmm::clear(Uphi); gmm::clear(CFM);
	
	return true;
}


bool  
root_problem3d1d::solve_time(size_type tempo)
{

	//Importing descriptors
	size_type Max_iter = PARAM.int_value("Max_iter");
	scalar_type minERR = PARAM.real_value("minERR");
	bool Saturated = PARAM.int_value("Saturated");

	// Starting Condition
	size_type N_iter = 0;	
	scalar_type Error  = minERR+1.0;
	double time = gmm::uclock_sec();	
	bool solved=false;
	size_type shift=0;
	size_type dofi=mf_Uvi[0].nb_dof(); 

	#ifdef M3D1D_VERBOSE_
	cout<<"Assembling non linear tissue matrix"<<endl;
	#endif
	
	assembly_nonlinear_mat(N_iter,tempo);

	#ifdef M3D1D_VERBOSE_
	cout<<"Solving iter=0"<<endl;
	#endif
	
	if(!(solve_iter(N_iter))) 
		GMM_ASSERT1(0," At the iteration "<<N_iter<<" at time step "<<dT*tempo<<" the solver stopped");

	if(Saturated)
		solved=true;

	#ifdef M3D1D_VERBOSE_
	cout<<endl;
	#endif

	// Update of iter, old pressure in the tissue
	N_iter++;
	gmm::copy(gmm::sub_vector(UM,
					gmm::sub_interval( dof.Ut(), dof.Pt() ) ), 
				Pt_old_ITER);
	while( Error>minERR  && N_iter<Max_iter && !Saturated ) {   
		#ifdef M3D1D_VERBOSE_
		cout<<"Solving iter "<<N_iter<<"  with minIncremental Error "<<minERR <<endl;
		cout<<"  Assembling nonlinear term at iter "<<N_iter<<endl;
		#endif
		assembly_nonlinear_mat(N_iter,tempo);

		//Solving at Iteration:= N_iter
		if(!(solve_iter(N_iter))) 
			GMM_ASSERT1(0," At the step "<<N_iter<<" the solver stopped");
		
		//Computing of relative Incremental Error
		Error=asm_L2_dist(mimt, 
				mf_Pt , gmm::sub_vector(UM,gmm::sub_interval(dof.Ut(), dof.Pt())), 
				mf_Pt , Pt_old_ITER )
			    /asm_L2_norm(mimt, mf_Pt,gmm::sub_vector(UM,gmm::sub_interval(dof.Ut(),dof.Pt())));

		#ifdef M3D1D_VERBOSE_
		cout<<"Increment at iter "<<N_iter<<" of time step "<<dT*tempo<<" : "<<Error<<endl<<endl;
		#endif

		//Update
		N_iter++;
		gmm::clear(Pt_old_ITER);	
		gmm::copy(gmm::sub_vector(UM,gmm::sub_interval( dof.Ut(), dof.Pt() ) ), Pt_old_ITER);
	}

	double time_end=gmm::uclock_sec() - time;

	if((Error>minERR )&&(Max_iter>1 )&& !Saturated)
	{
		#ifdef M3D1D_VERBOSE_
		cout<<"*---------------------------------------------------------------------------*"<<endl;
		cout<<" The method doesn't converged for minumum Increment Error="<<minERR<<" at Timestep= "<<dT*tempo<<endl;
		cout<<"    Increment Error="<<Error<<endl;
		cout<<"*---------------------------------------------------------------------------*"<<endl<<endl;
		#endif
		solved= false;
	} else{  
		#ifdef M3D1D_VERBOSE_
		cout<<"Method converged in "<<N_iter<<" iterations for time step "<<dT*tempo<<" with L2: error "<<Error<<"   *"<<endl;
		cout << "* Time to solve the time iteration : " << time_end << " seconds               *"<<endl;
		#endif
		solved= true;
	}

	return solved;
}

bool 
root_problem3d1d::solve(){
	#ifdef M3D1D_VERBOSE_
	cout << "Allocating Non linear Term: NLMtt..." << endl;
	#endif
	gmm::resize(NLMtt, dof.Ut(),dof.Ut()); gmm::clear(NLMtt);

	#ifdef M3D1D_VERBOSE_
	cout << "Allocating Auxiliar Term CM,Pt_old_ITER..." << endl;
	#endif

	gmm::resize(CM, dof.tot(),dof.tot());  gmm::clear(CM);
	gmm::resize(CFM,dof.tot());			   gmm::clear(CFM);
	gmm::resize(Pt_old_ITER,dof.Pt());     gmm::clear(Pt_old_ITER);

	bool TIME_STEP=PARAM.int_value("SOLVE_TIME_STEP");
	bool EXPORT_TIME_STEP=PARAM.int_value("VTK_EXPORT_TIME");

	double time_ = gmm::uclock_sec();	

	if(TIME_STEP){
		gmm::resize(NLRtt,dof.Pt(),dof.Pt()); 	gmm::clear(NLRtt);
		gmm::resize(NLFtt,dof.Pt());			gmm::clear(NLFtt);
		gmm::resize(Pt_Prev,dof.Pt());			gmm::clear(Pt_Prev);
		dT=PARAM.real_value("TIME_STEP","Give the increment time step value");
		maxT=PARAM.real_value("MAX_TIME","Give the maximum interval of time in which solve the problem");
	}
	else{
		dT=3;
		maxT=1;
	}

	size_type Tempo=0;
	bool solved=true;

	while((Tempo*dT<maxT)&&(solved)){
		solved=solve_time(Tempo);
		if(solved){
			gmm::clear(Pt_Prev);	
			gmm::copy(gmm::sub_vector(UM,gmm::sub_interval( dof.Ut(), dof.Pt() ) ), Pt_Prev);
			if(EXPORT_TIME_STEP && TIME_STEP){
				//string suff="time_"+std::on_string(tempo);
				export_vtk("time_"+std::to_string(Tempo)+"_");
			}
			Tempo++;
		}
	}
	double time_end = gmm::uclock_sec()-time_;	
	if(!solved){
		#ifdef M3D1D_VERBOSE_
		cout<<"Il progamma si è arrestato allo step "<<Tempo*dT<<" su "<<maxT<<"\n\n";
		#endif
	}
	else{
		#ifdef M3D1D_VERBOSE_
		cout<<"Il programma è riuscito a risolvere senza problemi la simulazione\n Tempo impiegato "<<time_end<<" secondi \n\n";
		#endif
	}

			 
	gmm::clear(AM);
	gmm::clear(FM);
	gmm::clear(Pt_old_ITER);
	gmm::clear(NLMtt);
	gmm::clear(CM);
	gmm::clear(CFM);
	gmm::clear(Pt_Prev);
	gmm::clear(NLFtt);
	gmm::clear(NLRtt);
	return solved;
}


////////// Export results into vtk files ///////////////////////////////
void 
root_problem3d1d::export_vtk(const string & suff)
{
  if (PARAM.int_value("VTK_EXPORT"))
  {
		#ifdef M3D1D_VERBOSE_
		cout << "Exporting the solution (vtk format) to " << descr.OUTPUT << " ..." << endl;
		#endif
		#ifdef M3D1D_VERBOSE_
		cout << "  Saving the results from the monolithic unknown vector ... " << endl;
		#endif
		// Array of unknown dof of the interstitial velocity
		vector_type Ut(dof.Ut()); 
		// Array of unknown dof of the interstitial pressure
		vector_type Pt(dof.Pt()); 
		// Array of unknown dof of the network velocity
		vector_type Uv(dof.Uv()); 
		// Array of unknown dof of the network pressure
		vector_type Pv(dof.Pv()); 
		gmm::copy(gmm::sub_vector(UM, 
			gmm::sub_interval(0, dof.Ut())), Ut);
		gmm::copy(gmm::sub_vector(UM, 
			gmm::sub_interval(dof.Ut(), dof.Pt())), Pt);
		gmm::copy(gmm::sub_vector(UM, 
			gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())), Uv);
		gmm::copy(gmm::sub_vector(UM, 
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())), Pv);

		#ifdef M3D1D_VERBOSE_
		// Save vessel solution for test-cases
		if (nb_branches==1){
			std::ofstream outUv("Uv.txt");
			outUv << gmm::col_vector(Uv);
			outUv.close();
			std::ofstream outPv("Pv.txt");
			outPv << gmm::col_vector(Pv);
			outPv.close();
		}
		#endif
		#ifdef M3D1D_VERBOSE_
		cout << "  Exporting Ut ..." << endl;
		#endif
		pfem pf_Ut = fem_descriptor(descr.FEM_TYPET);
		if(pf_Ut->is_lagrange()==0){ 
			/*
				There is no built-in export for non-lagrangian FEM.
				If this is the case, we need to project before exporting.
			 */
			#ifdef M3D1D_VERBOSE_
			cout << "    Projecting Ut on P1 ..." << endl;
			#endif
			mesh_fem mf_P1(mesht);
			mf_P1.set_qdim(bgeot::dim_type(DIMT)); 
			mf_P1.set_classical_finite_element(1);
			sparse_matrix_type M_RT0_P1(mf_P1.nb_dof(), dof.Ut());
			sparse_matrix_type M_P1_P1(mf_P1.nb_dof(), mf_P1.nb_dof());
			vector_type Ut_P1(mf_P1.nb_dof());
			asm_mass_matrix(M_RT0_P1, mimt, mf_P1, mf_Ut);
			asm_mass_matrix(M_P1_P1,  mimt, mf_P1, mf_P1);
			
			vector_type Utt(mf_P1.nb_dof());
			gmm::mult(M_RT0_P1, Ut, Utt);
			double cond1;
			gmm::SuperLU_solve(M_P1_P1, Ut_P1, Utt, cond1);

			vtk_export exp1(descr.OUTPUT+"Ut"+suff+".vtk");
			exp1.exporting(mf_P1);
			exp1.write_mesh();
			exp1.write_point_data(mf_P1, Ut_P1, "Ut");
		}	
		else {
			vtk_export exp_Ut(descr.OUTPUT+"Ut"+suff+".vtk");
			exp_Ut.exporting(mf_Ut);
			exp_Ut.write_mesh();
			exp_Ut.write_point_data(mf_Ut, Ut, "Ut");	 
		}
		#ifdef M3D1D_VERBOSE_
		cout << "  Exporting Pt ..." << endl;
		#endif
		vtk_export exp_Pt(descr.OUTPUT+"Pt"+suff+".vtk");
		exp_Pt.exporting(mf_Pt);
		exp_Pt.write_mesh();
		exp_Pt.write_point_data(mf_Pt, Pt, "Pt");

		#ifdef M3D1D_VERBOSE_
		cout << "  Exporting Uv ..." << endl;
		#endif
		size_type start = 0;
		size_type length = 0;
		for(size_type i=0; i<nb_branches; ++i){
			if(i>0) start += mf_Uvi[i-1].nb_dof();
			length = mf_Uvi[i].nb_dof();
			vtk_export exp_Uv(descr.OUTPUT+"Uv"+suff+std::to_string(i)+".vtk");
			exp_Uv.exporting(mf_Uvi[i]);
			exp_Uv.write_mesh();
			exp_Uv.write_point_data(mf_Uvi[i], 
				gmm::sub_vector(Uv, gmm::sub_interval(start, length)), "Uv"); 
		}

		#ifdef M3D1D_VERBOSE_
		cout << "  Exporting Pv ..." << endl;
		#endif
		vtk_export exp_Pv(descr.OUTPUT+"Pv"+suff+".vtk");
		exp_Pv.exporting(mf_Pv);
		exp_Pv.write_mesh();
		exp_Pv.write_point_data(mf_Pv, Pv, "Pv");

		#ifdef M3D1D_VERBOSE_
		cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
		#endif
     }
} // end of vtk

  
 
} /* end of namespace */

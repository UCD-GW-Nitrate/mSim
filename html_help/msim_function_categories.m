%% Functions by Category
% 
% <msim_help_main.html | main>   <msim_help_demos.html | Tutorials> 
% <msim_function_categories.html | Functions> <http://www.subsurface.gr | website> |
%
% mSim Toolbox
% Version 1.0 28-March-2014
%
%% mFlow
% * <Assemble_LHS_help.html |Assemble_LHS|> - Assemble the condunctunce matrix
% * <Assemble_RHS_help.html |Assemble_RHS|> - Assemble the right hand side
%
% * <calcBKB_help.html |caclBKB|> - Computes the product B'*K*B
% * <calcBKBline_Lin_help.html |calcBKBline_Lin|> - Computes the the product
% B'*K*B for linear line elements
% * <calcBKBline_quad_help.html |calcBKBline_quad|> - Computes the the product
% B'*K*B for linear line elements
% * <calcBKBtriang_Lin_help.html |calcBKBtriang_Lin|> - Computes the the product
% B'*K*B for linear trianlge elements
% * <calcBKBtriang_quad_help.html |calcBKBtriang_quad|> - Computes the the product
% B'*K*B for quadratic trianlge elements
% * <calcBKBQuad_Lin_help.html |calcBKBQuad_Lin|> - Computes the the product
% B'*K*B for linear quadrilateral elements
% * <calcBKBQuad_quad_9_help.html |calcBKBQuad_quad_9|> - Computes the the product
% B'*K*B for quadratic quadrilateral elements
% * <calcBKBPrism_Lin_help.html |calcBKBPrism_Lin|> - Computes the the product
% B'*K*B for linear prism elements
% * <calcBKBPrism_quad_help.html |calcBKBPrism_quad|> - Computes the the product
% B'*K*B for quadratic prism elements
% * <calcBKBHex_Lin_help.html |calcBKBHex_Lin|> - Computes the the product
% B'*K*B for linear hexahedral elements
% shape function derivatives of quadratic hexahedral elements
% * <calcBKBHex_quad_27_help.html |calcBKBHex_quad_27|> - Computes the the product
% B'*K*B for quadratic hexahedral elements
%
% * <initialize_BEM_help.html |initialize_BEM|> - Initialize the matrix that will hold the shape functions derivatives
% values
% * <integr_order_help.html |integr_order|> - Computes the
% gaussian integration points
% * <interp_with_shapefnc_help.html |interp_with_shapefnc|> - Interpolates
% using shape functions values
% 
% * <shapeDerivatives_help.html |shapeDerivatives|> - Computes the
% shape function derivatives
% * <shapeDerivLine_Lin_help.html |shapeDerivLine_Lin|> - Computes the
% shape function derivatives of linear line elements
% * <shapeDerivLine_quad_help.html |shapeDerivLine_quad|> - Computes the
% shape function derivatives of quadratic line elements
% * <shapeDerivTriang_Lin_help.html |shapeDerivTriang_Lin|> - Computes the
% shape function derivatives of linear trianlge elements
% * <shapeDerivTriang_quad_help.html |shapeDerivTriang_quad_help|> - Computes the
% shape function derivatives of quadratic trianlge elements
% * <shapeDerivQuad_Lin_help.html |shapeDerivQuad_Lin|> - Computes the
% shape function derivatives of linear quadrilateral elements
% * <shapeDerivQuad_quad_9_help.html |shapeDerivQuad_quad_9|> - Computes the
% shape function derivatives of quadratic quadrilateral elements
% * <shapeDerivPrism_Lin_help.html |shapeDerivPrism_Lin|> - Computes the
% shape function derivatives of linear prism elements
% * <shapeDerivPrism_quad_help.html |shapeDerivPrism_quad|> - Computes the
% shape function derivatives of quadratic prism elements
% * <shapeDerivHex_Lin_help.html |shapeDerivHex_Lin|> - Computes the
% shape function derivatives of linear hexahedral elements
% * <shapeDerivHex_quad_27_help.html |shapeDerivHex_quad_27|> - Computes the
% shape function derivatives of quadratic hexahedral elements
%
% * <shapefunctions_help.html |shapefunctions|> - Computes the
% shape functions
% * <solve_system_help.html |solve_system|> - Solves the linear system
%
%% mTrans
% * <ADE1Danalytical_help.html |ADE1Danalytical|> - Analytical solution of the ADE in 1D
% * <Assemble_LHS_std_help.html |Assemble_LHS_std|> - Assembles the
% matrices for ADE in steady state flow
%
% * <calcNvB_help.html |calcNvB|> - Computes the product N*v*B in vectorized manner
% * <calcNDBline_Lin_help.html |calcNDBline_Lin|> - Computes the product N*v*B for linear line elements
%
% * <calcNLN_help.html |calcNLN|> - Computes the product N'*L*N in
% vectorized manner for linear line elements
% * <DispCoeff_help.html |DispCoeff|> - Calculates the dispersion coefficients
% * <SteadyFlowTransport_help.html |SteadyFlowTransport|> - Solves the Advection 
% dispersion equation for steady state flow
%
%% mPart
% * <ParticleTracking_main_help.html |ParticleTracking_main|> - main
% particle tracking function
% options
% * <part_options_help.html |part_options|> - returns the particle tracking
% options
%
%% mUtil
% * <Build2Dmeshinfocpp_help.html |Build2Dmeshinfocpp|> - Generates the
% connectivity mesh
% * <Calc_Area_help.html |Calc_Area|> - Calculates the area of 2D mesh
% elements
% * <Calc_Barycenters_help.html |Calc_Barycenters|> - Calculates the
% volume of mesh elements
% * <Calc_volume_help.html |Calc_volume|> - Calculates the
% barycentes of mesh elements
% * <Centerfor2points_help.html |Centerfor2points|> - Computes the center 
% of a circle
% * <CSGobj_v2_help.html |CSGobj_v2|> - A class for the mesh
% generation
% * <Dist_Point_LineSegment_help.html |Dist_Point_LineSegment|> - Computes
% Point line distance
% * <extrude_mesh_help.html |extrude_mesh|> - Extrudes 2D mesh to 3D mesh
% * <find_elem_id_point_help.html |find_elem_id_point|> - Finds the id of the
% element that contains a given point
% * <helix_axis_help.html |helix_axis|> - Generates points along a helix
% with axis defined by two points
% * <iscw_help.html |iscw|> - Checks the orientation of mesh elements
% * <mapElemto2D_help.html |mapElemto2D|> - Transform elements (used
% internally)
% * <msim_compile_help.html |msim_compile|> - Compiles the C/C++ functions
% * <msim_mesh_options_help.html |msim_mesh_options|> - Default mesh
% options
% * <msim_polysplit_help.html |msim_polysplit|> - reimplementation of
% polysplit
% * <read_2D_Gmsh_help.html |read_2D_Gmsh|> - Reads mesh file generated by Gmsh
% * <readArcGisASCIIfile_help.html |readArcGisASCIIfile|> - Reads ARCGIS raster data in ascii format
% * <WriteAscii4Raster_help.html |WriteAscii4Raster|> - Writes Data to
% ASCII ARCGIS fromat
% * <WriteVtkMesh_help.html |WriteVtkMesh|> - Writes Data to Vtk format

%% mNPSAT
% * <ComputeURF_help.html |ComputeURF|> - Computes the unit response function for a given streamline
% * <ConvoluteURF_help.html |ConvoluteURF|> - Convolutes the unit response function with loading  functions
% * <distribute_particle_streams_help.html |distribute_particle_streams|> - Distribute particles around streams
% * <distribute_particle_wells_help.html |distribute_particle_wells|> - Distribute particles around wells
%%
% 
% <msim_help_main.html | main>   <msim_help_demos.html | Tutorials> 
% <msim_function_categories.html | Functions> <http://www.subsurface.gr | website> |
%
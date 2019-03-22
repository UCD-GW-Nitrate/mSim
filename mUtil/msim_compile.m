function msim_compile()
% msim_compile
%
% This will compile all the C/C++ mSim functions.
% Make sure that you call this function on the mSim root folder
% Next will go into various folders and compile the c/c++ functions
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : http://groundwater.ucdavis.edu/msim
% Date 5-Sep_2014 
% Department of Land Air and Water
% University of California Davis

mroot = msim_root;
if exist('OCTAVE_VERSION')
    cd([mroot '/mPart/'])
    mkoctfile Part_Track_oct.cpp
    
    cd([mroot '/mNPSAT/'])
    mkoctfile calcbtc_oct.cpp
    
else
    cd([mroot '/mPart/'])
    mex Part_Track_mat.cpp
    
    cd([mroot '/mNPSAT/'])
    mex calcbtc_mat.cpp
    
    cd([mroot '/mUtil/'])
    mex Build2Dmeshinfocpp.cpp
    mex read_2D_mesh_cpp.cpp
    
    
end
cd(msim_root);

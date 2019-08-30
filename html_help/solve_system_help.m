%% solve_system
% 
% <msim_help_main.html | main>   <msim_help_demos.html | Tutorials> 
% <msim_function_categories.html | Functions> <http://www.subsurface.gr | website> |
%
%
% This function solves the linear system of equations Kglo*H = F using default matlab
% methods. Note that this function takes care of the boundary conditions.
% It can also calculate the secondary variables which are the F terms on
% the constant head nodes.
% 
%
% Version : 1.0
%
% Author : George Kourakos
%
% email: giorgk@gmail.com
%
% web : <https://gwt.ucdavis.edu/research-tools-and-applications/msim
% https://gwt.ucdavis.edu/research-tools-and-applications/msim>
%
% Date : 18-Mar-2014 | update : 30-Aug-2019
%
% Department of Land Air and Water
%
% University of California Davis
%
%% Usage
% H = solve_system(Kglo, H, F)
%%%
% or
%%%
% [Hnew, Fnew] = solve_system(Kglo, H, F)
%% Input
% _*Kglo*_: System matrix 
%
% _*H*_: Solution vector. This function will have scalar values on the
% nodes associated with dirichlet boundary conditions and nan on the
% unknown dofs
%
% _*F*_: The right hand side
%
%% Output
% _*H*_: The solution of the system
%
%%
% <msim_help_main.html | main>   <msim_help_demos.html | Tutorials> 
% <msim_function_categories.html | Functions> <http://www.subsurface.gr | website> |
%
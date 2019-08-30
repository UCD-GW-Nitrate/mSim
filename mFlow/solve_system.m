function [Hnew, Fnew] = solve_system(Kglo,H,F)
%  Hnew = solve_system(Kglo, H, F)
%  or
% [Hnew, Fnew] = solve_system(Kglo, H, F)
%
% This function solves the linear system of equations Kglo*H = F using default matlab
% methods. Note that this function takes care of the boundary conditions.
% It can also calculate the secondary variables which are the F terms on
% the constant head nodes
%
% Input
% Kglo : System matrix 
% H    : Solution vector. This function will have scalar values on the
%        nodes associated with dirichlet boundary conditions and nan on the
%        unknown dofs
% F    : The right hand side
%
% Output
% Hnew    : The solution of the system
% Fnew    : The solution of the secondary variable
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : http://groundwater.ucdavis.edu/msim
% Date : 14-Jan-2013
% Department of Land Air and Water
% University of California Davis
%
% see also Assemble_LHS, Assemble_RHS

%find which nodes are known (cnst) and which are unknonw(var)
id_cnst=find(~isnan(H)); 
id_var=find(isnan(H));
%partition the system
KK=Kglo(id_var,id_var);
GG=Kglo(id_var,id_cnst);
DD=H(id_cnst);
B=F(id_var)-GG*DD;

%solution
x=KK\B;
%x=LISs(KK,B,zeros(length(B),1),'bicgstab',[]);

H(id_var,1)=x;
Hnew = H;
if nargout == 1
    return;
end

KK1=Kglo(id_cnst,id_cnst);
GG1=Kglo(id_cnst,id_var);
Fnew = nan(length(H),1);
Fnew(id_cnst,1) = KK1*DD + GG1*H(id_var,1);
Fnew(id_var) = F(id_var);

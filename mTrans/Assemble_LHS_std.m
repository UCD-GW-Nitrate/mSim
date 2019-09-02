function [Dglo Mglo c]= Assemble_LHS_std(p, MSH, A, V, rho_b, K_d, lambda, theta, Dm, CC, opt)
% [Dglo Mglo c]= Assemble_LHS_std(p, MSH, A, V, rho_b, K_d, lambda, theta, Dm, CC, opt)
%
% Assembles the global advection-dispersion and sorption matrix of the
% advection dispersion equation for steady state flow field.
%
% Inputs
% p     : coordinates of the Mesh nodes
% MSH   : Mesh structure
% A     : dispertivity [aL, aTH, aTV].
% V     : velocity
% rho_b : bulk density
% K_d   : equilibrium distribution coefficient
% lambda: decay constant
% theta : porosity
% Dm    : Molecular diffusion coefficient
% CC    : Dirichlet boundary conditions
% opt   : option structure with the following fields
%       .assemblemode : 'vect' for vectorize or 'nested' (just dont used the
%                       latter)
%       .capacmode : 'lumped' or 'consistent'
%       .int_ord : Gaussian intergration order. If the field does not exist a
%                  default value will be used according to element order
%                  and the dimension of the problem
%       .el_order : 'linear' Currently this is the only option for
%                    transport
%       .dim      : Dimension of the problem.
%       .el_type  : this is the element type.  
%
%
% Output
% Dglo  : Global advection dispersion matrix
% Mglo  : Global sorption matrix
% c     : The unkown concentration vector with nan values for unknown nodes 
%          scalar values for the dirichlet nodes
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : http://groundwater.ucdavis.edu/msim
% Date 28-Mar-2014
% Department of Land Air and Water
% University of California Davis
%


if ~isfield(opt, 'assemblemode')
    opt.assemblemode = 'vect';
end
  


Np=size(p,1);
switch opt.assemblemode
    case 'vect'
        [Dglo Mglo]=assemble_vect(p, MSH, A, V, rho_b, K_d, lambda, theta, Dm, opt);
    case 'nested'
        [Dglo Mglo]=assemble_nested_v2(p, MSH, D, v, opt);
end


c=nan(Np,1);
%Assign Constant Concentration 
if ~isempty(CC)
    c(CC(:,1))=CC(:,2);
end



function [Dglo Mglo]=assemble_vect(p, MSH, A, V, rho_b, K_d, lambda, theta, Dm, opt)
Nel=size(MSH,1);
Np=size(p,1);
pint=integr_order(opt);
Npint=size(pint,2);
BEM=initialize_BEM(Nel,opt);

if ~isfield(opt,'capacmode')
    opt.capacmode='lumped';
    consistent = false;
end

if strcmp(opt.capacmode,'consistent')
    Mel=BEM;
    consistent=true;
else
    opt.capacmode='lumped';
    consistent = false;
    Mel=zeros(Nel,1);
end

for ii=1:size(pint,1)
    [B, Jdet]=shapeDerivatives(p, MSH, pint(ii,1:Npint-1),opt);
    N=shapefunctions(pint(ii,1:Npint-1),opt);
    
    Ael=node2element(N,MSH,A,Np,Nel,opt);
    Vel=node2element(N,MSH,V,Np,Nel,opt);
    Del=DispCoeff(Ael,Vel,Dm);
    %assemble diffusion matrix
    BEM=BEM+bsxfun(@times,calcBKB(B,Del,opt,[]),Jdet).*pint(ii,Npint)+...
            bsxfun(@times,calcNvB(N,Vel,B,opt,[]),Jdet).*pint(ii,Npint);
        
    rho_b_el=node2element(N,MSH,rho_b,Np,Nel,opt);
    K_d_el=node2element(N,MSH,K_d,Np,Nel,opt);
    theta_el=node2element(N,MSH,theta,Np,Nel,opt);
    tempL=(1 + rho_b_el.*K_d_el./theta_el);
    
    if any(lambda~=0)
        lambda_el=node2element(N,MSH,lambda,Np,Nel,opt);
        L=lambda_el.*tempL;
        BEM=BEM+bsxfun(@times,calcNLN(N,L,opt),Jdet).*pint(ii,Npint);
    end
    if consistent
        Mel=Mel+bsxfun(@times,calcNLN(N,tempL,opt),Jdet).*pint(ii,Npint);
    else
        Mel=Mel+tempL.*pint(ii,Npint);
    end
    
end
if ~consistent
    Mel=Mel.*Jdet./sum(pint(:,Npint));
end
    


%fprintf('Assembling Conductance ... x out of %.0f\n',size(BEM,2))
Dglo=sparse(Np,Np);
Mglo=sparse(Np,Np);
Nsh=size(MSH,2);
for ii=1:Nsh
    for jj=1:Nsh
        i_lin=sub2ind([Nsh Nsh],jj,ii);
        %fprintf('                    %.0f\n',i_lin)
        Dglo=Dglo+sparse(MSH(:,ii),MSH(:,jj),BEM(:,i_lin),Np,Np);
        switch opt.capacmode
            case 'lumped'
                if ii==jj;
                    Mglo=Mglo+sparse(MSH(:,ii),MSH(:,jj),Mel,Np,Np);
                end
            case 'consistent'
                Mglo=Mglo+sparse(MSH(:,ii),MSH(:,jj),Mel(:,i_lin),Np,Np);
        end
    end
end





function Pel=node2element(N,MSH,Pnd,Np,Nel,opt)
Pel=zeros(Nel,size(Pnd,2));
if size(Pnd,1)==Np
    for jj=1:size(Pnd,2)
        Pel(:,jj)=interp_with_shapefnc(N,Pnd(:,jj),MSH,opt);
    end
else
    Pel=Pnd;
end

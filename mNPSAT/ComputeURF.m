function URF=ComputeURF(XYZ, Vxyz, opt)
% URF=ComputeURF(XYZ, Vxyz, opt)
% Computes the Unit Response Function for a given streamline. This function
% discretizes the streamline and solves the 1D Advection Dispersion
% Equation using a unit input of unit step at one side and returns the
% concentration at the other side of the streamline. Using the option
% structure parameters such as discretization size, lambda or longitudinal
% dispersivity can be configured.
%
%
% Input data:
% XYZ  : Streamline coordinates [Np x 3].
% Vxyz : Velocity of streamline points. 
% opt  : A structure with the following options
%        dx : discretization along streamline
%        dt : time discretization in years
%        Ttime : Total simulation time
%        Lmin : Minimum length for numerical solution. For streamlines with
%               length less than Lmin an analytical solution will be used. 
%        lambda.type : 1->lambda is scalar
%                    : 0->variable
%        lambda.val : if type == 1, lambda = lambda.val(1)
%                     if type == 0, lambda = lambda.val ( lambda.val has to
%                     be equal to Np)
%       aL : The longitudinal dispersivity is defined as function of
%            streamline length in the form of aL = alpha*L^beta, where L is
%            the streamline length and the two parameters are defined in a
%            structure with the following fields:
%       aL. alpha : is the alpha value. The default value if either the 
%                   aL or alpha field is not present in the structure is 
%                   0.32, which is based on leterature e.g. Neuman 1990
%                (https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/WR026i008p01749)
%       aL.beta : Similarly this is the beta of the function that describes the
%              longitudinal dispesivity. The default value is 0.83.
%
% Output data
% URF  : Unit response function
%
%
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : http://groundwater.ucdavis.edu/msim
% Date 28-Mar_2014
% Department of Land Air and Water
% University of California Davis
%
% IMPORTANT NOTE
% The first point in XYZ is the starting point of the backward particle 
% tracking near the CDS (e.g. the wells). 
%
% In FEM formulation we assume that the first point is near 
% the land surface, therefore we reverse the order inside the script.
%
%
%  |---v1----|---v2---|----v3-----|
%  -------------------------------------
%  x1       x2        x3         x4
%
if ~isfield(opt,'Lmin')
    opt.Lmin = 100;
end


T=(0:opt.dt:opt.Ttime)'*365;



Dm= 1.1578e-004;
P=[0;cumsum(sqrt(sum(diff(XYZ).^2,2)))];

if ~isfield(opt,'lambda')
    lambda = 0;
else
    if opt.lambda.type == 1
        lambda = opt.lambda.val(1);
    elseif opt.lambda.type == 0
        lambda = opt.lambda.val;
    end
end

if ~isfield(opt, 'aL')
    opt.aL.alpha = 0.32;
    opt.aL.beta = 0.83;
else
    if ~isfield(opt.aL, 'alpha')
        opt.aL.alpha = 0.32;
    end
    if ~isfield(opt.aL, 'beta')
        opt.aL.beta = 0.83;
    end
end


aL = opt.aL.alpha.*P(end).^opt.aL.beta;

if opt.dx > 0.1*aL
   opt.dx = 0.1*aL;
end

if P(end) < opt.Lmin
    lambda = mean(lambda);
    P=flipud(P(end)-P);
    dL = diff(P);
    v = flipud(sqrt(sum(Vxyz.^2,2)));
    age = [0;cumsum(dL./v(1:end-1))];
    v=P(end)./age(end);
    urf=ADE1Danalytical(P,T,v,1,aL,Dm,lambda,1);
    val = urf(:,end)';
    val(:,1)=[];
    negval=-1*val;
    negval=[zeros(1,1) negval];
    negval=negval(1,1:size(val,2)); 
    URF=val+negval;
    return
end
XYZ=flipud(XYZ);
Vxyz=flipud(Vxyz);
P=[0;cumsum(sqrt(sum(diff(XYZ).^2,2)))];


V=sqrt(sum(Vxyz.^2,2));
p_1d=0;Msh_1D=[];Np=1;Vel=[];
for i=1:length(P)-1
    Dx=P(i+1)-P(i);
    if Dx > opt.dx
        Nseg=round(Dx/opt.dx);
        x_temp=linspace(P(i),P(i+1),Nseg+1)';
        Np_seg=length(x_temp);
        el_temp=[(1:Np_seg-1)' (2:Np_seg)']+Np-1;
        p_1d=[p_1d;x_temp(2:end)];
        
    else
        p_1d=[p_1d;P(i+1)];
        if i == 1
            el_temp = [1 2];
        else
            el_temp = [Msh_1D(end,2) length(p_1d)];
        end
    end
    Msh_1D=[Msh_1D;el_temp];
    Np=length(p_1d);
    Vel=[Vel;V(i)*ones(size(el_temp,1),1)];
end
MSH(1,1).elem(1,1).type='Bndpnt';
MSH(1,1).elem(1,1).id=[1;Np];
MSH(2,1).elem(1,1).type='line';
MSH(2,1).elem(1,1).id=Msh_1D;
Nel=size(Msh_1D,1);
Np=length(p_1d);


opt.dim=1;
opt.el_type='line';
opt.el_order='linear';
opt.assemblemode='vect';
opt.capacmode='consistent';
CC=[1 1]; %boundary conditions

theta=ones(Nel,1);

K_d=0;
rho_b=1;%it doesnt matter since K_d=0

if isfield(opt,'lambda')
    if opt.lambda.type == 0
        lambda = interp1(P, lambda, p_1d, 'nearest');
    end
end
[Dglo Mglo c]= Assemble_LHS_std(p_1d, MSH(2,1).elem(1,1).id,...
    aL, Vel, rho_b, K_d, lambda, theta, Dm, CC, opt);

Cinit=zeros(Np,1);

F=sparse(Np,1);
wmega=0.5;%crank Nickolson scheme
C=SteadyFlowTransport(Mglo,Dglo,F,Cinit,T,c,wmega);
val=C(:,end)';
negval=-1*val;
negval=[zeros(1,1) negval];
negval=negval(1,1:size(val,2)); 
URF=val+negval;







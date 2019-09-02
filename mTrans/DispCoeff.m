function D=DispCoeff(A,V,Dm)
% D = DispCoeff(A, V, Dm)
%
% Calculates the dispersion coefficients for given longitudinal and
% transverse coefficient, velocity vectors and molecular diffusion
% coefficient
%
%
% Input
% A     :[aL aTH aTV] an array with the longitudinal and transverse
%         coefficients
% V     : [n x dim] velocity vectors. The code will assume the dimension of
%         the problem by the number of columns of velocity. For example if
%         two columns are given then two columns are expected for the A array
%         as well.
% Dm    : Molecular diffusion coefficient
%
%
% Output
% D     : Dispersivity coefficients.
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

%Calulate velocity magnitude
absV=sqrt(sum(V.^2,2));
if size(V,2)==1
    %In 1D
    %Dxx=a_L*v_x
    D=A.*absV+Dm;
elseif size(V,2)==2;
    %In 2D
    %Dxx=(a_T*v_y^2+a_L*v_x^2)/|v|
    %Dyy=(a_T*v_x^2+a_L*v_y^2)/|v|
    %Dxy=Dyx=(a_L-a_T)*v_x*v_y)/|v|
    D(:,1)=(A(:,2).*V(:,2).^2+A(:,1).*V(:,1).^2)./absV+Dm;
    D(:,2)=(A(:,1)-A(:,2)).*V(:,1).*V(:,2)./absV+Dm;
    D(:,3)=D(:,2);
    D(:,4)=(A(:,2).*V(:,1).^2+A(:,1).*V(:,2).^2)./absV+Dm;
elseif size(V,2)==3
    %In 3D
    D(:,1)=(A(:,1).*V(:,1).^2+A(:,2).*V(:,2).^2+A(:,3).*V(:,3).^2)/absV+Dm;%Dxx
    D(:,2)=(A(:,1)-A(:,2)).*V(:,1).*V(:,2)/absV+Dm;%Dxy
    D(:,3)=(A(:,1)-A(:,3)).*V(:,1).*V(:,3)/absV+Dm;%Dxz
    D(:,4)=D(:,2);
    D(:,5)=(A(:,1).*V(:,2).^2+A(:,2).*V(:,1).^2+A(:,3).*V(:,3).^2)/absV+Dm;%Dyy
    D(:,6)=(A(:,1)-A(:,3)).*V(:,2).*V(:,3)/absV+Dm;%Dyz
    D(:,7)=D(:,3);
    D(:,8)=D(:,6);
    D(:,9)=(A(:,1).*V(:,3).^2+A(:,3).*V(:,1).^2+A(:,3).*V(:,2).^2)/absV+Dm;%Dzz
end
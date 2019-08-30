function A = Calc_Area(p,MSH)
% V=Calc_Area(p,MSH)
%
% Computes the area of 2D mesh elements in a vectorized manner. It workds
% for triangular and quadrilateral elements of any order. However the sides
% of the elements must be striaght lines otherwise the calculation will be
% wrong.
%
% Input 
% p   :[Np x 2] coordinates
% MSH :[Nel, Nsh] Mesh element ids. Nel is the number of elements and
%       Nsh is the number of shapefunctions. The ids correspond to the row of the
%       array
% Output
% A : [Nel x 1] The area of each element
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : http://groundwater.ucdavis.edu/msim
% Date 22-Mar-2014 
% Department of Land Air and Water
% University of California Davis

[Nel Nsh] = size(MSH);

A=zeros(Nel,1);
if Nsh == 3 || Nsh == 6
    Np = 3;
elseif Nsh == 4 || Nsh >= 8
   Np = 4;
end
for i = 1:Np 
    if i == Np
        A = A + p(MSH(:,i),1).*p(MSH(:,1),2) - p(MSH(:,i),2).*p(MSH(:,1),1);
    else
        A = A + p(MSH(:,i),1).*p(MSH(:,i+1),2) - p(MSH(:,i),2).*p(MSH(:,i+1),1);
    end
end
A=A./2;
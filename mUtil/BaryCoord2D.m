function BC = BaryCoord2D( pnt, ND )
% BC = BaryCoord2D( pnt, ND ) 
%
% Calculates the barycentric coordinates of a trianlge or
%   quadrilateral. The function will return nan if the point is outside the
%   element
%
%   Input
% pnt   : [1x2] The coordinates of the point to calculate the barycentric
%               coordinates
% ND    : [3 or 4 x2] the coordinates of the triangle or the quadrilateral
%                     element
%
%
% Output 
% BC    : [3 or 4 x 1] The barycentric coordinates
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : https://gwt.ucdavis.edu/research-tools-and-applications/msim
% Date 22-Aug-2019 
% Department of Land Air and Water
% University of California Davis


% It make sense to calculate the barycentric coordinates of the point is
% inside the polygon. So the function will return nan if the point is
% outside
if ~inpolygon(pnt(1), pnt(2), ND(:,1), ND(:,2))
    BC= nan;
    return
end

if size(ND,1) == 3
    BC = BaryCoordsTriangle(pnt, ND);
elseif size(ND,1) == 4
    BC = BaryCoordsQuad(pnt, ND);
else
    error('Supported shapes are 2D triangles and quadrilaterals');
end


end

function BC = BaryCoordsTriangle(pnt, ND)
% https://gamedev.stackexchange.com/questions/33352/convert-point-to-its-barymetric-coordinates?noredirect=1&lq=1
% https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates/23745#23745?newreg=dbbc88792a27450abb5f03562392459f
v0 = ND(2,:) - ND(1,:);
v1 = ND(3,:) - ND(1,:);
v2 = pnt - ND(1,:);

den = v0(1)*v1(2) - v1(1)*v0(2);
v = (v2(1) * v1(2) - v1(1) * v2(2))/den;
w = (v0(1) * v2(2) - v2(1) * v0(2))/den;
u = 1 - v - w;

BC = [u;v;w];
end

function BC = BaryCoordsQuad(pnt, ND)
% https://numfactory.upc.edu/web/FiniteElements/Pract/P4-QuadInterpolation/html/QuadInterpolation.html
a = ND(1,:) - pnt;
b = ND(2,:) - ND(1,:);
c = ND(4,:) - ND(1,:);
d = ND(1,:) - ND(2,:) - ND(4,:) + ND(3,:);

A = cross([c,0], [d,0]); %must be 3D vectors
B = cross([c,0], [b,0]) + cross([a,0], [d,0]);
C = cross([a,0], [b,0]);

A=A(3);
B=B(3);
C=C(3);

% Check for unique solutions
if (abs(A)<1.e-14)
    u1= -C/B;
    u2=u1;
else
    % Check for non complex solutions
    if B^2-4*A*C >0
        u1=(-B+sqrt(B^2-4*A*C))/(2*A);
        u2=(-B-sqrt(B^2-4*A*C))/(2*A);
    else
        u1=-1000;
        u2=u1;
    end
end

mu=-10000;
if(u1>=0 && u1<=1)
    mu=u1;
end
if(u2>=0 && u2<=1)
    mu=u2;
end

A = cross([b,0], [d,0]); %must be 3D vectors
B = cross([b,0], [c,0]) + cross([a,0], [d,0]);
C = cross([a,0], [c,0]);

A=A(3);
B=B(3);
C=C(3);

if abs(A) < 1.e-14
    w1 = -C/B;
    w2 = w1;
else
    % Check for non complex solutions
    if B^2-4*A*C > 0
        w1 = (-B+sqrt(B^2-4*A*C))/(2*A);
        w2 = (-B-sqrt(B^2-4*A*C))/(2*A);
    else
        w1 = -1000;
        w2 = w1;
    end
end

lambda = -10000;
if(w1>=0 && w1<=1)
    lambda = w1;
end
if(w2>=0 && w2<=1)
    lambda = w2;
end

BC(1,1) = (1-mu)*(1-lambda);
BC(2,1) = lambda*(1-mu);
BC(3,1) = mu*lambda;
BC(4,1) = (1-lambda)*mu;
end


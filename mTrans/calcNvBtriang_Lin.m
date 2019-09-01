function NvB=calcNvBtriang_Lin(N,v,B,ii)
% NvB = calcNvBtriang_Lin(N, v, B, ii)
% Computes the product N x v x B for 2D triangular Linear elements in vectorized manner
% This is used internally.
% 
% Input 
% N     : shape functions
% v     : velocity
% B     : shape function derivatives
% ii    : in case of nested calculation this defines which loop will be
%        calculated. 
%
% Output
% NvB   : the product N x v x B computed in vectorized manner
%
% see also calcNvB
%
% Version : 1.1
% Author : George Kourakos
% email: gkourakos@ucdavis.edu
% web : https://gwt.ucdavis.edu/research-tools-and-applications/msim
% Date 1-Sep-2019
% Department of Land Air and Water
% University of California Davis
%
% How to compute the products using symbolic toolbox
% syms n1 n2 n3
% syms b1 b2 b3 b4 b5 b6
% syms vx vy
% N=[n1 n1;n2 n2; n3 n3];
% V=[vx 0; 0 vy];
% B=[b1 b2 b3;b4 b5 b6];
% NvB = N*V*B
%
% cnt=0;
% for i=1:3
%     for j=1:3
%         cnt=cnt+1;
%         fprintf(['NvB(:,%d) = ' char(NvB(i,j)) ';\n'],cnt);
%     end
% end

% see also calcNvB, Assemble_LHS_std

vx = v(:,1); vy = v(:,2);
if isempty(ii)
    NvB(:,1) = B(:,1).*N(1).*vx + B(:,4).*N(1).*vy;
    NvB(:,2) = B(:,2).*N(1).*vx + B(:,5).*N(1).*vy;
    NvB(:,3) = B(:,3).*N(1).*vx + B(:,6).*N(1).*vy;
    NvB(:,4) = B(:,1).*N(2).*vx + B(:,4).*N(2).*vy;
    NvB(:,5) = B(:,2).*N(2).*vx + B(:,5).*N(2).*vy;
    NvB(:,6) = B(:,3).*N(2).*vx + B(:,6).*N(2).*vy;
    NvB(:,7) = B(:,1).*N(3).*vx + B(:,4).*N(3).*vy;
    NvB(:,8) = B(:,2).*N(3).*vx + B(:,5).*N(3).*vy;
    NvB(:,9) = B(:,3).*N(3).*vx + B(:,6).*N(3).*vy;
else
    if ii == 1; NvB(:,1) = B(:,1).*N(1).*vx + B(:,4).*N(1).*vy; end
    if ii == 2; NvB(:,2) = B(:,2).*N(1).*vx + B(:,5).*N(1).*vy; end
    if ii == 3; NvB(:,3) = B(:,3).*N(1).*vx + B(:,6).*N(1).*vy; end
    if ii == 4; NvB(:,4) = B(:,1).*N(2).*vx + B(:,4).*N(2).*vy; end
    if ii == 5; NvB(:,5) = B(:,2).*N(2).*vx + B(:,5).*N(2).*vy; end
    if ii == 6; NvB(:,6) = B(:,3).*N(2).*vx + B(:,6).*N(2).*vy; end
    if ii == 7; NvB(:,7) = B(:,1).*N(3).*vx + B(:,4).*N(3).*vy; end
    if ii == 8; NvB(:,8) = B(:,2).*N(3).*vx + B(:,5).*N(3).*vy; end
    if ii == 9; NvB(:,9) = B(:,3).*N(3).*vx + B(:,6).*N(3).*vy; end
end


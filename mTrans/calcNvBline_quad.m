function NvB=calcNvBline_quad(N,v,B,ii)
% NvB = calcNvBline_quad(N, v, B, ii)
% Computes the product N x v x B for 1D quadratic elements in vectorized manner
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
% syms n1 n2 n3 b1 b2 b3 vx
% N=[n1;n2;n3];
% V=vx;
% B=[b1 b2 b3];
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

vx=v;
if isempty(ii)
    NvB(:,1) = B(:,1).*N(1)*vx;
    NvB(:,2) = B(:,2).*N(1)*vx;
    NvB(:,3) = B(:,3).*N(1)*vx;
    NvB(:,4) = B(:,1).*N(2)*vx;
    NvB(:,5) = B(:,2).*N(2)*vx;
    NvB(:,6) = B(:,3).*N(2)*vx;
    NvB(:,7) = B(:,1).*N(3)*vx;
    NvB(:,8) = B(:,2).*N(3)*vx;
    NvB(:,9) = B(:,3).*N(3)*vx;
else
    if ii == 1; NvB(:,1) = B(:,1).*N(1)*vx; end
    if ii == 2; NvB(:,2) = B(:,2).*N(1)*vx; end
    if ii == 3; NvB(:,3) = B(:,3).*N(1)*vx; end
    if ii == 4; NvB(:,4) = B(:,1).*N(2)*vx; end
    if ii == 5; NvB(:,5) = B(:,2).*N(2)*vx; end
    if ii == 6; NvB(:,6) = B(:,3).*N(2)*vx; end
    if ii == 7; NvB(:,7) = B(:,1).*N(3)*vx; end
    if ii == 8; NvB(:,8) = B(:,2).*N(3)*vx; end
    if ii == 9; NvB(:,9) = B(:,3).*N(3)*vx; end
end


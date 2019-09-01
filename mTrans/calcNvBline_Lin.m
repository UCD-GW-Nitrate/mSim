function NvB=calcNvBline_Lin(N,v,B,ii)
% NvB = calcNDBline_Lin(N, v, B, ii)
% Computes the product N x v x B for 1D linear elements in vectorized manner
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
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : http://groundwater.ucdavis.edu/msim
% Date 28-Mar-2014
% Department of Land Air and Water
% University of California Davis
%

vx=v;
if isempty(ii)
    NvB(:,1)= B(:,1).*N(1).*vx;
    NvB(:,2)= B(:,2).*N(1).*vx;
    NvB(:,3)= B(:,1).*N(2).*vx;
    NvB(:,4)= B(:,2).*N(2).*vx;
else
    if ii==1; NvB(:,1)= B(:,1).*N(1).*vx; end
    if ii==2; NvB(:,2)= B(:,2).*N(1).*vx; end
    if ii==3; NvB(:,3)= B(:,1).*N(2).*vx; end
    if ii==4; NvB(:,4)= B(:,2).*N(2).*vx; end
end

% syms n1 n2 b1 b2 vx
% N=[n1;n2];
% V=vx;
% B=[b1 b2];
% N*V*B
function BKB=calcBKBtriang_Lin(B,K,ii)
% BKB = calcBKBtriang_Lin(B, K, ii)
%
% Computes the product B'*K*B in vectorized manner for linear triangle
% elements
%
% Input
% B : [Nel x 6] contains the contributions of each element to the final
%      conductance matrix
% K : [Nel x Nanis] Hydraulic conductiviy element values. The number of columns 
%     is defined by the anisotropy. Maximum number is 3.
% ii : In case of nested assembly this indicates the iteration. In
%       vectorized assembly this is empty 
%
% Output
% BKB the product B'*K*B
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : http://groundwater.ucdavis.edu/msim
% Date 18-Dec-2012
% Department of Land Air and Water
% University of California Davis
%
% how to compute the product using matlab symbolic toolbox
% syms b1 b2 b3
% syms b4 b5 b6
% syms kx ky
% B=[b1 b2 b3; b4 b5 b6];
% BT = [b1 b4;b2 b5; b3 b6];
% BKB = BT*[kx kxy; kyx ky]*B
%
%
% cnt=0;
% for i=1:3
%     for j=1:3
%         cnt=cnt+1;
%         fprintf(['BKB(:,%d) = ' char(BKB(i,j)) ';\n'],cnt);
%     end
% end
%
% see also calcBKB, Assemble_LHS

Nel = size(B,1);
if size(K,2) == 1
    kx = K(:,1);
    ky = K(:,1);
    kxy = zeros(Nel,1);
    kyx = zeros(Nel,1);
elseif size(K,2) == 2
    kx = K(:,1);
    ky = K(:,2);
    kxy = zeros(Nel,1);
    kyx = zeros(Nel,1);
elseif size(K,2) == 4
    kx = K(:,1);
    kxy = K(:,2);
    kyx = K(:,3);
    ky = K(:,4);
end

if isempty(ii)
	BKB(:,1) = B(:,1).*(B(:,1).*kx + B(:,4).*kyx) + B(:,4).*(B(:,1).*kxy + B(:,4).*ky);
    BKB(:,2) = B(:,2).*(B(:,1).*kx + B(:,4).*kyx) + B(:,5).*(B(:,1).*kxy + B(:,4).*ky);
    BKB(:,3) = B(:,3).*(B(:,1).*kx + B(:,4).*kyx) + B(:,6).*(B(:,1).*kxy + B(:,4).*ky);
    BKB(:,4) = B(:,1).*(B(:,2).*kx + B(:,5).*kyx) + B(:,4).*(B(:,2).*kxy + B(:,5).*ky);
    BKB(:,5) = B(:,2).*(B(:,2).*kx + B(:,5).*kyx) + B(:,5).*(B(:,2).*kxy + B(:,5).*ky);
    BKB(:,6) = B(:,3).*(B(:,2).*kx + B(:,5).*kyx) + B(:,6).*(B(:,2).*kxy + B(:,5).*ky);
    BKB(:,7) = B(:,1).*(B(:,3).*kx + B(:,6).*kyx) + B(:,4).*(B(:,3).*kxy + B(:,6).*ky);
    BKB(:,8) = B(:,2).*(B(:,3).*kx + B(:,6).*kyx) + B(:,5).*(B(:,3).*kxy + B(:,6).*ky);
    BKB(:,9) = B(:,3).*(B(:,3).*kx + B(:,6).*kyx) + B(:,6).*(B(:,3).*kxy + B(:,6).*ky);
else
	if ii == 1; BKB(:,1) = B(:,1).*(B(:,1).*kx + B(:,4).*kyx) + B(:,4).*(B(:,1).*kxy + B(:,4).*ky); return;end
    if ii == 2; BKB(:,2) = B(:,2).*(B(:,1).*kx + B(:,4).*kyx) + B(:,5).*(B(:,1).*kxy + B(:,4).*ky); return;end
    if ii == 3; BKB(:,3) = B(:,3).*(B(:,1).*kx + B(:,4).*kyx) + B(:,6).*(B(:,1).*kxy + B(:,4).*ky); return;end
    if ii == 4; BKB(:,4) = B(:,1).*(B(:,2).*kx + B(:,5).*kyx) + B(:,4).*(B(:,2).*kxy + B(:,5).*ky); return;end
    if ii == 5; BKB(:,5) = B(:,2).*(B(:,2).*kx + B(:,5).*kyx) + B(:,5).*(B(:,2).*kxy + B(:,5).*ky); return;end
    if ii == 6; BKB(:,6) = B(:,3).*(B(:,2).*kx + B(:,5).*kyx) + B(:,6).*(B(:,2).*kxy + B(:,5).*ky); return;end
    if ii == 7; BKB(:,7) = B(:,1).*(B(:,3).*kx + B(:,6).*kyx) + B(:,4).*(B(:,3).*kxy + B(:,6).*ky); return;end
    if ii == 8; BKB(:,8) = B(:,2).*(B(:,3).*kx + B(:,6).*kyx) + B(:,5).*(B(:,3).*kxy + B(:,6).*ky); return;end
    if ii == 9; BKB(:,9) = B(:,3).*(B(:,3).*kx + B(:,6).*kyx) + B(:,6).*(B(:,3).*kxy + B(:,6).*ky); return;end
end
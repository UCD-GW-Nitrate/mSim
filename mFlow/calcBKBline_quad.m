function BKB=calcBKBline_quad(B,K,ii)
% BKB = calcBKBline_quad(B, K, ii)
% Computes the product B'*K*B in vectorized manner for quadratic line
% elements
%
% Input
% B : [Nel x 3] contains the contributions of each element to the final
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
% Date : 18-Dec-2012
% Department of Land Air and Water
% University of California Davis
%
% how to compute the product using matlab symbolic toolbox
% syms b1 b2 b3
% syms kx
% B=[b1 b2 b1]
% BT = [b1;b2;b2]
% BKB = BT*kx*B
%
% see also calcBKB, Assemble_LHS

k=K(:,1);
if isempty(ii)
	BKB(:,1)=B(:,1).^2.*k;
	BKB(:,2)=B(:,1).*B(:,2).*k;
	BKB(:,3)=B(:,1).*B(:,3).*k;
	BKB(:,4)=B(:,1).*B(:,2).*k;
	BKB(:,5)=B(:,2).^2.*k;
	BKB(:,6)=B(:,2).*B(:,3).*k;
	BKB(:,7)=B(:,1).*B(:,3).*k;
	BKB(:,8)=B(:,2).*B(:,3).*k;
	BKB(:,9)=B(:,3).^2.*k;
else
	if ii==1;BKB=B(:,1).^2.*k;return;end
	if ii==2;BKB=B(:,1).*B(:,2).*k;return;end
	if ii==3;BKB=B(:,1).*B(:,3).*k;return;end
	if ii==4;BKB=B(:,1).*B(:,2).*k;return;end
	if ii==5;BKB=B(:,2).^2.*k;return;end
	if ii==6;BKB=B(:,2).*B(:,3).*k;return;end
	if ii==7;BKB=B(:,1).*B(:,3).*k;return;end
	if ii==8;BKB=B(:,2).*B(:,3).*k;return;end
	if ii==9;BKB=B(:,3).^2.*k;return;end
end
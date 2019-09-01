function NLN = calcNLNtriang_Lin(N, L)


% How to compute the products using symbolic toolbox
% syms n1 n2 n3
% syms L
% N1 = [n1; n2; n3];
% N2 = transpose(N1);
% NLN = N1*L*N2
%
% cnt=0;
% for i=1:3
%     for j=1:3
%         cnt=cnt+1;
%         fprintf(['NLN(:,%d) = ' char(NLN(i,j)) ';\n'],cnt);
%     end
% end

NLN(:,1) = L.*N(:,1).^2;
NLN(:,2) = L.*N(:,1).*N(:,2);
NLN(:,3) = L.*N(:,1).*N(:,3);
NLN(:,4) = L.*N(:,1).*N(:,2);
NLN(:,5) = L.*N(:,2).^2;
NLN(:,6) = L.*N(:,2).*N(:,3);
NLN(:,7) = L.*N(:,1).*N(:,3);
NLN(:,8) = L.*N(:,2).*N(:,3);
NLN(:,9) = L.*N(:,3).^2;

% see also calcNvB, Assemble_LHS_std
function [B Jdet]=shapeDerivLine_Lin(p,MSH,~)
% [B Jdet] = shapeDerivLine_Lin(p, MSH, ~)
%
% Computes the shape function derivatives and the determinant for
% of 1D linear isoparametric elements.
% This function is used internally from shapeDerivatives.
% 
% Input
% p   : [Np x 3] Coodrinates of nodes [x1 y1 z1; x2 y2 z2;...xn yn zn]
%        Np : Number of nodes
%
% MSH : [Nel x Np_el] id of elements. Each row correspond to an element
%        Nel : Number of elements
%        Np_el : Number of nodes to define the element
% Output
% B    : Shape function derivatives
% Jdet : The determinant of the Jacobian Matrix 
%
%
% Shape functions of linear line elements
% N1=0.5*(1-ksi)
% N2=0.5*(1+ksi)
%
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : http://groundwater.ucdavis.edu/msim
% Date 18-Dec_2012
% Department of Land Air and Water
% University of California Davis
%
% see also shapeDerivatives

%linearize the coordinates
[b n m]=unique([MSH(:,1);MSH(:,2)]);
p_tmp=p(b,:);
p_tmp_min=min(p_tmp,[],1);
p_tmp_max=max(p_tmp,[],1);
Len=sum(sqrt((p_tmp_max-p_tmp_min).^2));
% Since we use this function often in octave where nanmin is missing we
% will add the following hack
if exist ("OCTAVE_VERSION", "builtin")
    % Many years later I can't remember why I had to use nanmean instead of
    % mean. It seems impossible to have nan values here, unless the domain
    % is zero. (25-Sep-2019)
    p_param = mean(bsxfun(@rdivide,bsxfun(@minus,p_tmp,p_tmp_min),(p_tmp_max-p_tmp_min)),2);
else
    p_param = nanmean(bsxfun(@rdivide,bsxfun(@minus,p_tmp,p_tmp_min),(p_tmp_max-p_tmp_min)),2);
end

p_real=p_param.*Len;
p_real=reshape(p_real(m),size(MSH,1),2);

x1=p_real(:,1); 
x2=p_real(:,2);

test=x2<x1;
x2(test)=p_real(test,1);
x1(test)=p_real(test,2);


%Derivatives of shape functions wrt the parametric variable ksi
gn1=-1/2;
gn2= 1/2;

Jdet=gn1.*x1+gn2.*x2;

B=[gn1./Jdet gn2./Jdet];
%Jdet=x2-1;%in 1D Jdet simplifies to element lenght
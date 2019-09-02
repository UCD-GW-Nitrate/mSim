function [B, Jdet]=shapeDerivatives(p, MSH, n, opt)
% [B Jdet]=shapeDerivatives(p, MSH, n, opt)
%
% Computes the shape function derivatives and the determinant of the
% Jacobian matrix.
%
% Input
% p   : [Np x 3] Coodrinates of nodes [x1 y1 z1; x2 y2 z2;...xn yn zn]
%        Np : Number of nodes
%
% MSH : [Nel x Np_el] id of elements. Each row correspond to an element
%        Nel : Number of elements
%        Np_el : Number of nodes to define the element
% n   : the intergration point where the derivatives will be evaluated
% opt : Structure with the following fields
%      dim : Dimension of the problem : 1, 2 or 3.
%      el_type : type of element e.g. triangle, prism etc... 
%      el_order order of the element: linear or quadratic 
% 
% Output
% B    : Shape function derivatives
% Jdet : The determinant of the Jacobian Matrix 
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : https://gwt.ucdavis.edu/research-tools-and-applications/msim
% Date : 2-Sep-2019
% Department of Land Air and Water
% University of California Davis
%


if opt.dim==1
    switch opt.el_order
        case 'linear'
            [B Jdet]=shapeDerivLine_Lin(p,MSH,n);
        case 'quadratic'
            [B Jdet]=shapeDerivLine_quad(p,MSH,n);
    end
elseif opt.dim==2
    switch opt.el_type
        case 'triangle'
            switch opt.el_order
                case 'linear'
                    if ~isfield(opt,'project')
                        opt.project=0;
                    end
                    [B Jdet]=shapeDerivTriang_Lin(p,MSH,n,opt.project);
                case 'quadratic'
                    if ~isfield(opt,'project')
                        opt.project=0;
                    end
                    [B Jdet]=shapeDerivTriang_quad(p,MSH,n,opt.project);
            end
        case 'quad'
            switch opt.el_order
                case 'linear'
                    if ~isfield(opt,'project')
                        opt.project=0;
                    end
                    [B Jdet]=shapeDerivQuad_Lin(p,MSH,n,opt.project);
                case 'quadratic_9'
                    if ~isfield(opt,'project')
                        opt.project=0;
                    end
                    [B Jdet]=shapeDerivQuad_quad_9(p,MSH,n,opt.project);
            end   
    end
elseif opt.dim==3
    switch opt.el_type
        case 'prism'
            switch opt.el_order
                case 'linear'
                    [B Jdet]=shapeDerivPrism_Lin(p,MSH,n);
                case 'quadratic'
                    [B Jdet]=shapeDerivPrism_quad(p,MSH,n);
            end
        case 'hex'
            switch opt.el_order
                case 'linear'
                    [B Jdet]=shapeDerivHex_Lin(p,MSH,n);
                case 'quadratic_27'
                    [B Jdet]=shapeDerivHex_quad_27(p,MSH,n);
            end
            
    end
end

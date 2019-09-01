function NvB=calcNvB(N,v,B,opt,ii)
% NvB = calcNvB(N, v, B, opt, ii)
% 
% Computes the product N x v x B in vectorized manner, where N are the
% shape function values, v is the velocity and B the shape function
% derivatives. This function is used internally. Based on the options
% selects the appropriate function to perform the vectorized computations
%
% Input:
% N     : shape functions
% v     : velocity
% B     : shape function derivatives
% opt   : option structure. for details see Assemble_LHS_std
% ii    : in case of nested calculation this defines which loop well be
%        calculated. (Currently this is not used)
%
% Output
% NvB   : the product N x v x B computed in vectorized manner
%
% see also Assemble_LHS_std, calcNDBline_Lin
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : http://groundwater.ucdavis.edu/msim
% Date 28-Mar-2014
% Department of Land Air and Water
% University of California Davis
%

if opt.dim == 1
    switch opt.el_order
        case 'linear'
            NvB = calcNvBline_Lin(N,v,B,ii);
        case 'quadratic'
            NvB = calcNvBline_quad(N,v,B,ii);
    end
elseif opt.dim == 2
    switch opt.el_type
        case 'triangle'
            switch opt.el_order
                case 'linear'
                    NvB = calcNvBtriang_Lin(N,v,B,ii);
                case 'quadratic'
                    NvB = calcNvBtriang_quad(N,v,B,ii);
            end
        case 'quad'
            switch opt.el_order
                case 'linear'
                    NvB = calcNvBQuad_Lin(N,v,B,ii);
                case 'quadratic_9'
                    NvB = calcNvBQuad_quad(N,v,B,ii);
            end
    end
end
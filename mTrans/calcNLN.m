function NLN=calcNLN(N, L, opt)
% NLN = calcNLN(N, L)
% Computes the product N' x L x N for 1D linear elements in vectorized manner
% This is used internally. Currently the computations are performed within
% this function. However in the future this will be used as wraper wich
% will call the appropriate function depending the element type and order.
% 
% Documentation to be completed

old = false;

if old
    Nsh=length(N);
    Nel=size(L,1);
    NLN=zeros(Nel,Nsh^2);
    for ii=1:Nsh
        for jj=1:Nsh
            i_lin=sub2ind([Nsh Nsh],jj,ii);
            NLN(:,i_lin)=N(ii).*L.*N(jj);
        end
    end
else
    if opt.dim == 1
        switch opt.el_order
            case 'linear'
                Nsh=length(N);
                Nel=size(L,1);
                NLN=zeros(Nel,Nsh^2);
                for ii=1:Nsh
                    for jj=1:Nsh
                        i_lin=sub2ind([Nsh Nsh],jj,ii);
                        NLN(:,i_lin)=N(ii).*L.*N(jj);
                    end
                end
            case 'quadratic'
                error('calcNLN has not been implemeted for 1D quadratic elements'); 
        end
    elseif opt.dim == 2
        switch opt.el_type
            case 'triangle'
                switch opt.el_order
                    case 'linear'
                        NLN = calcNLNtriang_Lin(N,L);
                    case 'quadratic'
                        error('calcNLN has not been implemeted for 2D triangle quadratic elements');
                end
            case 'quad'
                switch opt.el_order
                    case 'linear'
                        error('calcNLN has not been implemeted for 2D quadrilateral linear elements');
                    case 'quadratic_9'
                        error('calcNLN has not been implemeted for 2D quadrilateral quadratic elements');
                end
        end
    end
end

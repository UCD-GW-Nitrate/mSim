function Vel = calcVelocityField(p, MSH, H, K, P, opt)

% Calculate the shape function derivatives for all elements at the center 0.5 0.5

switch opt.el_type
    case 'line'
        warning(' Center point for LINE elements Not tested yet');
        n = [0.5];
    case 'triangle'
        n = [0.5 0.5];
    case 'quad'
        warning(' Center point for QUAD elements Not tested yet');
        n = [0 0];
    case 'prism'
        warning(' Center point for PRISM elements Not tested yet');
        n = [0.5 0.5 0];
    case 'hex'
        warning(' Center point for HEX elements Not tested yet');
        n = [0 0 0];
end

[B, ~] = shapeDerivatives(p, MSH, n, opt);

dH = zeros(size(MSH,1),opt.dim);

ndof = size(MSH,2);
for ii = 1:ndof
    dH(:,1) = dH(:,1) + B(:,ii).*H(MSH(:,ii),1);
    if opt.dim > 1
        dH(:,2) = dH(:,2) + B(:,ii + ndof).*H(MSH(:,ii),1);
    end
    if opt.dim > 2
        dH(:,3) = dH(:,3) + B(:,ii + 2*ndof).*H(MSH(:,ii),1);
    end
end
Vel = -bsxfun(@times, bsxfun(@times,K,dH), P);
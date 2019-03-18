function [Xs, Ys] = msim_polysplit(X, Y)
% [Xs, Ys] = msim_polysplit(X, Y)
% This is equivalent to the matlab polysplit command. The only reason that
% we reimplement this, it's to avoid the dependency on the mapping
% toolbox. In addition as for now the Octave is missing this function
% therefore many functions of mSim couldn't run.
%
% Input:
% X, Y : Vectors of X and Y coordinates. These can be the X and Y fields of
%        the readshapefile output for example.
%
% Output:
% Xs, Ys: a cell array with the parts of the X and Y inputs split at nan
%        locations
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : http://subsurface.gr/software/msim/
% Date 18-Mar-2019
% Department of Land Air and Water
% University of California Davis

X = X(:);
Y = Y(:);
% find nans
i_nanX = find(isnan(X));
i_nanY = find(isnan(Y));
if sum(i_nanX - i_nanY) ~= 0
    error('The nan in the X and Y must be identical');
end

if isempty(i_nanX)
   Xs{1,1} = X;
   Ys{1,1} = Y;
   return;
end

istart = 1;
cnt = 1;
for ii = 1:length(i_nanX)
    iend = i_nanX(ii)-1;
    
    if iend - istart >= 0
        Xs{cnt,1} = X(istart:iend);
        Ys{cnt,1} = Y(istart:iend);
        cnt = cnt + 1;
    end
    istart = i_nanX(ii)+1;
end



function S = mesh2poly( p, msh, props, nameprops)
% S = mesh2poly( p, msh)
% or
% S = mesh2poly( p, msh, props, nameprops)
%
% Converts an element mesh to a shapefile structure variable.  
%
%   Input
% p   : [Np x 2] a matrix that holds the coordinates of the element nodes.
%
% msh : [Nel x 3 or 4] contains the connectivity ids.
%
% props: [Nel x Nprop] (optional) contains property values for
%        each mesh element
%
% nameprops: (optional) is a cell variable with the names of the
%            properties. The lenght of this variable is Nprop.
%
%
% Output 
% S : A shapefile structure variable
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : https://gwt.ucdavis.edu/research-tools-and-applications/msim
% Date 09-Sep-2019 
% Department of Land Air and Water
% University of California Davis

Nel = size(msh,1);

S(Nel,1).Geometry = 'Polygon';
S(Nel,1).X = [];
S(Nel,1).Y = [];
if nargin == 4
    Nprop = length(nameprops);
    if Nprop ~= size(props,2)
        error('The size of the properties must be equal to the property names');
    end
    for ii = 1:length(nameprops)
        S(Nel,1).(nameprops{ii}) = [];
    end
end

for ii = 1:Nel
    S(ii,1).Geometry = 'Polygon';
    S(ii,1).X = [p([msh(ii,:) msh(ii,1)],1)' nan];
    S(ii,1).Y = [p([msh(ii,:) msh(ii,1)],2)' nan];
    if nargin == 4
        for jj = 1:Nprop
            S(ii,1).(nameprops{jj}) = props(ii,jj);
        end
    end
end

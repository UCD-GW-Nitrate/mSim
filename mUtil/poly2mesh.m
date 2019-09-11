function [p, msh, prop] = poly2mesh( S, type, thres )
% [p, MSH] = poly2mesh( S )
%
% Converts a shapefile like structure of polygons to a mesh. As the
% purpose of this function is to convert a shapefile into a mesh that will
% be used eventually for finite element computations, it make sense that the
% polygon shapes be of the same type e.g triangles or quadrilaterals and 
% of course they should not contain islands or be multipart polygons. 
% The code will still work but the results will be interesting at least. 
%
%   Input
% S   : is a structure variable with the following fields
%     Geometry : This should always be 'Polygon'. 
%     X : The X coordinates of the polygon followed by nan at the end 
%           i.e. [X1, X2, X3, X1, nan]
%     Y : The Y coordinates of the polygon followed by nan at the end 
%           i.e. [Y1, Y2, Y3, Y1, nan]
%
% This structure assumes that the first and last point of each polygon are
% identical
%
% type : this is either triangle or quad.
%
% thres : is the search threshold. If the distance between two points is
% smaller than threshold then the two points are considered identical.
%
%
% Output 
% p    : [Np x 3] The coordinate of the nodes
%
% msh  : [Nel x 3 or 4] are the ids of the polygon nodes in the p variable
%
% prop : [Nelx1] is a structure that contains the remaining properties of 
%        the shapefile
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : https://gwt.ucdavis.edu/research-tools-and-applications/msim
% Date 09-Sep-2019 
% Department of Land Air and Water
% University of California Davis

    if nargin < 3
        thres = 0.001;
    end

    Nel = length(S);
    switch type
        case 'triangle'
            msh = zeros(Nel,3);
        case 'quad'
            msh = zeros(Nel,4);
    end
    p = nan(2*Nel,2); % allocate more space than what we need
    cnt = 0;

    input_fields =  fields(S);
    % delete the Geometry, X and Y fields
    dlt = [];
    for ii = 1:length(input_fields)
        if strcmp(input_fields{ii},'Geometry') || strcmp(input_fields{ii},'X') ||...
                strcmp(input_fields{ii},'Y')
            dlt = [dlt ii];
        end
    end
    input_fields(dlt,:) = [];
    if isempty(input_fields)
        prop = [];
    else
        for ii = 1:length(input_fields)
            prop(Nel,1).(input_fields{ii}) = [];
        end
    end


    for ii = 1:Nel
        [Xs, Ys] = msim_polysplit(S(ii,1).X, S(ii,1).Y);
        % This assumes
        for jj = 1:length(Xs{1})-1
            id = point_exist(Xs{1,1}(jj), Ys{1,1}(jj), p, cnt, thres);
            if id < 0
                cnt = cnt + 1;
                p(cnt,:) = [Xs{1,1}(jj) Ys{1,1}(jj)];
                msh(ii,jj) = cnt;
            else
                msh(ii,jj) = id;
            end
        end
        for jj = 1:length(input_fields)
            prop(ii,1).(input_fields{jj}) = S(ii,1).(input_fields{jj});
        end
    end
    % delete the extra points
    p(cnt+1:end,:) = [];
end

% This function searches the first cnt points in p to see if the point x,y
% already exists and returns the id. If the point is not in p returns -9
function id = point_exist(x, y, p, cnt, thres)
    id = -9;
    if cnt > 0
        [c, d] = min(sqrt((p(1:cnt,1) - x).^2 + (p(1:cnt,2) - y).^2)) ;
        if c < thres
            id = d;
        end
    end
end



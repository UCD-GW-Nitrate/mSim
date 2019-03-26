function y = dewater_obj_fun(x, L, U, P, R, Q, GGDD, H, F, id_var, p, trimesh, cc, id_cnstr)
% x is a vector 6*3 of the decision variables.
% we reshape it so that they correspond to X Y Q
x = reshape(x, 6,3);
% distance from center
dst = sqrt((x(:,1) - 2500).^2 + (x(:,2) - 2500).^2);
% find if the wells violate the location contraints
id_in = dst < 650;
id_out = dst > 2000;

if any([id_in id_out])
    y = 100000*(sum(650-dst(id_in)) + sum(dst(id_out) - 2000)); 
else
    FLUX_point = [];
    weights = zeros(size(trimesh,2), 1);
    for ii = 1:size(x, 1)
        dst = sqrt((x(ii,1) - cc(:,1)).^2 + (x(ii,2) - cc(:,2)).^2);
        [~, d] = sort(dst);
        for jj = 1:length(d)
            in = inpolygon(x(ii,1), x(ii,2), p(trimesh(d(jj),:),1), p(trimesh(d(jj),:),2));
            if in
                node_distances = sqrt((x(ii,1) - p(trimesh(d(jj),:),1)).^2 + (x(ii,2) - p(trimesh(d(jj),:),2)).^2);
            if any(node_distances < 0.001)
                weights(:,1) = 0;
                weights(node_distances < 0.001,:) = 1;
            else
                weights = (1./node_distances)./sum(1./node_distances);
            end
            FLUX_point = [FLUX_point; trimesh(d(jj),:)', x(ii,3)* weights];
            break 
            end
        end
    end
    Fall = F + sparse(FLUX_point(:,1),1,FLUX_point(:,2), length(H), 1);
    B = Fall(id_var)-GGDD;
    X = Q * (U \ (L \ (P * (R \ B))));
    H(id_var,1) = X;
    H_constraint = H(id_cnstr);
    id_violate =  H_constraint > 30;
    if any(id_violate) 
        y = 100*sum((id_violate)) + abs(sum(x(:,3)));
    else
       y = abs(sum(x(:,3))); 
    end
end


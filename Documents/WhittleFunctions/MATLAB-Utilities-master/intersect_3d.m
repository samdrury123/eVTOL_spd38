function p = intersect_3d(v1,v2,n,c)
% INTERSECT_3D  Calculate intersection coordinates between a plane and a line

% Measure perpendicular distances between plane and points
d1 = dot(v1,n) - c; d2 = dot(v2,n) - c;

% Calculate coordinates of intersection
p = v1 + (v2 - v1) * d1 / (d1 - d2);


end


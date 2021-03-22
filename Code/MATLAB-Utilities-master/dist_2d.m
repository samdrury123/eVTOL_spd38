function s = dist_2d(xy,normalise)
% DIST_2D  Calculate cumilative distance through a polyline

% Default to not normalise
if exist('normalise','var') == 0
    normalise = 0;
end

% Difference in distances in final dimension
s = [0 ; cumsum(sum(diff(xy,1,1).^2,2).^0.5,1)];

% Normalise to 1 if required
if normalise == 1
    s = s/s(end);
end

end
function n = norm_2d(xy,is_loop)
% NORM_2D  Calculate the 2D normal vectors to a polyline

% Default to not a loop
if exist('is_loop','var') == 0
    is_loop = 0;
end

% Add points to the line if a loop
if is_loop == 1
    xy = xy([end-4:end-1 1:end 2:5],:);
end

% Get normal direction from differences between points
n = [diff(xy(:,2)) -diff(xy(:,1))];

% Check the size of the normal list
if size(n,1) > 1
    
    % Take midpoint averages of normals and add on end points
    n = 0.5 * (n(2:end,:) + n(1:end-1,:)); 
    n = [n(1,:) ; n ; n(end,:)]; 
else
    
    % Take the single normal for both points
    n = repmat(n,[2 1]);
    
end

% Normalise to unit length
n = n ./ repmat(sum(n.^2,2).^0.5,[1 2]);

% Remove added points
if is_loop == 1
    n = n(5:end-4,:);
end


end
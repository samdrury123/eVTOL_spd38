function x = smooth_2d(x,is_loop,renormalise,span,method)
% SMOOTH_2D  Smooth a vector of coordinates

% Default to not a loop
if exist('is_loop','var') == 0
    is_loop = 0;
end

% Default to not renormalise
if exist('renormalise','var') == 0
    renormalise = 0;
end

% Default to 5 point span
if exist('span','var') == 0
    span = 5;
end

% Default to moving average method 
if exist('method','var') == 0
    method = 'moving';
end

% Add points to the line if a loop
if is_loop == 1
    x = x([end-(span-1):end-1 1:end 2:span],:);
end

% Apply moving average
for i = 1:size(x,2); x(:,i) = smooth(x(:,i),span,method); end

% Normalise if requested
if renormalise == 1
    x = x ./ repmat(sum(x.^2,2).^0.5,[1 size(x,2)]);
end

% Remove added points
if is_loop == 1
    x = x(span:end-(span-1),:);
end

end
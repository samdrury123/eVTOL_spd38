function xy = generate_block(xy_1,xy_2,s,n_smooth,plot_stuff)
% GENERATE_BLOCK  Algegraic block generation with elliptic smoothing

% Default number of smoothing iterations
if exist('n_smooth','var') == 0
    n_smooth = 10;
end

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Grid size
ni = length(s); nj = size(xy_1,1);

% Linear grid interpolation 
xy = repmat(reshape(s,[ni 1 1]),[1 nj 2]) .* repmat(reshape(xy_2 - xy_1,[1 nj 2]),[ni 1 1]) + ...
    repmat(reshape(xy_1,[1 nj 2]),[ni 1 1]); 

% Plot default grid
if plot_stuff == 1
    figure(); hold on; grid on; box on; axis equal;
    plot(xy_1(:,1),xy_1(:,2),'.-')
    plot(xy_2(:,1),xy_2(:,2),'.-')
    mesh(xy(:,:,1),xy(:,:,2),zeros(ni,nj),zeros(ni,nj,3))
end

% Overwrite coordinates with normals to line
si = [zeros(1,nj) ; cumsum(sum(diff(xy,1,1).^2,3).^0.5,1)];
xy_1b = xy_1 + repmat(diff(si(1:2,:),1,1)',[1 2]) .* norm_2d(xy_1);
xy_2b = xy_2 - repmat(diff(si(end-1:end,:),1,1)',[1 2]) .* norm_2d(xy_2);

% Smooth normal offset lines
for n = 1:2
    xy_1b(:,n) = smooth(xy_1b(:,n),4);
    xy_2b(:,n) = smooth(xy_2b(:,n),4);
end


end
function xyz = rolling_sphere(xyz,rad,plot_stuff)
% ROLLING_SPHERE  Roll a ball on a surface and change coordinates to reduce curvature

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Size of the array
[ni,nj,~] = size(xyz); i_mid = round(ni/2); j_mid = round(nj/2);

% Plot the current array
if plot_stuff == 1
    figure(); hold on; grid on; box on; axis equal;
    mesh(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),zeros(ni,nj,3))
end

% Define quadrants of points to operate on
jl{1} = j_mid:nj-1; il{1} = i_mid:ni-1; 
jl{2} = j_mid-1:-1:2; il{2} = i_mid:ni-1; 
jl{3} = j_mid:nj-1; il{3} = i_mid-1:-1:2; 
jl{4} = j_mid-1:-1:2; il{4} = i_mid-1:-1:2; 

% Assemble quadrants into a single list
is = []; js = [];
for m = 1:length(jl)
    nil = length(il{m}); njl = length(jl{m});
    is = [is ; repmat(il{m}',[njl 1])]; js = [js ; reshape(repmat(jl{m},[nil 1]),[nil*njl 1])];
end

% Loop over all points in the surface
for m = 1:length(is)
    
    % Current point
    i = is(m); j = js(m);
    
    % Calculate normals near the current point
    a2 = xyz(i:i+1,j:j+1,:) - xyz(i-1:i,j-1:j,:);
    a1 = xyz(i-1:i,j:j+1,:) - xyz(i:i+1,j-1:j,:);
    n = cross(a1,a2); n = n ./ repmat(sum(n.^2,3).^0.5,[1 1 3]);
    
    % Average normals near the point
    n = 0.25 * (n(1,1,:) + n(2,1,:) + n(2,2,:) + n(1,2,:));
    n = n ./ sum(n.^2,3).^0.5;
    
    % Offset to sphere centre
    xyz_cen = xyz(i,m,:) + rad * n;
    
    % Check distance from all points to the current one
    d = reshape(pdist2(reshape(xyz,[ni*nj 3]),reshape(xyz_cen,[1 3])),[ni nj]);
    
    % Look for points which are inside the sphere
    q = d < rad;
    if sum(q(:)) > 0
        
        % Displace the points by the normal
        
        
        
    end
end




end
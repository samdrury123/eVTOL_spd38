function qi = interp_3D(X,Y,Z,Q,xi,yi,zi)
% Linearly interpolate functions from a general non-uniform 3D grid

% Pull out all cell corner coordinates along with functions
u_xyz = cat(3,corners(X),corners(Y),corners(Z));
u_q = corners(Q);

% Centre coordinates of cells
c(:,1,1) = sum(u_xyz(:,:,1),2) / 8;
c(:,1,2) = sum(u_xyz(:,:,2),2) / 8;
c(:,1,3) = sum(u_xyz(:,:,3),2) / 8;

% Calculate radius based on centre that encloses all nodes - pythagoras on
% differences between cell centre coord and cell node coords
c(:,1,4) = max( (u_xyz(:,:,1) - repmat(c(:,1,1),[1 8])).^2  + ...
    (u_xyz(:,:,2) - repmat(c(:,1,2),[1 8])).^2  + ...
    (u_xyz(:,:,3) - repmat(c(:,1,3),[1 8])).^2  ,[],2);

% All diagonal vectors
a = cat(2, u_xyz(:,1,:) - u_xyz(:,5,:),...
    u_xyz(:,1,:) - u_xyz(:,7,:),...
    u_xyz(:,3,:) - u_xyz(:,8,:),...
    u_xyz(:,1,:) - u_xyz(:,6,:),...
    u_xyz(:,2,:) - u_xyz(:,8,:),...
    u_xyz(:,4,:) - u_xyz(:,8,:));
b = cat(2, u_xyz(:,2,:) - u_xyz(:,3,:),...
    u_xyz(:,4,:) - u_xyz(:,3,:),...
    u_xyz(:,7,:) - u_xyz(:,5,:),...
    u_xyz(:,4,:) - u_xyz(:,2,:),...
    u_xyz(:,5,:) - u_xyz(:,6,:),...
    u_xyz(:,6,:) - u_xyz(:,7,:));

% Face centred normals, coordinates and function values
nq = size(Q,4); nu = size(u_xyz,1);
f = zeros(nu, 6, 6+nq);
f(:,:,4:6) = faces(u_xyz);

% Point normals to cell centre by dotting with cell centre to face centre
% vectors
f(:,:,1:3) = -cross(a,b) .* sign(repmat(dot(cross(a,b),...
    f(:,:,4:6) - repmat(c(:,:,1:3),[1 6 1]),3),[1 1 3]));

% Put in face centred function values
f(:,:,7:end) = faces(u_q);

clear u_q; clear u_xyz; clear a; clear b;

% Reshape desired coordinates into a list
size_i = size(xi); ni = numel(xi);
xi = reshape(xi,[ni 1]); yi = reshape(yi,[ni 1]); zi = reshape(zi,[ni 1]);

% Initialise matrices to store bounding cell data
f_in = nan(ni,6,size(f,3));

% Loop over all points and find bounding cell
for n = 1:ni
    % Determine which spheres point falls in
    p1 = ((xi(n) - c(:,:,1)).^2 + (yi(n) - c(:,:,2)).^2 + (zi(n) - c(:,:,3)).^2 ).^0.5 ...
        < c(:,:,4);

    % Resolve current point into coords of all selected face normals and
    % check for positive values
    dx = f(p1,:,4) - xi(n); dy = f(p1,:,5) - yi(n); dz = f(p1,:,6) - zi(n); 
    
    p2 = sum(dx.*f(p1,:,1) + dy.*f(p1,:,2) + dz.*f(p1,:,3) < 0,2) == 6;
    
    % Intepret logical array of cell centre method and spherical method 
    p1(p1 == 1) = p2;
    
    % Record "in" cell
    if sum(p1) == 1
        f_in(n,:,:) = f(p1,:,:);
    end
end

% Interpolate function properties in cell with barycentric coodinates weighted on
% volume of square based pyramid and face centred values of the functions
dx = repmat(xi,[1 6]) - f_in(:,:,4);
dr = repmat(yi,[1 6]) - f_in(:,:,5);
drt = repmat(zi,[1 6]) - f_in(:,:,6);

v = (dx.*f_in(:,:,1) + dr.*f_in(:,:,2) + drt.*f_in(:,:,3)) / 3;

% Calculate weights based on opposite pyramid
w = flipdim(v,2) ./ repmat(sum(v,2),[1 6]);

% Function values at new coordinates
qi = sum(repmat(w,[1 1 nq]) .* f_in(:,:,7:end),2); 
    
% Reshape back to original requested shape
qi = reshape(qi,[size_i nq]);

end

function x_faces = faces(x)
% Pull out face averages from corner list

x_faces = cat(2, x(:,2,:) + x(:,3,:) + x(:,1,:) + x(:,5,:),...
    x(:,4,:) + x(:,3,:) + x(:,1,:) + x(:,7,:),...
    x(:,7,:) + x(:,5,:) + x(:,3,:) + x(:,8,:),...
    x(:,4,:) + x(:,2,:) + x(:,1,:) + x(:,6,:),...
    x(:,5,:) + x(:,6,:) + x(:,2,:) + x(:,8,:),...
    x(:,6,:) + x(:,7,:) + x(:,4,:) + x(:,8,:)) / 4;

end

function x_corners = corners(x)
% Pull out corner coordinates and assemble into a list

% Reshape 4D inputs - such as interpolating functions
x_size = size(x); x_size = [x_size 1];
x = reshape(x,[x_size(1:3) 1 x_size(4)]);

% Assemble all corners along the 4th dimension
x_corners = cat(4,x(1:end-1,1:end-1,1:end-1,:,:),x(2:end,1:end-1,1:end-1,:,:),...
    x(1:end-1,2:end,1:end-1,:,:), x(1:end-1,1:end-1,2:end,:,:),...
    x(2:end,2:end,1:end-1,:,:), x(2:end,1:end-1,2:end,:,:),...
    x(1:end-1,2:end,2:end,:,:), x(2:end,2:end,2:end,:,:));

% Reshape corners into a list
x_corners = reshape(x_corners,[(x_size(1)-1)*(x_size(2)-1)*(x_size(3)-1) 8 x_size(4)]);

end
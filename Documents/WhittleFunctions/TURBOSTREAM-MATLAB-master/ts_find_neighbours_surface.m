function [v_neighbours] = ts_find_neighbours_surface(u, surface, surface_type, connectivity, coincident_flag)
% Search around a node and return neighbouring node information

% Find full I J K co-ordinates of surface of interest 
I = surface{u(1),u(2)}(u(3),u(4),1);
J = surface{u(1),u(2)}(u(3),u(4),2);
K = surface{u(1),u(2)}(u(3),u(4),3);

% Read connectivity list
C = connectivity{u(1)}{I+1,J+1,K+1};

% Read x, r, rt coords of interest
coords_u = surface{u(1), u(2)}(u(3), u(4),4:6);

% Determine coords of co-incident connected nodes
u_connected = zeros(size(C,1),4);
coincident = zeros(size(C,1),1);
for n = 1:size(C,1)
    nxSurf = surface{C(n,1)+1,C(n,2)};
    nxsurf_type = surface_type{C(n,1)+1,C(n,2)};

    nxIJK = [C(n,3) C(n,4) C(n,5)];
    xy = nxIJK(nxsurf_type(1:2));
    xy = xy+1; % Convert notation for matlab 1:nk+1 indexing
    % Adjust for surfaces where nStart ~= 0
    xy(1) = xy(1) - nxSurf(1,1,nxsurf_type(1));
    xy(2) = xy(2) - nxSurf(1,1,nxsurf_type(2));
    
    u_connected(n,:) = [C(n,1)+1 C(n,2) xy];
    
    coords_u_connected = nxSurf(xy(1),xy(2),4:6);

    coincident(n) = isequal(coords_u_connected,coords_u);
end

% Allow movement along diagonals
A = [ 0 0 1 0 ; 0 0 1 1 ; 0 0 0 1 ; 0 0 -1 1 ;...
    0 0 -1 0 ; 0 0 -1 -1 ; 0 0 0 -1 ; 0 0 1 -1 ];

A_u = repmat(A,6,1);
for n = 1:6
    A_u(8*n-7:8*n,2) = ones(8,1) * (n - u(2));
end

% Initialise neighbours array
if coincident_flag == 1
    coincident_extras = sum(coincident);
else
    coincident_extras = sum(coincident) * size(A,1);
end
v_neighbours = zeros(coincident_extras + sum(coincident==0) + size(A_u,1), 4);

% Put in points surrounding u, including moving the surfaces
n = size(A_u,1);
v_neighbours(1:n,:) = A_u + repmat(u,size(A_u,1),1);

% Put in non-coincident connected points
v_neighbours(n+1:n+sum(coincident==0),:) = u_connected(coincident==0,:);
n = n+sum(coincident==0);

% Put in coincident points
if coincident_flag == 1
    v_neighbours(n+1:n+sum(coincident),:) = u_connected(coincident == 1,:);
else
    % Put in points surrounding coincident ones
    coinc_index = find(coincident == 1);
    for m = 1:sum(coincident)
        v_neighbours(n+m*size(A,1)-(size(A,1)-1):n+m*size(A,1),:) = A + repmat(u_connected(coinc_index(m),:), size(A,1), 1);
    end
end

% Remove neighbours that don't exist in surface
for m = size(v_neighbours,1): -1 : 1
    v = v_neighbours(m,:);
    if (v(3)<=size(surface{v(1),v(2)},1) && v(3)>0 &&...
        v(4)<=size(surface{v(1),v(2)},2) && v(4)>0) == 0
        v_neighbours(m,:) = [];
    end
end

% Remove neighbours that have jumped too far accross surfaces
for m = size(v_neighbours,1): -1 : 1
    v = v_neighbours(m,:);
    n_jump = abs(surface{v(1),v(2)}(v(3),v(4),1) - surface{u(1),u(2)}(u(3),u(4),1))...
        + abs(surface{v(1),v(2)}(v(3),v(4),2) - surface{u(1),u(2)}(u(3),u(4),2))...
        + abs(surface{v(1),v(2)}(v(3),v(4),3) - surface{u(1),u(2)}(u(3),u(4),3));
    if v(1) == u(1) && n_jump > 2
        v_neighbours(m,:) = [];
    end
end

end
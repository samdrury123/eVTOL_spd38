function [v_neighbours] = ts_find_neighbours(u, coords, connectivity,...
    mwall, surface, surface_type, diag_switch, surface_switch, coincident_switch, return_surface, input_surface)
% Search around a node and return adjacent nodes.
% Input either [bid+1 I J K] or [bid+1 s I J K]

% If diag_switch = 1 then allow movement on diagonals
% If surface_switch = 1 then only return surface nodes
% If coincident_switch = 1 then only return coincident points rather than
% those around it
% If return_surface = 1 then format of v is [bid+1 s x+1 y+1]
% If input_surface = 1 then convert u into full [bid I J K] format

global nblade

if input_surface == 1
    u(1,:) = [u(1) reshape(surface{u(1), u(2)}(u(3),u(4),1:3),1,3)];
end

% Isolate bid, coords and connectivity list
b = u(1); I = u(2); J = u(3); K = u(4);
C = connectivity{b}{I+1,J+1,K+1};

% Read x, y, z coords of interest
coords_u = squeeze(coords{b}(I+1,J+1,K+1,:));
coords_u_2 = coords_u + [0 ; 0 ; single((2*pi())/nblade) * coords_u(2)];
coords_u_3 = coords_u - [0 ; 0 ; single((2*pi())/nblade) * coords_u(2)];

% Determine coords of co-incident connected nodes
u_connected = zeros(size(C,1),4);
coincident = zeros(size(C,1),1);
pitch = zeros(size(C,1),1);

for n = 1:size(C,1)
    nxbid = C(n,1); nxI = C(n,3); nxJ = C(n,4); nxK = C(n,5); 
    XYZ = squeeze(coords{nxbid+1}(nxI+1,nxJ+1,nxK+1,:));
       
    u_connected(n,:) = [nxbid+1 nxI nxJ nxK];
        
    coords_u_connected = XYZ;

    coincident(n) = isequal(coords_u_connected, coords_u) || sum(isalmost(coords_u_connected, coords_u_2, 1e-6))==3 ...
        || sum(isalmost(coords_u_connected, coords_u_3, 1e-6))==3;
    
    % Determine pitch index to enable co-ordinate change later
    if coincident(n) == 1
        if sum(isalmost(coords_u_connected, coords_u_2, 1e-6)) == 3
            pitch(n) = -1;
        elseif sum(isalmost(coords_u_connected, coords_u_2, 1e-6)) == 3
            pitch(n) = 1;
        end
    end
end

% Allow movement along diagonals or not
if diag_switch == 1
    A = [ 0 0 1 0 ; 0 0 1 1 ; 0 0 0 1 ; 0 0 -1 1 ;...
    0 0 -1 0 ; 0 0 -1 -1 ; 0 0 0 -1 ; 0 0 1 -1 ...
    
    0 1 1 0 ; 0 1 1 1 ; 0 1 0 1 ; 0 1 -1 1 ;...
    0 1 -1 0 ; 0 1 -1 -1 ; 0 1 0 -1 ; 0 1 1 -1 ...
    
    0 -1 1 0 ; 0 -1 1 1 ; 0 -1 0 1 ; 0 -1 -1 1 ;...
    0 -1 -1 0 ; 0 -1 -1 -1 ; 0 -1 0 -1 ; 0 -1 1 -1  ;...
    0 1 0 0 ; 0 -1 0 0];
else
    A = [ 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1 ; 0 -1 0 0 ; 0 0 -1 0 ; 0 0 0 -1 ];
end

% Initialise neighbours array
if coincident_switch == 1
    coincident_extras = sum(coincident);
else
    coincident_extras = sum(coincident) * size(A,1);
end
v_neighbours = zeros(coincident_extras + sum(coincident==0) + size(A,1), 5);

% Put in points surrounding u
n = size(A,1);
v_neighbours(1:n,:) = [A + repmat(u,size(A,1),1) zeros(size(A,1),1)];

% Put in non-coincident connected points
v_neighbours(n+1:n+sum(coincident==0),:) = [u_connected(coincident==0,:) zeros(sum(coincident==0),1)];
n = n+sum(coincident==0);

% Put in coincident points
if coincident_switch == 1
    v_neighbours(n+1:n+sum(coincident),:) = [u_connected(coincident == 1,:) pitch];
else
    % Put in points surrounding coincident ones
    coinc_index = find(coincident == 1);
    for m = 1:sum(coincident)
        v_neighbours(n+m*size(A,1)-(size(A,1)-1):n+m*size(A,1),:) = ...
            [(A + repmat(u_connected(coinc_index(m),:), size(A,1), 1)) ...
            repmat(pitch(coinc_index(m)) , size(A,1), 1)];
    end
end

% Remove neighbours that don't exist in domain at all
for n = size(v_neighbours,1):-1:1
    v = v_neighbours(n,:);
    try coords{v(1)}(v(2)+1,v(3)+1,v(4)+1,:);
    catch
        v_neighbours(n,:) = [];
    end
end

% Remove neighbours that don't exist in surface or not
if surface_switch == 1
    for m = size(v_neighbours,1): -1 : 1
        v = v_neighbours(m,:);
        if mwall{v(1)}(v(2)+1,v(3)+1,v(4)+1) == 1
            v_neighbours(m,:) = [];
        end
    end
    if return_surface == 1
        v_neighbours = ts_identify_surface(v_neighbours,surface,surface_type);
    end
end

end
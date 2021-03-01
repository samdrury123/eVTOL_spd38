function [surface, surface_type] = ts_find_faces(g)
% Function to find and plot mesh of surfaces in a given data structure

% Get bids
[bids] = ts_get_bids(g);
% Get coords
[~, coords_cart coords_pol primary] = ts_get_basics(g);

% Initialise cell arrays
surface = cell(length(bids),6);
surface_type = cell(length(bids),6);

% Loop over blocks, check six faces in turn and assign to a cell in 
% surface and surface_type
for bid = bids
    if isempty(coords_pol{bid+1}) == 0
        % I = 0 face
        Imin = 0; Imax = 0;
        Jmin = 0; Jmax = size(coords_pol{bid+1},2)-1;
        Kmin = 0; Kmax = size(coords_pol{bid+1},3)-1;

        [surface{bid+1,1} surface_type{bid+1,1}] =...
            ts_create_surface_array(Imin:Imax, Jmin:Jmax, Kmin:Kmax);

        % I = nI face
        Imin = size(coords_pol{bid+1},1)-1; Imax = size(coords_pol{bid+1},1)-1;
        Jmin = 0; Jmax = size(coords_pol{bid+1},2)-1;
        Kmin = 0; Kmax = size(coords_pol{bid+1},3)-1;

        [surface{bid+1,2} surface_type{bid+1,2}] =...
            ts_create_surface_array(Imin:Imax, Jmin:Jmax, Kmin:Kmax);

        % J = 0 face
        Imin = 0; Imax = size(coords_pol{bid+1},1)-1;
        Jmin = 0; Jmax = 0;
        Kmin = 0; Kmax = size(coords_pol{bid+1},3)-1;

        [surface{bid+1,3} surface_type{bid+1,3}] =...
            ts_create_surface_array(Imin:Imax, Jmin:Jmax, Kmin:Kmax);

        % J = nJ face
        Imin = 0; Imax = size(coords_pol{bid+1},1)-1;
        Jmin = size(coords_pol{bid+1},2)-1; Jmax = size(coords_pol{bid+1},2)-1;
        Kmin = 0; Kmax = size(coords_pol{bid+1},3)-1;

        [surface{bid+1,4} surface_type{bid+1,4}] =...
            ts_create_surface_array(Imin:Imax, Jmin:Jmax, Kmin:Kmax);

        % K = 0 face
        Imin = 0; Imax = size(coords_pol{bid+1},1)-1;
        Jmin = 0; Jmax = size(coords_pol{bid+1},2)-1;
        Kmin = 0; Kmax = 0;

        [surface{bid+1,5} surface_type{bid+1,5}] =...
            ts_create_surface_array(Imin:Imax, Jmin:Jmax, Kmin:Kmax);

        % K = nK face
        Imin = 0; Imax = size(coords_pol{bid+1},1)-1;
        Jmin = 0; Jmax = size(coords_pol{bid+1},2)-1;
        Kmin = size(coords_pol{bid+1},3)-1; Kmax = size(coords_pol{bid+1},3)-1;

        [surface{bid+1,6} surface_type{bid+1,6}] =...
            ts_create_surface_array(Imin:Imax, Jmin:Jmax, Kmin:Kmax);

    end        
    
end

% Loop over blocks again and assign cartesian and polar coordinates to
% surface data set
for bid = bids
    block_coords_cart = coords_cart{bid+1};
    block_coords_pol = coords_pol{bid+1};
    for s = 1:size(surface(bid+1,:),2)
        if isempty(surface{bid+1,s}) == 0
            surface_indices = surface{bid+1,s};
            surface_data = zeros(size(surface_indices,1),size(surface_indices,2),9);
            for x = 1:size(surface_indices,1)
                for y = 1:size(surface_indices,2)
                    % Read indices from surface array
                    node_indices = transpose(squeeze(surface_indices(x,y,:)));
                    I = node_indices(1); J = node_indices(2); K = node_indices(3);
                    % Read coords at given indices from coords arrays
                    % [ I J K x r rt x y z ]
                    node_coords_pol = transpose(squeeze(block_coords_pol(I+1, J+1, K+1, :)));
                    node_coords_cart = transpose(squeeze(block_coords_cart(I+1, J+1, K+1, :)));
                    surface_data(x,y,:) = double([node_indices node_coords_pol node_coords_cart ]);
                end
            end

            surface{bid+1,s} = surface_data;
            clear surface_data;
        end
    end
end
    

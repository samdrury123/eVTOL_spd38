function [Tree, T_surf, T_near] = ts_build_tree(g, phi_near)
% Function to build KD tree from given data in polar coordinates

% Set up indices ready
for bid = ts_get_bids(g)
    g{bid+1}.I = repmat(double(reshape(0:g{bid+1}.attribute.ni-1,[],1,1)),[1 g{bid+1}.attribute.nj g{bid+1}.attribute.nk]);
    g{bid+1}.J = repmat(double(reshape(0:g{bid+1}.attribute.nj-1,1,[],1)),[g{bid+1}.attribute.ni 1 g{bid+1}.attribute.nk]);
    g{bid+1}.K = repmat(double(reshape(0:g{bid+1}.attribute.nk-1,1,1,[])),[g{bid+1}.attribute.ni g{bid+1}.attribute.nj 1]);
end

% Initialise list of coords
T_surf = [];

% Loop over blocks and concatenate coords and indices into
% a list to feed into KD-Tree build function
for bid = ts_get_bids(g)
    nblade = double(g{bid+1}.bv.nblade);

    % Pull out surfaces only
    index = g{bid+1}.mwall == 0;
    x = double(g{bid+1}.x(index));
    r = double(g{bid+1}.r(index));
    rt = double(g{bid+1}.rt(index));
    
    I = g{bid+1}.I(index);
    J = g{bid+1}.J(index);
    K = g{bid+1}.K(index);
    BID = ones(size(I)) * bid;

    % Add co-ords of +/- a pitch with the same indices.
    temp = [BID I J K x r rt];
    n = size(temp,1);
    temp = repmat(temp,3,1);
    temp(n+1:2*n,end) = rt + (2 * pi() * r) / nblade;
    temp(2*n+1:3*n,end) = rt - (2 * pi() * r) / nblade;

    temp = [temp, [zeros(n,1) ; ones(n,1) ; -ones(n,1)] ones(3*n,1) * nblade];

    T_surf = [T_surf ; temp];
end

Tpol = single(T_surf(:,5:7));

Tree = kdtree(Tpol);

% Create list of all points in mesh too
% Initialise list of coords
T_near = [];

% Loop over blocks and concatenate coords and indices into
% a list to feed into KD-Tree build function
for bid = ts_get_bids(g)
    nblade = double(g{bid+1}.bv.nblade);
    
    % Pull out points near surfaces only based on phi
    index = g{bid+1}.phi < phi_near & g{bid+1}.mwall == 1;
    x = double(g{bid+1}.x(index));
    r = double(g{bid+1}.r(index));
    rt = double(g{bid+1}.rt(index));
    
    phi = double(g{bid+1}.phi(index));
    I = g{bid+1}.I(index);
    J = g{bid+1}.J(index);
    K = g{bid+1}.K(index);
    BID = ones(size(I)) * bid;

    % Add co-ords of +/- a pitch with the same indices.
    temp = [BID I J K x r rt phi];
    n = size(temp,1);
    temp = repmat(temp,3,1);
    temp(n+1:2*n,7) = rt + (2 * pi() * r) / nblade;
    temp(2*n+1:3*n,7) = rt - (2 * pi() * r) / nblade;

    temp = [temp(:,1:7), [zeros(n,1) ; ones(n,1) ; -ones(n,1)] temp(:,8) ones(3*n,1) * nblade];

    T_near = [T_near ; temp];
end

end
function g_crit = ts_get_crit(g,bids,plot_stuff)
% TS_GET_CRIT  Find critical points on the surface of a TURBOSTREAM mesh

% Default to plot critical points
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

% Default to analyse all blocks
if exist('bids','var') == 0 || isempty(bids) == 1
    bids = ts_get_bids(g);
end

% Number of iterations of Newton-Raphson process
nn = 20;

% Plot the surfaces of the mesh
if plot_stuff == 1
    C = ts_plot_surface(g,[0 0 0],[],1,1,bids);
else
    C = ts_get_surface(g);
end

% Pre-allocate critical point cell array
g_crit = cell(3,1);
varnames = {'x' 'r' 'rt' 'gam' 'lam'};
for n = 1:length(g_crit)
    for v = 1:length(varnames)
        g_crit{n}.(varnames{v}) = [];
    end
end

% Loop over all surface cuts
C = ts_secondary(C);
for m = 1:length(C)

    % Record cartesian velocities and coordinates
    xyz = cat(3,C{m}.x,C{m}.y,C{m}.z); Vxyz = cat(3,C{m}.Vx,C{m}.Vy_rel,C{m}.Vz_rel);
    xyz_av = 0.25 * (xyz(1:end-1,1:end-1,:) + xyz(2:end,1:end-1,:) + xyz(1:end-1,2:end,:) + ...
        xyz(2:end,2:end,:));
    xyz_corner = cat(4,xyz(1:end-1,1:end-1,:),xyz(2:end,1:end-1,:),xyz(1:end-1,2:end,:),...
        xyz(2:end,2:end,:));    

    % Calculate and normalise diagonal vectors in 3D
    a = xyz(2:end,2:end,:) - xyz(1:end-1,1:end-1,:); a = a ./ repmat(sum(a.^2,3).^0.5,[1 1 3]);
    b = xyz(1:end-1,2:end,:) - xyz(2:end,1:end-1,:); b = b ./ repmat(sum(b.^2,3).^0.5,[1 1 3]);

    % Resolve velocities into diagonal vectors
    V_corner = cat(4,Vxyz(1:end-1,1:end-1,:),Vxyz(2:end,1:end-1,:),Vxyz(1:end-1,2:end,:),...
        Vxyz(2:end,2:end,:));
    Va = squeeze(sum(V_corner .* repmat(a,[1 1 1 4]),3));
    Vb = squeeze(sum(V_corner .* repmat(b,[1 1 1 4]),3));

    % Generate initial guess in centre of cell and matrix to store critical point flags
    [ni,nj,~] = size(Va); p = ones(ni,nj,2) * 0.5; q = zeros(ni,nj,nn);

    % Apply Newton-Raphson method in two dimensions
    for n = 1:nn

        % Derivatives for current point
        dudx = Va(:,:,2) - Va(:,:,1) + p(:,:,2) .* (Va(:,:,1) - Va(:,:,3) - Va(:,:,2) + Va(:,:,4));
        dudy = Va(:,:,3) - Va(:,:,1) + p(:,:,1) .* (Va(:,:,1) - Va(:,:,3) - Va(:,:,2) + Va(:,:,4));
        dvdx = Vb(:,:,2) - Vb(:,:,1) + p(:,:,2) .* (Vb(:,:,1) - Vb(:,:,3) - Vb(:,:,2) + Vb(:,:,4));
        dvdy = Vb(:,:,3) - Vb(:,:,1) + p(:,:,1) .* (Vb(:,:,1) - Vb(:,:,3) - Vb(:,:,2) + Vb(:,:,4));

        % Assemble Jacobian for current point
        J = cat(3,dudx,dudy,dvdx,dvdy);

        % Invert Jacobian
        detJ = J(:,:,1) .* J(:,:,4) - J(:,:,2) .* J(:,:,3);
        inJ = cat(3,J(:,:,4),-J(:,:,2),-J(:,:,3),J(:,:,1)) ./ repmat(detJ,[1 1 4]);

        % Calculate velocities at current point
        u = Va(:,:,1) + p(:,:,1) .* (Va(:,:,2) - Va(:,:,1)) + p(:,:,2) .* (Va(:,:,3) - Va(:,:,1)) + ...
            p(:,:,1) .* p(:,:,2) .* (Va(:,:,1) - Va(:,:,3) - Va(:,:,2) + Va(:,:,4));
        v = Vb(:,:,1) + p(:,:,1) .* (Vb(:,:,2) - Vb(:,:,1)) + p(:,:,2) .* (Vb(:,:,3) - Vb(:,:,1)) + ...
            p(:,:,1) .* p(:,:,2) .* (Vb(:,:,1) - Vb(:,:,3) - Vb(:,:,2) + Vb(:,:,4));        

        % Find next point
        p = p - cat(3,sum(inJ(:,:,1:2) .* cat(3,u,v),3),sum(inJ(:,:,3:4) .* cat(3,u,v),3));

        % Record critical point locations
        q(:,:,n) = all(p > 0 & p < 1,3) & u < 1 & v < 1; 
        
    end

    % Record only points which are converged to lie within a cell
    q = all(q(:,:,end-4:end),3) == 1;
%     xyz_crit = reshape(xyz_av(repmat(q,[1 1 3])),[],3);
    xyz_crit = xyz_corner(:,:,:,1) + repmat(p(:,:,1),[1 1 3]) .* (xyz_corner(:,:,:,2) - ...
        xyz_corner(:,:,:,1)) + repmat(p(:,:,2),[1 1 3]) .* (xyz_corner(:,:,:,3) - ...
        xyz_corner(:,:,:,1)) + repmat(p(:,:,1),[1 1 3]) .* repmat(p(:,:,2),[1 1 3]) .* ...
        (xyz_corner(:,:,:,1) - xyz_corner(:,:,:,3) - xyz_corner(:,:,:,2) + xyz_corner(:,:,:,4));
    xyz_crit = reshape(xyz_crit(repmat(q,[1 1 3])),[],3);

    % Record Eigenvalues to determine class of point
    lam = dudx(q) + dvdy(q); gam = dudx(q) .* dvdy(q) - dudy(q) .* dvdx(q);

    % Group critical points by type: saddles, nodes and foci
    t{1} = gam <= 0; t{2} = gam > 0 & gam.^2 <= 4 * lam; t{3} = gam > 0 & gam.^2 > 4 * lam;

    % Append critical points to cell array
    for n = 1:length(g_crit)

        % Extract coordinates
        xyz = xyz_crit(t{n},:); xrrt = bl_cart2pol(xyz);
        g_crit{n}.x = [g_crit{n}.x ; xrrt(:,1)]; 
        g_crit{n}.r = [g_crit{n}.r ; xrrt(:,2)]; 
        g_crit{n}.rt = [g_crit{n}.rt ; xrrt(:,3)];

        % Record Eigenvalues
        g_crit{n}.gam = [g_crit{n}.gam ; gam(t{n})]; 
        g_crit{n}.lam = [g_crit{n}.lam ; lam(t{n})]; 

    end
end

% Update attributes in cell array and plot points
l = zeros(size(g_crit));
for n = 1:length(g_crit)
    
    % Record attributes
    g_crit{n}.attribute.ni = size(g_crit{n}.x,1); g_crit{n}.attribute.nj = 1; 
    g_crit{n}.attribute.nk = 1; g_crit{n}.attribute.bid = n-1; 

    % Plot points
    if plot_stuff == 1
        [y,z] = pol2cart(g_crit{n}.rt ./ g_crit{n}.r, g_crit{n}.r);
        l(n) = plot3(g_crit{n}.x,z,y,'.','markersize',20);
    end
end

% Add colour legend
if plot_stuff == 1
    legend(l,'Saddles','Nodes','Foci');    
end

    
end


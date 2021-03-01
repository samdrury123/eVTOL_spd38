function [g_bump, newfile] = ts_write_bump(g, u_bump, Lp, dp, normal, neighbour, T_surf, T_near,...
    phi_near, flow_type, directory, file_prefix, M, plot_stuff, hh_order)
% Function to add a bump to a geometry and write it to a hdf5 file
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end
if exist('hh_order','var') == 0
    hh_order = 5;
end

% Convert u_bump to 1:nk+1 notation from 0:nk (include bid)
bid_bump = u_bump(1);
I_bump = u_bump(2);
J_bump = u_bump(3);
K_bump = u_bump(4);

xrrt_bump = [g{bid_bump+1}.x(I_bump+1,J_bump+1,K_bump+1);...
    g{bid_bump+1}.r(I_bump+1,J_bump+1,K_bump+1);...
    g{bid_bump+1}.rt(I_bump+1,J_bump+1,K_bump+1)];

% Find all points on surface within a distance of Lp/2 in every dimension around u_bump
range = Lp/2;
index_vals_1 = find(T_surf(:,5) > xrrt_bump(1) - range & T_surf(:,5) < xrrt_bump(1) + range &...
    T_surf(:,6) > xrrt_bump(2) - range & T_surf(:,6) < xrrt_bump(2) + range &...
    T_surf(:,7) > xrrt_bump(3) - range & T_surf(:,7) < xrrt_bump(3) + range);

% Narrow down by computing straight line distance and comparing to Lp/2
xrrt_targets = T_surf(index_vals_1,5:7);
dist = sum((xrrt_targets - repmat(xrrt_bump.',length(index_vals_1),1)).^2,2);
index_vals_2 = find(dist < (range)^2);
[~, index_vals_3] = sort(dist(index_vals_2));
index_surfs = index_vals_1(index_vals_2(index_vals_3));

xrrt_surf = T_surf(index_surfs,5:7); normal_surf = zeros(size(xrrt_surf));
d_xrrt_surf = zeros(length(index_surfs),1);

% Find all points near to surface within a distance of Lp/2 in every dimension around u_bump
index_vals_1 = find(T_near(:,5) > xrrt_bump(1) - range & T_near(:,5) < xrrt_bump(1) + range &...
    T_near(:,6) > xrrt_bump(2) - range & T_near(:,6) < xrrt_bump(2) + range &...
    T_near(:,7) > xrrt_bump(3) - range & T_near(:,7) < xrrt_bump(3) + range);

% Narrow down by computing straight line distance and comparing to Lp/2
xrrt_targets = T_near(index_vals_1,5:7);
dist = sum((xrrt_targets - repmat(xrrt_bump.',length(index_vals_1),1)).^2,2);
index_vals_2 = find(dist < (range)^2);
[~, index_vals_3] = sort(dist(index_vals_2));
index_nears = index_vals_1(index_vals_2(index_vals_3));

% Initialise new data array
g_bump = g;
disp([num2str(length(index_surfs)) ' Points in Range on Surface']);
out_bound = 0;

% Create arrays to drive A* search [ Q/visited  e/distance  h/straight_dist
% f/A*total  previous bid    previous I    previous J    previous K]
if exist([directory file_prefix 'bump_' num2str(bid_bump) '_' num2str(I_bump) '_'...
    num2str(J_bump) '_' num2str(K_bump) '_' num2str(Lp) '_A.mat'],'file') == 0
    A = cell(size(g));
    for bid = ts_get_bids(g)
        ni = g{bid+1}.attribute.ni;
        nj = g{bid+1}.attribute.nj;
        nk = g{bid+1}.attribute.nk;

        A{bid+1} = repmat(reshape([1 inf inf inf 0 0 0 0],1,1,1,8),[ni nj nk 1]);
    end
    A_calculated = 0;
else
    load([directory file_prefix 'bump_' num2str(bid_bump) '_' num2str(I_bump) '_'...
    num2str(J_bump) '_' num2str(K_bump) '_' num2str(Lp) '_A.mat']);
    A_calculated = 1;
end

for n = 1:length(index_surfs);

    % Take target point from list of values found from range query
    index_surf = index_surfs(n);
    u_target = T_surf(index_surf,1:4);
    bid_target = u_target(1);
    I_target = u_target(2);
    J_target = u_target(3);
    K_target = u_target(4);

    xrrt_target = [g{bid_target+1}.x(I_target+1,J_target+1,K_target+1);...
        g{bid_target+1}.r(I_target+1,J_target+1,K_target+1);...
        g{bid_target+1}.rt(I_target+1,J_target+1,K_target+1)];

    % Find distances in the surface from the centre of the bump to
    % the target using A* algorithm & backtrack along the path found
    % Calculate heuristic function (shortest straight line distance first
    if A_calculated == 0
        xrrt_vs = T_surf(index_surfs,5:7);
        hs = sum((xrrt_vs - repmat(xrrt_target.',length(index_surfs),1)).^2,2).^0.5;
        for m = 1:length(index_surfs)
            v = T_surf(index_surfs(m),1:4);
            h = hs(m);
            f = h + A{v(1)+1}(v(2)+1,v(3)+1,v(4)+1,2);
            A{v(1)+1}(v(2)+1,v(3)+1,v(4)+1,4) = f;
            A{v(1)+1}(v(2)+1,v(3)+1,v(4)+1,3) = h;
        end

        [A path] = ts_A_star(g, u_bump, u_target, neighbour, A, Lp/2, 0);
    end

    dist = A{u_target(1)+1}(u_target(2)+1,u_target(3)+1,u_target(4)+1,2);

    if dist<Lp/2
        if plot_stuff == 1
            disp(['Distance from bump to target = ' num2str(dist)]);
        end
        % Calculate bump height and assign target coords & normal
        d = dp * (cos( pi() * dist / Lp))^hh_order;

        target_normal = squeeze(normal{bid_target+1}(I_target+1,J_target+1,K_target+1,:));

%         % Perturb mesh by bump height * unit normal and write into new
%         % coords, determine closest points normal to the surface and move
%         % by decreasing amount
%         
%         v_normals = normal_neighbours(g,u_target);
% 
%         % Loop over normal points and move mesh
%         for m = 1:size(v_normals,1)
%             v = v_normals(m,:);
%             
%             xrrt_v = [g{v(1)+1}.x(v(2)+1,v(3)+1,v(4)+1);...
%                 g{v(1)+1}.r(v(2)+1,v(3)+1,v(4)+1);...
%                 g{v(1)+1}.rt(v(2)+1,v(3)+1,v(4)+1)];
%             q = xrrt_v - xrrt_target;
% 
%             z = dot(q,target_normal);
% 
%             % Perturb mesh by d*exp(-z/d)
%             xrrt_v = xrrt_v + d * target_normal * 2^((-2 * z) / d);
% 
%             % Update bumped coordinate cell array
%             g_bump{v(1)+1}.x(v(2)+1,v(3)+1,v(4)+1) = xrrt_v(1);
%             g_bump{v(1)+1}.r(v(2)+1,v(3)+1,v(4)+1) = xrrt_v(2);
%             g_bump{v(1)+1}.rt(v(2)+1,v(3)+1,v(4)+1) = xrrt_v(3);
% 
%         end

        % Perturb surface of mesh by bump height * unit normal.
        g_bump{bid_target+1}.x(I_target+1,J_target+1,K_target+1) = ...
            g{bid_target+1}.x(I_target+1,J_target+1,K_target+1) + d * target_normal(1);
        g_bump{bid_target+1}.r(I_target+1,J_target+1,K_target+1) = ...
            g{bid_target+1}.r(I_target+1,J_target+1,K_target+1) + d * target_normal(2);
        g_bump{bid_target+1}.rt(I_target+1,J_target+1,K_target+1) = ...
            g{bid_target+1}.rt(I_target+1,J_target+1,K_target+1) + d * target_normal(3);
        
        d_xrrt_surf(n,1) = d; normal_surf(n,:) = target_normal;

        out_bound = 0;
    else
        out_bound = out_bound + 1;
        % After 10% of the list has been covered without a target in
        % actual distance break from search
        if out_bound > length(index_surfs) * 0.1
            disp('Assume following points lie out of bounds')
            break
        end
    end
end

% Write A into .mat file for future use
save([directory file_prefix 'bump_' num2str(bid_bump) '_' num2str(I_bump) '_'...
    num2str(J_bump) '_' num2str(K_bump) '_' num2str(Lp) '_A.mat'],'A','-v7.3');

% Look at all close mesh points, find the nearest on the surface and move by a decaying
% distance determined by phi

% xrrt_surf = T_surf(index_surfs,5:7); xrrt_near = T_near(index_near,5:7);
% n_surf = length(index_surfs); n_near = length(index_near);
% xrrt_surf = repmat(reshape(xrrt_surf,[n_surf 1 3]),[1 n_near 1]); 
% xrrt_near = repmat(reshape(xrrt_near,[1 n_near 3]),[n_surf 1 1]);

% Pull out near and surface coordinates and construct field to interpolate
% phi
n_extrap = 12; l_extrap = 5e-4;
xrrt_near = T_near(index_nears,5:7); xrrt_surf = T_surf(index_surfs,5:7);

% Change coordinates back by a pitch
xrrt_surf(:,3) = xrrt_surf(:,3) - T_surf(index_surfs,8) .* (2 * pi ./ T_surf(index_surfs,9)) .* xrrt_surf(:,2);
xrrt_near(:,3) = xrrt_near(:,3) - T_near(index_nears,8) .* (2 * pi ./ T_near(index_nears,10)) .* xrrt_near(:,2);

n_surf = length(index_surfs); n_near = length(index_nears);
F_phi = TriScatteredInterp(xrrt_near(:,1),xrrt_near(:,2),xrrt_near(:,3),T_near(index_nears,9));

% Calculate coordinates of n_extrap new points directly normal to the
% surface
d_xrrt_extrap = repmat(d_xrrt_surf,[1 n_extrap]); 
l_extrap = linspace(0,l_extrap,n_extrap); 
normal_extrap = repmat(reshape(normal_surf,[n_surf 1 3]),[1 n_extrap 1]);
xrrt_extrap = repmat(reshape(xrrt_surf,[n_surf 1 3]),[1 n_extrap 1]) +...
    repmat(l_extrap,[n_surf 1 3]) .* normal_extrap;

% Reshape new region of extrapolation into lists
d_xrrt_extrap = reshape(d_xrrt_extrap,[],1);
normal_extrap = reshape(normal_extrap,[],3);
xrrt_extrap = reshape(xrrt_extrap,[],3);

% Re calculate shift based on wall distance at extrapolated points
phi_extrap = F_phi(xrrt_extrap(:,1),xrrt_extrap(:,2),xrrt_extrap(:,3));
d_xrrt_extrap = d_xrrt_extrap .* (cos(pi * phi_extrap / phi_near) + 1) / 2;

% Calculate normals and shift at actual near grid points
F_d_xrrt = TriScatteredInterp(xrrt_extrap(:,1),xrrt_extrap(:,2),xrrt_extrap(:,3),d_xrrt_extrap);
F_x_normal = TriScatteredInterp(xrrt_extrap(:,1),xrrt_extrap(:,2),xrrt_extrap(:,3),normal_extrap(:,1));
F_r_normal = TriScatteredInterp(xrrt_extrap(:,1),xrrt_extrap(:,2),xrrt_extrap(:,3),normal_extrap(:,2));
F_rt_normal = TriScatteredInterp(xrrt_extrap(:,1),xrrt_extrap(:,2),xrrt_extrap(:,3),normal_extrap(:,3));

d_xrrt_near = F_d_xrrt(xrrt_near(:,1),xrrt_near(:,2),xrrt_near(:,3));
x_normal_near = F_x_normal(xrrt_near(:,1),xrrt_near(:,2),xrrt_near(:,3));
r_normal_near = F_r_normal(xrrt_near(:,1),xrrt_near(:,2),xrrt_near(:,3));
rt_normal_near = F_rt_normal(xrrt_near(:,1),xrrt_near(:,2),xrrt_near(:,3));

% Move actual grid points and update the whole cell array
normal_near = [x_normal_near r_normal_near rt_normal_near];
normal_near(isnan(normal_near(:,1)),:) = 0;
d_xrrt_near(isnan(d_xrrt_near)) = 0;

xrrt_near = xrrt_near + repmat(d_xrrt_near, [1 3]) .* normal_near;

for n = 1:n_near
    index_near = index_nears(n);
    u = T_near(index_near,1:4);
    bid = u(1);    I = u(2);    J = u(3);    K = u(4);
    
    g_bump{bid+1}.x(I+1,J+1,K+1) = xrrt_near(n,1);
    g_bump{bid+1}.r(I+1,J+1,K+1) = xrrt_near(n,2);
    g_bump{bid+1}.rt(I+1,J+1,K+1) = xrrt_near(n,3);
end

% Plot out both surfaces in polar coordinates to inspect bump
if plot_stuff == 1
    ts_plot_surface(g, [0 0 0],1);
    ts_plot_surface(g_bump, [1 0 0],1);
    ts_plot_surface(g, [0 0 0],1);
    plot3(T_surf(index_surfs,5),T_surf(index_surfs,7),T_surf(index_surfs,6),'k.')    
end

% Set new filename based on position and size of the bump
newfile = [file_prefix 'M0' num2str(100*M) '_t' num2str(flow_type) '_bump_'...
    num2str(bid_bump) '_' num2str(I_bump) '_'...
    num2str(J_bump) '_' num2str(K_bump) '_' num2str(Lp) '_'...
    num2str(dp) '_input.hdf5'];

% Write bumped cell array to turbostream input file  
delete([directory newfile]);
ts_write_hdf5(g_bump, [directory newfile]);
ts_write_xdmf(g_bump, strrep([directory newfile],'.hdf5','.xdmf'));

end

function v = normal_neighbours(g,u)
% Determine points to be shifted in normal direction dependant on position
% of walls

% Set number of points away from wall
n_points = 5;

% Determine whether it lies near block edges or walls
nijk = [g{u(1)+1}.attribute.ni g{u(1)+1}.attribute.nj g{u(1)+1}.attribute.nk];
IJK = u(2:4);

edge_dim = [];
for dim = 1:3
    if IJK(dim) == 0 || IJK(dim) == nijk(dim)-1
        edge_dim = [edge_dim dim];
    end
end

% size edge_dim is one, then the point is in the centre of the block, 5
% points away from the surface should be sufficient to accommodate the bump
if length(edge_dim) == 1
    IJK_v = zeros(n_points,3);
    IJK_v(:,edge_dim) = 0:n_points-1;
    if IJK(edge_dim) == 0
        v = repmat(u,n_points,1) + [zeros(n_points,1) IJK_v];
    else
        v = repmat(u,n_points,1) - [zeros(n_points,1) IJK_v];
    end
    
% Determine whether it is in the corner of a wall, or just on the edge
% of the block    
elseif length(edge_dim) == 2

    IJK_v1 = zeros(1,3);
    IJK_v2 = zeros(1,3);
    IJK_v1(edge_dim(1)) = 1;
    IJK_v2(edge_dim(2)) = 1;
    
    if IJK(edge_dim(1)) == 0
        v1 = u + [0 IJK_v1];
    else
        v1 = u - [0 IJK_v1];
    end
    
    if IJK(edge_dim(2)) == 0
        v2 = u + [0 IJK_v2];
    else
        v2 = u - [0 IJK_v2];
    end    
        
    % The point in question is in the corner of a wall, leave all other
    % points alone
    if g{v1(1)+1}.mwall(v1(2)+1,v1(3)+1,v1(4)+1) == 0 &&...
        g{v2(1)+1}.mwall(v2(2)+1,v2(3)+1,v2(4)+1) == 0
        v = u;
        return
    end
    
    % The point is in the corner of a block but not a wall, determine
    % normal neighbours appropriately
    if g{v1(1)+1}.mwall(v1(2)+1,v1(3)+1,v1(4)+1) == 0
        edge_dim = edge_dim(2);
    else
        edge_dim = edge_dim(1);
    end
    
    IJK_v = zeros(n_points,3);
    IJK_v(:,edge_dim) = 0:n_points-1;
    if IJK(edge_dim) == 0
        v = repmat(u,n_points,1) + [zeros(n_points,1) IJK_v];
    else
        v = repmat(u,n_points,1) - [zeros(n_points,1) IJK_v];
    end

% The point is a vertex of a block, is it also in a wall corner or not
elseif length(edge_dim) == 3
    IJK_v1 = zeros(1,3);
    IJK_v2 = zeros(1,3);
    IJK_v3 = zeros(1,3);
    IJK_v1(edge_dim(1)) = 1;
    IJK_v2(edge_dim(2)) = 1;
    IJK_v3(edge_dim(3)) = 1;

    if IJK(edge_dim(1)) == 0
        v1 = u + [0 IJK_v1];
    else
        v1 = u - [0 IJK_v1];
    end

    if IJK(edge_dim(2)) == 0
        v2 = u + [0 IJK_v2];
    else
        v2 = u - [0 IJK_v2];
    end  
    
    if IJK(edge_dim(3)) == 0
        v3 = u + [0 IJK_v3];
    else
        v3 = u - [0 IJK_v3];
    end       
    
    % Is it a wall corner?
    if sum([g{v1(1)+1}.mwall(v1(2)+1,v1(3)+1,v1(4)+1) == 0, ... 
        g{v2(1)+1}.mwall(v2(2)+1,v2(3)+1,v2(4)+1) == 0, ...
        g{v3(1)+1}.mwall(v3(2)+1,v3(3)+1,v3(4)+1) == 0]) > 2
        v = u;
        return
    end        
    
    % Which direction is normal to the wall
    if g{v1(1)+1}.mwall(v1(2)+1,v1(3)+1,v1(4)+1) == 0 && ...
        g{v2(1)+1}.mwall(v2(2)+1,v2(3)+1,v2(4)+1) == 0
        edge_dim = edge_dim(3);
    elseif g{v1(1)+1}.mwall(v1(2)+1,v1(3)+1,v1(4)+1) == 0 && ...
        g{v3(1)+1}.mwall(v3(2)+1,v3(3)+1,v3(4)+1) == 0
        edge_dim = edge_dim(2);
    else
        edge_dim = edge_dim(1);
    end
    
    IJK_v = zeros(n_points,3);
    IJK_v(:,edge_dim) = 0:n_points-1;
    if IJK(edge_dim) == 0
        v = repmat(u,n_points,1) + [zeros(n_points,1) IJK_v];
    else
        v = repmat(u,n_points,1) - [zeros(n_points,1) IJK_v];
    end
        
end

end
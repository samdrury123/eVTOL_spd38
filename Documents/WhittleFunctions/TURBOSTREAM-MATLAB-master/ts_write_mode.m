function [g_mode filename] = ts_write_mode_2(g, f, dp, limit, normal, neighbour, T_surf, T_near,...
    flow_type, directory, file_prefix, M, mode_prefix, plot_stuff, smooth_f)
% Add a mode shape to a blade geometry to target a loss generation
% mechanism

if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end
if exist('smooth_f','var') == 0
    smooth_f = 1;
end

%% Set Limits on f
% Limit the scope of interest based on indices or coordinates
% Set boundaries to have loss ratio of 1 for smooth edges
varnames = fieldnames(limit);
for v = 1:length(varnames)
    varname = varnames{v};
    for n = 1:size(limit.(varname),1)
        if limit.(varname)(n,2) == 1
            f.Loss(f.(varname) < limit.(varname)(n,1)) = 1;
        elseif limit.(varname)(n,2) == -1
            f.Loss(f.(varname) > limit.(varname)(n,1)) = 1; 
        end
    end
end
% Set areas of loss ratio to 0 for total removal later
for v = 1:length(varnames)
    varname = varnames{v};
    for n = 1:size(limit.(varname),1)
        if limit.(varname)(n,2) == 2
            f.Loss(f.(varname) < limit.(varname)(n,1)) = 0;
        elseif limit.(varname)(n,2) == -2
            f.Loss(f.(varname) > limit.(varname)(n,1)) = 0;
        end
    end
end

% Remove unwanted areas of mode
f.bid(f.Loss == 0) = [];
f.I(f.Loss == 0) = [];
f.J(f.Loss == 0) = [];
f.K(f.Loss == 0) = [];
f.x(f.Loss == 0) = [];
f.r(f.Loss == 0) = [];
f.rt(f.Loss == 0) = [];
f.Loss(f.Loss == 0) = [];

%% Interpolate intensity - Q
% Initialise mode intensity parameter - Q
Q = cell(size(g));
for bid = ts_get_bids(g) 
    Q{bid+1} = zeros(size(g{bid+1}.x));
end

% Construct interpolant for Q, Q is loss minus 1
% Find new tangential coordinates and normal vectors - a
a = zeros(size(f.Loss,1),3); rt1 = zeros(size(f.Loss)); rt2 = zeros(size(f.Loss));
bids = f.bid; I = f.I; J = f.J; K = f.K; q = f.Loss - 1;
x = f.x; r = f.r; rt = f.rt;
for n = 1:length(q)
    nblade = double(g{bids(n)+1}.bv.nblade);

    rt1(n) = double(g{bids(n)+1}.rt(I(n)+1,J(n)+1,K(n)+1)) + (2*pi*r(n)) / nblade;
    rt2(n) = double(g{bids(n)+1}.rt(I(n)+1,J(n)+1,K(n)+1)) - (2*pi*r(n)) / nblade;

    a(n,:) = double(reshape(normal{bids(n)+1}(I(n)+1,J(n)+1,K(n)+1,:),[1 3]));
end

% Repeat by +/- a pitch 
x = repmat(x,[3 1]); r = repmat(r,[3 1]); rt = [rt ; rt1 ; rt2];
q = repmat(q,[3 1]); a = repmat(a,[3 1]);

% Repeat by +/- normals
if smooth_f == 1
    d_wall = double(g{1}.r(1,2,1) - g{1}.r(1,1,1));
    x = repmat(x,[3 1]) + [d_wall * a(:,1) ; zeros(size(x)) ; - d_wall * a(:,1)];
    r = repmat(r,[3 1]) + [d_wall * a(:,2) ; zeros(size(r)) ; - d_wall * a(:,2)];
    rt = repmat(rt,[3 1]) + [d_wall * a(:,3) ; zeros(size(rt)) ; - d_wall * a(:,3)];
    q = repmat(q,[3 1]);

    F_Q = TriScatteredInterp(x,r,rt,q);
end
% Interpolate new intensities
for bid = unique(bids).'
    [I,J,K] = ind2sub(size(g{bid+1}.mwall), find(g{bid+1}.mwall == 0));
    I = I-1; J = J-1; K = K-1;

    x_interp = zeros(length(I),1); r_interp = zeros(length(I),1); 
    rt_interp = zeros(length(I),1);
    for n = 1:length(I)
        x_interp(n) = g{bid+1}.x(I(n)+1,J(n)+1,K(n)+1);
        r_interp(n) = g{bid+1}.r(I(n)+1,J(n)+1,K(n)+1);
        rt_interp(n) = g{bid+1}.rt(I(n)+1,J(n)+1,K(n)+1);
    end

    if smooth_f == 1
        Q_interp = F_Q(x_interp,r_interp,rt_interp);

        for n = 1:length(Q_interp)
            Q{bid+1}(I(n)+1,J(n)+1,K(n)+1) = Q_interp(n);
        end

        % Replace those outside range with no movement
        Q{bid+1}(isnan(Q{bid+1})) = 0;
    end
end
if smooth_f == 0;
    for n = 1:length(f.Loss)
        Q{f.bid(n)+1}(f.I(n)+1,f.J(n)+1,f.K(n)+1) = f.Loss(n) - 1;
    end
end
    

% if plot_stuff == 1
%     H.surface = figure();
%     ts_plot_surface(g, [0 0 0], H.surface)
%     ts_plot_surface(g, Q, H.surface, 0, 2)
%     plot3(f.x,f.rt,f.r,'k.')
% %     plot3(x,rt,r,'b.')
%     colorbar
%     caxis([-1 1])
% end

%% Smooth intensity - Q
% Pull out points on all surfaces that are near those with non-zero Q
% intensity
d_length = max(max(abs(g{1}.x(1:end-1,:,1) - g{1}.x(2:end,:,1)))) * 2;
T_move = [];
for bid = ts_get_bids(g)
    [i,j,k] = ind2sub(size(Q{bid+1}), find(Q{bid+1} ~= 0));
    i = i-1; j = j-1; k = k-1;
    
    x = zeros(size(i)); r = zeros(size(i)); rt = zeros(size(i)); % q = zeros(size(i));
    for n = 1:length(i)
        x(n) = g{bid+1}.x(i(n)+1,j(n)+1,k(n)+1);
        r(n) = g{bid+1}.r(i(n)+1,j(n)+1,k(n)+1);
        rt(n) = g{bid+1}.rt(i(n)+1,j(n)+1,k(n)+1);
%         q(n) = Q{bid+1}(i(n)+1,j(n)+1,k(n)+1);
    end
    
    T_move = [T_move ; ones(size(i))*bid i j k x r rt];
end

index_Q = find(...
    T_surf(:,5) > min(T_move(:,5)) - d_length & T_surf(:,5) < max(T_move(:,5)) + d_length &...
    T_surf(:,6) > min(T_move(:,6)) - d_length & T_surf(:,6) < max(T_move(:,6)) + d_length &...
    T_surf(:,7) > min(T_move(:,7)) - d_length & T_surf(:,7) < max(T_move(:,7)) + d_length);
x = T_surf(index_Q,5);  r = T_surf(index_Q,6);  rt = T_surf(index_Q,7);

% Get intensities at surface points
% if smooth_f == 1
    q = zeros(size(x));  bidIJK = T_surf(index_Q,1:4);
    for n = 1:length(index_Q)
        q(n) = Q{bidIJK(n,1) + 1}(bidIJK(n,2)+1,bidIJK(n,3)+1,bidIJK(n,4)+1);
    end

    % Call smoothing routine
    % [q c] = smooth_q(x,r,rt,q,d_length,12);
    q = smooth_IJK(bidIJK,q,neighbour,2,0);

    % Reassign points back into Q
    for n = 1:length(index_Q)
        m = sum(repmat(bidIJK(n,:),[length(index_Q) 1]) == bidIJK ,2) == 4;
        Q{bidIJK(n,1) + 1}(bidIJK(n,2)+1,bidIJK(n,3)+1,bidIJK(n,4)+1) = sum(q(m))/sum(m);
    end
% end
    
% Non-dimensionalise Q
Q_max = 0; Q_min = 0;
for n = 1:length(Q)
    Q_max = max(max(max(max(Q{n}))),Q_max);
    Q_min = min(min(min(min(Q{n}))),Q_min);
end

Q_max = max(abs(Q_min),Q_max);
for n = 1:length(Q)
    Q{n} = Q{n} / Q_max;
end

% Plot output of intensities
if plot_stuff == 1
    H.surface = figure();
    ts_plot_surface(g, [0 0 0], H.surface)
    ts_plot_surface(g, Q, H.surface, 0, 2)
    plot3(f.x,f.rt,f.r,'k.')
%     plot3(x,rt,r,'b.')
    colorbar
    caxis([-1 1])
end

%% Define normals away from wall and smooth field
% Find points on the wall by their non-zero Q
T_move = [];
for bid = ts_get_bids(g)
    [i,j,k] = ind2sub(size(Q{bid+1}), find(Q{bid+1} ~= 0));
    i = i-1; j = j-1; k = k-1;

    x = zeros(size(i)); r = zeros(size(i)); rt = zeros(size(i)); 
    a = zeros([length(i) 3]); d = zeros(size(i));
    for n = 1:length(i)
        x(n) = g{bid+1}.x(i(n)+1,j(n)+1,k(n)+1);
        r(n) = g{bid+1}.r(i(n)+1,j(n)+1,k(n)+1);
        rt(n) = g{bid+1}.rt(i(n)+1,j(n)+1,k(n)+1);
        a(n,:) = reshape(normal{bid+1}(i(n)+1,j(n)+1,k(n)+1,:),[1 3]);
        d(n) = Q{bid+1}(i(n)+1,j(n)+1,k(n)+1) * dp;
    end

    T_move = [T_move ; ones(size(i))*bid i j k x r rt];
end

% Find points, normals and displacements on the surface by their proximity
% to previous points with non-zero Q
index_surf = find(...
    T_surf(:,5) > min(T_move(:,5)) - d_length & T_surf(:,5) < max(T_move(:,5)) + d_length &...
    T_surf(:,6) > min(T_move(:,6)) - d_length & T_surf(:,6) < max(T_move(:,6)) + d_length &...
    T_surf(:,7) > min(T_move(:,7)) - d_length & T_surf(:,7) < max(T_move(:,7)) + d_length);

n_surf = length(index_surf);
normal_surf = zeros(n_surf,3); d_xrrt_surf = zeros(n_surf,1);
xrrt_surf = T_surf(index_surf,5:7); bidIJK_surf = T_surf(index_surf,1:4);
pitch_surf = T_surf(index_surf,8);
% xrrt_surf = xrrt_surf(pitch_surf == 0,:); bidIJK_surf = bidIJK_surf(pitch_surf == 0,:);
% normal_surf = normal_surf(pitch_surf == 0,:); d_xrrt_surf = d_xrrt_surf(pitch_surf == 0,:);

for n = 1:n_surf
    normal_surf(n,:) = reshape(normal{bidIJK_surf(n,1)+1}...
        (bidIJK_surf(n,2)+1,bidIJK_surf(n,3)+1,bidIJK_surf(n,4)+1,:),[1 3]);
    d_xrrt_surf(n) = Q{bidIJK_surf(n,1)+1}...
        (bidIJK_surf(n,2)+1,bidIJK_surf(n,3)+1,bidIJK_surf(n,4)+1,:) * dp;
end

% Remove points that are not moved much
tol = 1e-8;
% n_remove = abs(d_xrrt_surf) < tol;
% xrrt_surf(n_remove,:) = [];
% normal_surf(n_remove,:) = [];
% bidIJK_surf(n_remove,:) = [];
% d_xrrt_surf(n_remove,:) = [];
% n_surf = size(xrrt_surf,1);

% Find points away from the surface but near to those with non-zero Q
index_near = find(...
    T_near(:,5) > min(T_move(:,5)) - d_length & T_near(:,5) < max(T_move(:,5)) + d_length &...
    T_near(:,6) > min(T_move(:,6)) - d_length & T_near(:,6) < max(T_move(:,6)) + d_length &...
    T_near(:,7) > min(T_move(:,7)) - d_length & T_near(:,7) < max(T_move(:,7)) + d_length);

xrrt_near = T_near(index_near,5:7); bidIJK_near = T_near(index_near,1:4); 
phi_near = T_near(index_near,9); pitch_near = T_near(index_near,8);
xrrt_near = xrrt_near(pitch_near == 0,:); bidIJK_near = bidIJK_near(pitch_near == 0,:);
phi_near = phi_near(pitch_near == 0,:); n_near = size(xrrt_near,1);

% Define normals on the wall by using a surface of constant phi
phis = logspace(-8,log10(max(phi_near)*1.2),10); 
F_phi = TriScatteredInterp(xrrt_near(:,1),xrrt_near(:,2),xrrt_near(:,3),phi_near);
d = 1e-4;

% Get triangulated surfaces and normals to each triangle
FV = cell(length(phis),1); 
for n = 1:length(phis)
    FV{n}.xrrt = []; FV{n}.normal = []; 
    for bid = unique(bidIJK_near(:,1)).'
        fv = isosurface(g{bid+1}.x,g{bid+1}.r,g{bid+1}.rt,g{bid+1}.phi,phis(n));
        xrrt = (fv.vertices(fv.faces(:,1),:) + fv.vertices(fv.faces(:,2),:) +...
            fv.vertices(fv.faces(:,3),:)) / 3; xrrt = double(xrrt);
        FV{n}.xrrt = [FV{n}.xrrt ; xrrt];
        normal_f = cross(fv.vertices(fv.faces(:,1),:) - fv.vertices(fv.faces(:,2),:),...
            fv.vertices(fv.faces(:,2),:) - fv.vertices(fv.faces(:,3),:));
        normal_f = normal_f ./ repmat(sum(normal_f.^2,2).^0.5, [1 3]); normal_f = double(normal_f);
        
        % Make sure normal_f points outwards using tri scattered interp
        phi_1 = F_phi(xrrt(:,1),xrrt(:,2),xrrt(:,3));
        phi_2 = F_phi(xrrt(:,1)+d*normal_f(:,1), xrrt(:,2)+d*normal_f(:,2), xrrt(:,3)+d*normal_f(:,3));
        s = phi_2 > phi_1; s(s~=1) = -1;
        normal_f = - normal_f .* repmat(s,[1 3]);
        
        FV{n}.normal = [FV{n}.normal ; normal_f];
    end
end

% normal_fs inside surface
xrrt = FV{1}.xrrt - 1e-5 * FV{1}.normal; fv = cell(1);
normal_f = FV{1}.normal; fv{1}.xrrt = xrrt; fv{1}.normal = normal_f;
FV = [fv ; FV];

% Vector of all normals and triangle centres
normal_fv = []; xrrt_fv = [];
for n = 1:length(FV)
    normal_fv = [normal_fv ; FV{n}.normal];
    xrrt_fv = [xrrt_fv ; FV{n}.xrrt];
end

% Delete those not in volume
i = xrrt_fv(:,1) > min(T_move(:,5)) - d_length & xrrt_fv(:,1) < max(T_move(:,5)) + d_length &...
    xrrt_fv(:,2) > min(T_move(:,6)) - d_length & xrrt_fv(:,2) < max(T_move(:,6)) + d_length &...
    xrrt_fv(:,3) > min(T_move(:,7)) - d_length & xrrt_fv(:,3) < max(T_move(:,7)) + d_length;
xrrt_fv = xrrt_fv(i,:); normal_fv = normal_fv(i,:);

% Calculate distances to move at each triangle
d_xrrt_fv = zeros(size(xrrt_fv,1),1); n_fv = size(xrrt_fv,1);
% [~, m] = min(abs( sum((repmat(reshape(xrrt_fv,[n_fv 1 3]), [1 1 n_surf]) -...
%     repmat(reshape(xrrt_surf,[1 n_surf 3]), [n_fv 1 1])).^2,3).^0.5),[],1);

for n = 1:size(xrrt_fv,1)
%     i = abs(xrrt_surf(:,1) - xrrt_fv(n,1)) < d_length & abs(xrrt_surf(:,2) - xrrt_fv(n,2)) < d_length &...
%         abs(xrrt_surf(:,3) - xrrt_fv(n,3)) < d_length;
%     d_xrrt_temp = d_xrrt_surf(i);  xrrt_temp = xrrt_surf(i,:); n_temp = sum(i);
    
    % Closest straight distance
%     [~, m] = min(abs(sum( (repmat(xrrt_fv(n,:),n_temp,1) - xrrt_temp).^2, 2)));
    [~, m] = min(abs(sum( (repmat(xrrt_fv(n,:),n_surf,1) - xrrt_surf).^2, 2)));
    
    % Closest point to a line defined by the normal
%     [~, m] = min(abs(sum( ((repmat(xrrt_fv(n,:),n_surf,1) - xrrt_surf) - ...
%         repmat(dot(repmat(xrrt_fv(n,:),n_surf,1) - xrrt_surf,repmat(normal_fv(n,:),n_surf,1),2),1,3) .*...
%         repmat(normal_fv(n,:),n_surf,1)).^2, 2).^0.5));
%     d_xrrt_fv(n) = d_xrrt_temp(m);
    d_xrrt_fv(n) = d_xrrt_surf(m);
    
end

% Remove triangles which don't move much
n_remove = abs(d_xrrt_fv) < tol;
xrrt_fv(n_remove,:) = []; normal_fv(n_remove,:) = []; d_xrrt_fv(n_remove,:) = [];
n_fv = size(xrrt_fv,1); 

% Interpolants of normals
% F_nx = TriScatteredInterp(xrrt_fv(:,1),xrrt_fv(:,2),xrrt_fv(:,3),normal_fv(:,1));
% F_nr = TriScatteredInterp(xrrt_fv(:,1),xrrt_fv(:,2),xrrt_fv(:,3),normal_fv(:,2));
% F_nrt = TriScatteredInterp(xrrt_fv(:,1),xrrt_fv(:,2),xrrt_fv(:,3),normal_fv(:,3));

% % Reassign back into normal cell array to check by plotting a displaced
% % surface
% if plot_stuff == 1
%     normal_smooth = normal;
%     for n = 1:n_surf
%         normal_smooth{bidIJK_surf(n,1)+1}(bidIJK_surf(n,2)+1,bidIJK_surf(n,3)+1,bidIJK_surf(n,4)+1,:) =...
%             reshape(normal_surf(n,:),[1 1 1 3]);
%     end
%     figure()
%     hold on
%     ts_check_normals(g,normal_smooth)
% end

% Extrapolate wall normals, smooth at each step to build up field of
% normals away from the surface
% n_extrap = 20;
% xrrt_extrap = zeros(n_surf,n_extrap,3);
% xrrt_extrap(:,1,:) = reshape(xrrt_surf,[n_surf 1 3]);
% for n = 2:n_extrap
%     % Interpolate normals at extrapolated points
%     n_x = F_nx(reshape(xrrt_extrap(:,n-1,1),[],1,1), reshape(xrrt_extrap(:,n-1,2),[],1,1),...
%         reshape(xrrt_extrap(:,n-1,3),[],1,1));
%     n_r = F_nr(reshape(xrrt_extrap(:,n-1,1),[],1,1), reshape(xrrt_extrap(:,n-1,2),[],1,1),...
%         reshape(xrrt_extrap(:,n-1,3),[],1,1));
%     n_rt = F_nrt(reshape(xrrt_extrap(:,n-1,1),[],1,1), reshape(xrrt_extrap(:,n-1,2),[],1,1),...
%         reshape(xrrt_extrap(:,n-1,3),[],1,1));
%     
%     a = (n_x.^2 + n_r.^2 + n_rt.^2);
%     xrrt_extrap(:,n,:) = xrrt_extrap(:,n-1,:) + d_length/n_extrap *...
%         (cat(3,n_x,n_r,n_rt) ./ repmat(a,[1 1 3]));
%     
% end
% 
% % Plot positions of new points
% if plot_stuff == 1
%     ts_plot_surface(g, [0 0 0], figure())
%     plot3(reshape(xrrt_extrap(:,:,1),[],1),reshape(xrrt_extrap(:,:,3),[],1),reshape(xrrt_extrap(:,:,2),[],1),'r.')
% end
% 
% n_x = F_nx(xrrt_near(:,1),xrrt_near(:,2),xrrt_near(:,3));
% n_r = F_nr(xrrt_near(:,1),xrrt_near(:,2),xrrt_near(:,3));
% n_rt = F_nrt(xrrt_near(:,1),xrrt_near(:,2),xrrt_near(:,3));

% Plot positions of near point and positions of triangle centres to check
% coverage
if plot_stuff == 1
    ts_plot_surface(g, [0 0 0])
    plot3(xrrt_near(:,1),xrrt_near(:,3),xrrt_near(:,2),'r.');
    plot3(xrrt_fv(:,1),xrrt_fv(:,3),xrrt_fv(:,2),'b.');
end

% Get normals for near points by choosing n nearest triangles
d_xrrt_near = zeros(n_near,1); 
normal_near = zeros(n_near,3);

for n = 1:n_near
    m = abs(xrrt_near(n,1) - xrrt_fv(:,1)) < d_length/2 & ...
        abs(xrrt_near(n,2) - xrrt_fv(:,2)) < d_length/2 &...
        abs(xrrt_near(n,3) - xrrt_fv(:,3)) < d_length/2;
    if sum(m) > 0
        normal_near(n,:) = sum(normal_fv(m,:),1);
        d_xrrt_near(n) = sum(d_xrrt_fv(m)) / sum(m);
    end
    
%     d_xrrt_temp = d_xrrt_fv(i);  xrrt_temp = xrrt_fv(i,:); n_temp = sum(i);
%     normal_temp = normal_fv(i,:);
    
%     [~, m] = sort(abs(sum( (xrrt_fv(i,:) - repmat(xrrt_near(n,:),sum(i),1)).^2, 2)));
%     [~, m] = min(abs(sum( (xrrt_fv - repmat(xrrt_near(n,:),n_fv,1)).^2, 2)));
%     normal_near(n,:) = sum(normal_fv(m(1:n_tri),:),1);
%     d_xrrt_near = d_xrrt_fv(m(1:n_tri)) / n_tri;
%     i = double(i);
%     i(m(1:n_tri)) = 2;
%     normal_near(n,:) = sum(normal_fv(i == 2,:),1);
%     d_xrrt_near(n) = sum(d_xrrt_temp(i == 2)) / n_tri;
%     normal_near(n,:) = normal_fv(m,:);
%     d_xrrt_near(n) = d_xrrt_fv(m);
end

% Re normalise near normals
normal_near = normal_near ./ repmat(sum(normal_near.^2,2).^0.5,[1 3]);

% Extrapolate distance to shift, re-interpolate to near coordinates
% d_xrrt_extrap = repmat(d_xrrt_surf,[1 n_extrap]);

% Clear nans from extrapolated points
% xrrt_extrap = reshape(xrrt_extrap,[],3); 
% d_xrrt_extrap = reshape(d_xrrt_extrap,[],1);
% i = ~isnan(reshape(xrrt_extrap(:,1),[],1));
% xrrt_extrap = xrrt_extrap(i,:); d_xrrt_extrap = d_xrrt_extrap(i);
% F_d = TriScatteredInterp(xrrt_extrap(:,1),xrrt_extrap(:,2),...
%     xrrt_extrap(:,3),d_xrrt_extrap);
% 
% d_xrrt_near = F_d(xrrt_near(:,1),xrrt_near(:,2),xrrt_near(:,3));

% Surface normals
% n_x = F_nx(xrrt_surf(:,1),xrrt_surf(:,2),xrrt_surf(:,3));
% n_r = F_nr(xrrt_surf(:,1),xrrt_surf(:,2),xrrt_surf(:,3));
% n_rt = F_nrt(xrrt_surf(:,1),xrrt_surf(:,2),xrrt_surf(:,3));
% n_x(isnan(n_x)) = 0; n_r(isnan(n_r)) = 0; n_rt(isnan(n_rt)) = 0;
% 
% [n_x c] = smooth_IJK(bidIJK_surf,n_x,neighbour,10,0);
% [n_r c] = smooth_IJK(bidIJK_surf,n_r,neighbour,10,0,c);
% [n_rt c] = smooth_IJK(bidIJK_surf,n_rt,neighbour,10,0,c);

d_xrrt_surf = zeros(n_surf,1); 
normal_surf = zeros(n_surf,3);

for n = 1:n_surf
    m = abs(xrrt_surf(n,1) - xrrt_fv(:,1)) < d_length/2 & ...
        abs(xrrt_surf(n,2) - xrrt_fv(:,2)) < d_length/2 &...
        abs(xrrt_surf(n,3) - xrrt_fv(:,3)) < d_length/2;
    if sum(m) > 0
        normal_surf(n,:) = sum(normal_fv(m,:),1);
        d_xrrt_surf(n) = sum(d_xrrt_fv(m)) / sum(m);
    end
end

% Re normalise surface normals
normal_surf = normal_surf ./ repmat(sum(normal_surf.^2,2).^0.5,[1 3]);

% Calculate distance to shift on near points as a function of phi
phi_max = max(phi_near);
d_xrrt_near = d_xrrt_near .* cos((pi * phi_near) / (2 * phi_max));

% Nan removal from normals and shifts
d_xrrt_near(sum(isnan(normal_near),2) > 0) = 0;
normal_near(sum(isnan(normal_near),2) > 0,:) = 1;
d_xrrt_surf(sum(isnan(normal_surf),2) > 0) = 0;
normal_surf(sum(isnan(normal_surf),2) > 0,:) = 1;

%% Apply mode shape to mesh

% Find points that have significant movement
% d_xrrt_max = max(d_xrrt_near); d_xrrt_lim = 0;%d_xrrt_max / 10000;
% index_near = find(abs(d_xrrt_near) > d_xrrt_lim & pitch_near == 0);
% index_surf = find(abs(d_xrrt_surf) > d_xrrt_lim & pitch_surf == 0);

% xrrt_move = [xrrt_near(index_near,:) ; xrrt_surf(index_surf,:)];
% normal_move = [normal_near(index_near,:) ; normal_surf(index_surf,:)];
% d_xrrt_move = [d_xrrt_near(index_near) ; d_xrrt_surf(index_surf)];
% bidIJK_move = [bidIJK_near(index_near,:) ; bidIJK_surf(index_surf,:)];
xrrt_move = [xrrt_near ; xrrt_surf];
normal_move = [normal_near ; normal_surf];
d_xrrt_move = [d_xrrt_near ; d_xrrt_surf];
bidIJK_move = [bidIJK_near ; bidIJK_surf];
n_move = length(d_xrrt_move);

% Loop over all moved points and apply d_xrrt * normal
g_mode = g; xrrt_mode = xrrt_move - normal_move .* repmat(d_xrrt_move,[1 3]);

for n = 1:n_move
    bid = bidIJK_move(n,1); 
    I = bidIJK_move(n,2); J = bidIJK_move(n,3); K = bidIJK_move(n,4);
    
    g_mode{bid+1}.x(I+1,J+1,K+1) = xrrt_mode(n,1);
    g_mode{bid+1}.r(I+1,J+1,K+1) = xrrt_mode(n,2);
    g_mode{bid+1}.rt(I+1,J+1,K+1) = xrrt_mode(n,3);
end

% Write output
filename = [file_prefix 'M0' num2str(100*M) '_t' num2str(flow_type) '_mode_'...
    mode_prefix num2str(dp) '_input.hdf5'];

delete([directory filename]);
ts_write_hdf5(g_mode, [directory filename]);
ts_write_xdmf(g_mode, strrep([directory filename],'.hdf5','.xdmf'));

end


function [q c] = smooth_q(x,r,rt,q,range,n_smooth,m,c)
% Take a list of nodes with 3d coordinates and smooth a function via a moving average
% approach in three dimensions

% Calculate connectivity between points
if exist('c','var') == 0
    if exist('m','var') == 0
        m = 24;
    end
    c = zeros(length(x),m);

    for n = 1:length(x)
        % Select nodes in a range
        index_1 = find(x > x(n) - range & x < x(n) + range &...
            r > r(n) - range & r < r(n) + range &...
            rt > rt(n) - range & rt < rt(n) + range);

        xrrt_1 = [x(index_1) r(index_1) rt(index_1)];

        % Calculate straight line distance and choose m closest
        d = sum((xrrt_1 - repmat([x(n) r(n) rt(n)],[length(index_1) 1])).^2,2);

        [~, index_2] = sort(d);
        
        index_3 = index_1(index_2(1:m));

        c(n,:) = index_3.';
    end
else
    m = size(c,2);
end

% Smooth result a number of times
if exist('n_smooth','var') == 0
    n_smooth = 12;
end
for n = 1:n_smooth
    q_connect = q(c);
    
    q = sum(q_connect,2) / m;
end

end

function [q c] = smooth_IJK(bidIJK,q,neighbour,n_smooth,fix_0,c)
% Take a list of nodes with bid, IJK data and smooth a function via a
% connectivity approach in three dimensions

% Create connectivity array for list of points
n_q = length(q); 
if exist('c','var') == 0
    c = cell(n_q,1);
    for n = 1:n_q
        c_temp = neighbour{bidIJK(n,1)+1}{bidIJK(n,2)+1,bidIJK(n,3)+1,bidIJK(n,4)+1};

        c{n} = [];
        for m = 1:size(c_temp,1)
            c{n} = [c{n} ; find(sum(repmat(c_temp(m,:),[n_q 1]) == bidIJK,2) == 4)];
        end
        c{n} = unique(c{n});
    end
end

% Record points with zero intensity
if fix_0 == 1
    q_0 = q == 0;
end

% Loop over all points in list and average all neighbours as defined by c
for n = 1:n_smooth
    q_new = q;
    
    for m = 1:n_q
        q_new(m) = mean(q(c{m}));
    end

    % Reset zero intensities
    if fix_0 == 1
        q_new(q_0) = 0;
    end
    
    q = q_new;
end


end
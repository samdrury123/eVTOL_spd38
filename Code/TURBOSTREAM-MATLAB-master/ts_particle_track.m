function paths = ts_particle_track(g,seeds,xr_lim,forwards,cell_fraction,n_max,normal,order)
% Track particles forwards or backwards through the domain and return the
% paths

if exist('forwards','var') == 0
    forwards = 1;
end
if exist('cell_fraction','var') == 0
    cell_fraction = 0.8;
end
if exist('n_max','var') == 0
    n_max = 1500;
end
if exist('order','var') == 0
    order = 1;
end

% Limit bids range in format [x_min r_min ; x_max r_max];
x_min = xr_lim(1,1); x_max = xr_lim(2,1); 
r_min = xr_lim(1,2); r_max = xr_lim(2,2);
bids = [];
for bid = ts_get_bids(g)
    if min(min(min(g{bid+1}.x))) < x_max &&...
        max(max(max(g{bid+1}.x))) > x_min &&...
        min(min(min(g{bid+1}.r))) < r_max &&...
        max(max(max(g{bid+1}.r))) > r_min
        bids = [bids bid];
    end
end

% Align wall velocities with wall surfaces
for bid = bids
    Q = g{bid+1}.mwall == 0;

    n_x = normal{bid+1}(:,:,:,1); n_xq = n_x(Q);
    n_r = normal{bid+1}(:,:,:,2); n_rq = n_r(Q);
    n_rt = normal{bid+1}(:,:,:,3); n_rtq = n_rt(Q);

    Vx = g{bid+1}.Vx(Q); Vr = g{bid+1}.Vr(Q); Vt = g{bid+1}.Vt_rel(Q);

    mag = Vx .* n_xq + Vr .* n_rq + Vt .* n_rtq;

    g{bid+1}.Vx(Q) = Vx - mag.*n_xq;
    g{bid+1}.Vr(Q) = Vr - mag.*n_rq;
    g{bid+1}.Vt_rel(Q) = Vt - mag.*n_rtq;

    % Align wall velocities one node away from surface with surface
    [i,j,k] = ind2sub(size(g{bid+1}.mwall), find(g{bid+1}.mwall == 0));
    i = i-1; j = j-1; k = k-1;
    v0 = zeros(6,3); v1 = repmat(size(g{bid+1}.mwall)-1,[6 1]); 
    for n = 1:length(i)
        v = [i(n)+1 j(n) k(n) ; i(n)-1 j(n) k(n) ;...
        i(n) j(n)+1  k(n) ; i(n) j(n)-1  k(n) ;...
        i(n) j(n) k(n)+1 ; i(n) j(n) k(n)-1];
        v(sum(v < v0 | v > v1,2) > 0,:) = [];
        for m = 1:size(v,1)
            if g{bid+1}.mwall(v(m,1)+1,v(m,2)+1,v(m,3)+1) == 1
                Vx = g{bid+1}.Vx(v(m,1)+1,v(m,2)+1,v(m,3)+1); 
                Vr = g{bid+1}.Vr(v(m,1)+1,v(m,2)+1,v(m,3)+1); 
                Vt = g{bid+1}.Vt_rel(v(m,1)+1,v(m,2)+1,v(m,3)+1);

                mag = Vx .* n_x(i(n)+1,j(n)+1,k(n)+1) + Vr .* n_r(i(n)+1,j(n)+1,k(n)+1)...
                    + Vt .* n_rt(i(n)+1,j(n)+1,k(n)+1);

                g{bid+1}.Vx(v(m,1)+1,v(m,2)+1,v(m,3)+1) = Vx - mag * n_x(i(n)+1,j(n)+1,k(n)+1);
                g{bid+1}.Vr(v(m,1)+1,v(m,2)+1,v(m,3)+1) = Vr - mag * n_r(i(n)+1,j(n)+1,k(n)+1);
                g{bid+1}.Vt_rel(v(m,1)+1,v(m,2)+1,v(m,3)+1) = Vt - mag * n_rt(i(n)+1,j(n)+1,k(n)+1);
            end
        end
    end
end

% Pull out all cell corner coordinates along with velocities
varnames = {'x' 'r' 'rt' 'Vx' 'Vr' 'Vt'};
n_u = zeros(length(bids),1);
for n = 1:length(bids)
    n_u(n) = (g{bids(n)+1}.attribute.ni-1) * (g{bids(n)+1}.attribute.nj-1) *...
        (g{bids(n)+1}.attribute.nk-1);
end

u = zeros(sum(n_u),8,length(varnames)); n_u = [0 ; cumsum(n_u)];
for n = 1:length(bids)
    bid = bids(n);
    for v = 1:length(varnames)
        varname = varnames{v};
        
        u(n_u(n)+1:n_u(n+1),:,v) = reshape(cat(4,g{bid+1}.(varname)(1:end-1,1:end-1,1:end-1),...
            g{bid+1}.(varname)(2:end,1:end-1,1:end-1),...
            g{bid+1}.(varname)(1:end-1,2:end,1:end-1),...
            g{bid+1}.(varname)(1:end-1,1:end-1,2:end),...
            g{bid+1}.(varname)(2:end,2:end,1:end-1),...
            g{bid+1}.(varname)(2:end,1:end-1,2:end),...
            g{bid+1}.(varname)(1:end-1,2:end,2:end),...
            g{bid+1}.(varname)(2:end,2:end,2:end)),[],8);
    end
    nblade = double(g{bid+1}.bv.nblade);
end

% Apply limits in x and r to individual cells
r_lim = u(:,1,2) > r_min & u(:,1,2) < r_max;
x_lim = u(:,1,1) > x_min & u(:,1,1) < x_max;
u((r_lim & x_lim) == 0,:,:) = [];

% Include data either side of passage boundaries by pitch shifting some
% coordinates either side
n_u = size(u,1);
rt_mean = sum(u(:,1,3))/n_u;
q_1 = u(:,1,3) < rt_mean; q_2 = u(:,1,3) > rt_mean;
n_1 = sum(q_1); n_2 = sum(q_2);
u = cat(1,u,u(q_1,:,:),u(q_2,:,:));

u(n_u+1:n_u+n_1,:,3) = u(n_u+1:n_u+n_1,:,3) + u(n_u+1:n_u+n_1,:,2)*2*pi/nblade;
u(n_u+n_1+1:n_u+n_1+n_2,:,3) = u(n_u+n_1+1:n_u+n_1+n_2,:,3) - u(n_u+n_1+1:n_u+n_1+n_2,:,2)*2*pi/nblade;

% u = repmat(u,[3 1 1]);
% u(n_u+1:2*n_u,:,3) = u(n_u+1:2*n_u,:,3) + u(n_u+1:2*n_u,:,2)*2*pi/nblade;
% u(2*n_u+1:3*n_u,:,3) = u(2*n_u+1:3*n_u,:,3) - u(2*n_u+1:3*n_u,:,2)*2*pi/nblade;

% Pitch data
c = zeros(size(u,1),1,5);
c(:,1,5) = [zeros(n_u,1) ; ones(n_1,1) ; -ones(n_2,1)];
n_u = size(u,1);

% Calculate cell centred coordinates
c(:,1,1) = sum(u(:,:,1),2) / 8;
c(:,1,2) = sum(u(:,:,2),2) / 8;
c(:,1,3) = sum(u(:,:,3),2) / 8;

% Calculate radius based on centre that encloses all nodes - pythagoras on
% differences between cell centre coord and cell node coords
c(:,1,4) = max( (u(:,:,1) - repmat(c(:,1,1),[1 8])).^2  + ...
    (u(:,:,2) - repmat(c(:,1,2),[1 8])).^2  + ...
    (u(:,:,3) - repmat(c(:,1,3),[1 8])).^2  ,[],2);

% Calculate all 6 face normals
xrrt = u(:,:,1:3); VxVrVt = u(:,:,4:6); u(:,:,4:6) = [];
% Diagonal vectors for cross product
a = cat(2, xrrt(:,1,:) - xrrt(:,5,:),...
    xrrt(:,1,:) - xrrt(:,7,:),...
    xrrt(:,3,:) - xrrt(:,8,:),...
    xrrt(:,1,:) - xrrt(:,6,:),...
    xrrt(:,2,:) - xrrt(:,8,:),...
    xrrt(:,4,:) - xrrt(:,8,:));
b = cat(2, xrrt(:,2,:) - xrrt(:,3,:),...
    xrrt(:,4,:) - xrrt(:,3,:),...
    xrrt(:,7,:) - xrrt(:,5,:),...
    xrrt(:,4,:) - xrrt(:,2,:),...
    xrrt(:,5,:) - xrrt(:,6,:),...
    xrrt(:,6,:) - xrrt(:,7,:));

% Face centres
f = zeros(n_u,6,9);
f(:,:,4:6) = cat(2, xrrt(:,2,:) + xrrt(:,3,:) + xrrt(:,1,:) + xrrt(:,5,:),...
    xrrt(:,4,:) + xrrt(:,3,:) + xrrt(:,1,:) + xrrt(:,7,:),...
    xrrt(:,7,:) + xrrt(:,5,:) + xrrt(:,3,:) + xrrt(:,8,:),...
    xrrt(:,4,:) + xrrt(:,2,:) + xrrt(:,1,:) + xrrt(:,6,:),...
    xrrt(:,5,:) + xrrt(:,6,:) + xrrt(:,2,:) + xrrt(:,8,:),...
    xrrt(:,6,:) + xrrt(:,7,:) + xrrt(:,4,:) + xrrt(:,8,:)) / 4;

% Point normals to cell centre by dotting with cell centre to face centre
% vectors
f(:,:,1:3) = -cross(a,b) .* sign(repmat(dot(cross(a,b),...
    f(:,:,4:6) - repmat(c(:,:,1:3),[1 6 1]),3),[1 1 3]));

% Face centred velocities
f(:,:,7:9) = cat(2, VxVrVt(:,2,:) + VxVrVt(:,3,:) + VxVrVt(:,1,:) + VxVrVt(:,5,:),...
    VxVrVt(:,4,:) + VxVrVt(:,3,:) + VxVrVt(:,1,:) + VxVrVt(:,7,:),...
    VxVrVt(:,7,:) + VxVrVt(:,5,:) + VxVrVt(:,3,:) + VxVrVt(:,8,:),...
    VxVrVt(:,4,:) + VxVrVt(:,2,:) + VxVrVt(:,1,:) + VxVrVt(:,6,:),...
    VxVrVt(:,5,:) + VxVrVt(:,6,:) + VxVrVt(:,2,:) + VxVrVt(:,8,:),...
    VxVrVt(:,6,:) + VxVrVt(:,7,:) + VxVrVt(:,4,:) + VxVrVt(:,8,:)) / 4;
clear VxVrVt

% Calculate representative edge vectors
dl = zeros(n_u,3,3);
dl(:,:,1) = squeeze(xrrt(:,2,:) - xrrt(:,1,:));
dl(:,:,2) = squeeze(xrrt(:,3,:) - xrrt(:,1,:));
dl(:,:,3) = squeeze(xrrt(:,4,:) - xrrt(:,1,:));
clear xrrt

% And their magnitudes
mdl = (sum(dl.^2,2).^0.5);

%% Track particles

% Function to track particles through a passage from an inlet plane seed
% co-ordinate - use euler integration method and adaptive timestepping based on
% from Darmofal and Haimes
% Input seeds as a list of [X R RT]

% Initialise path matrix
paths = nan(size(seeds,1),4,n_max+10); paths(:,1:3,1) = seeds;
paths(:,4,1) = 1:size(paths,1); finished_paths = paths;

n_rt_total = size(seeds,1);
d = 5e-4 * ones(n_rt_total,1);
if forwards == 1
    x_out = x_max;
    x_in = x_min;
else
    x_out = x_min;
    x_in = x_max; 
end

% Recall arrays with node coordinates, cell coordinates, face data and
% edge data
% Nodes - u - [x r rt]
% Faces - f - [n_x n_r n_rt x r rt Vx Vr Vt]
% Centres - c - [x r rt radius pitch]
% Lines - dl - [dl1 dl2 dl3] & mdl - [mdl1 mdl2 mdl3]

z = ones(size(paths,1),1);

for n = 1:n_max
    % Get coordinates for unfinished paths
    x = paths(:,1,n); r = paths(:,2,n); rt = paths(:,3,n);
    disp([n sum(z) (max(x) - x_in) / (x_out - x_in) (min(x) - x_in) / (x_out - x_in)])
    
    n_rt = sum(z);
    
    % q is a matrix describing narrowing down of range to find bounding
    % cells
    
    % Find nodes within largest bounded dimensions of current points
    X = c(:,:,1); R = c(:,:,2); RT = c(:,:,3);
    qX = X<(max(x)+max(d)) & X>(min(x)-max(d));
    qR = R<(max(r)+max(d)) & R>(min(r)-max(d));
    qRT = RT<(max(rt)+max(d)) & RT>(min(rt)-max(d));
    
    q0 = qX & qR & qRT;

    % Initialise matrices to store cell data
    f_in = nan(n_rt,6,9); p_in = zeros(n_rt,1); 
    dl_in = nan(n_rt,3,3); mdl_in = nan(n_rt,1,3);
    
    % Find which cell it lies within, first tolerancing on an individual point, then
    % by discriminating on whether it is in the sphere described by a cell's 
    % maximum radius and finally by checking the face normals are all positive
    f_1 = f(q0,:,:); p_1 = c(q0,:,5); dl_1 = dl(q0,:,:); mdl_1 = mdl(q0,:,:);
    rad_2_1 = c(q0,:,4);
    
    X_c = X(q0); R_c = R(q0); RT_c = RT(q0);
    X_f = f_1(:,:,4); R_f = f_1(:,:,5); RT_f = f_1(:,:,6); 
    nX_f = f_1(:,:,1); nR_f = f_1(:,:,2); nRT_f = f_1(:,:,3); 

    for m = 1:n_rt
        % Tolerance in 3 dimensions based on a single coordinate
%         qX = X_c<(x(m)+d(m)) & X_c>(x(m)-d(m));
%         qR = R_c<(r(m)+d(m)) & R_c>(r(m)-d(m));
%         qRT = RT_c<(rt(m)+d(m)) & RT_c>(rt(m)-d(m));
%         
%         q2 = qX & qR & qRT;
%         
        % Find cell centred spheres that a point can exist within
%         dx = x(m) - X_c(q1); dr = r(m) - R_c(q1); drt = rt(m) - RT_c(q1);
%         q2 = (dx.^2 + dr.^2 + drt.^2) < rad_2_1(q1);
        dx = x(m) - X_c; dr = r(m) - R_c; drt = rt(m) - RT_c;
        q2 = (dx.^2 + dr.^2 + drt.^2) < rad_2_1;

        % Resolve current point into coords of face normals and check all
        % are positive for in cell flag
%         X_f_2 = X_f(q1,:); R_f_2 = R_f(q1,:); RT_f_2 = RT_f(q1,:); 
%         nX_f_2 = nX_f(q1,:); nR_f_2 = nR_f(q1,:); nRT_f_2 = nRT_f(q1,:); 
        dx = X_f(q2,:) - x(m); dr = R_f(q2,:) - r(m); drt = RT_f(q2,:) - rt(m);
        q3 = sum(dx.*nX_f(q2,:) + dr.*nR_f(q2,:) + drt.*nRT_f(q2,:) < 0,2) == 6;
%         dx = X_f_2(q2,:) - x(m); dr = R_f_2(q2,:) - r(m); drt = RT_f_2(q2,:) - rt(m);
%         q3 = sum(dx.*nX_f_2(q2,:) + dr.*nR_f_2(q2,:) + drt.*nRT_f_2(q2,:) < 0,2) == 6;

        % Translate logical arrays back up to narrow down whole domain to
        % one cell
%         q2(q2 == 1) = q3; q1(q1 == 1) = q2;
        q2(q2 == 1) = q3;
        
        % Record information about bounding cell
%         if sum(q1) == 1
%             f_in(m,:,:) = f_1(q1,:,:); p_in(m,:,:) = p_1(q1,:);
%             dl_in(m,:,:) = dl_1(q1,:,:); mdl_in(m,:,:) = mdl_1(q1,:,:);
%         end
        if sum(q2) == 1
            f_in(m,:,:) = f_1(q2,:,:); p_in(m,:,:) = p_1(q2,:);
            dl_in(m,:,:) = dl_1(q2,:,:); mdl_in(m,:,:) = mdl_1(q2,:,:);
        end    
    end
    
    % Interpolate velocity in cell with barycentric coodinates weighted on
    % volume of square based pyramid and face centred values of velocity.
    dx = repmat(x,[1 6]) - f_in(:,:,4);
    dr = repmat(r,[1 6]) - f_in(:,:,5);
    drt = repmat(rt,[1 6]) - f_in(:,:,6);
    
    Vol = (dx.*f_in(:,:,1) + dr.*f_in(:,:,2) + drt.*f_in(:,:,3)) / 3;
    
    % Calculate weights based on opposite pyramid
    W = flipdim(Vol,2) ./ repmat(sum(Vol,2),[1 6]);

    % Interpolate velocity in cell with tri linear interpolation adapted
    % for non regular cells
    
    
    % Velocities at coordinates
    Vx = sum(W.*f_in(:,:,7),2); 
    Vr = sum(W.*f_in(:,:,8),2); 
    Vt = sum(W.*f_in(:,:,9),2);

    % Pitch shift if out of bounds
    rt = rt - r.*p_in*2*pi/nblade;

    % Calculate time for a throughflow in each direction
    dt1 = mdl_in(:,1,1).^2 ./ abs((dl_in(:,1,1) .* Vx + dl_in(:,2,1) .* Vr + dl_in(:,3,1) .* Vt));
    dt2 = mdl_in(:,1,2).^2 ./ abs((dl_in(:,1,2) .* Vx + dl_in(:,2,2) .* Vr + dl_in(:,3,2) .* Vt));
    dt3 = mdl_in(:,1,3).^2 ./ abs((dl_in(:,1,3) .* Vx + dl_in(:,2,3) .* Vr + dl_in(:,3,3) .* Vt));
    
    % Set timestep and next search length    
    DT = cell_fraction * min(cat(2,dt1,dt2,dt3),[],2);
    d = 3 * max(mdl_in,[],3);

    % Find new particle position by Euler intergration
    if order == 1
        if forwards == 1
            paths(:,1,n+1) = x + DT.*Vx;
            paths(:,2,n+1) = r + DT.*Vr;
            paths(:,3,n+1) = rt + DT.*Vt;
        else
            paths(:,1,n+1) = x - DT.*Vx;
            paths(:,2,n+1) = r - DT.*Vr;
            paths(:,3,n+1) = rt - DT.*Vt;      
        end
    elseif order == 2 && n ~= 1
        if forwards == 1
            paths(:,1,n+1) = x + DT.*(2*Vx - Vx_old);
            paths(:,2,n+1) = r + DT.*(2*Vr - Vr_old);
            paths(:,3,n+1) = rt + DT.*(2*Vt - Vt_old);
        else
            paths(:,1,n+1) = x - DT.*(2*Vx - Vx_old);
            paths(:,2,n+1) = r - DT.*(2*Vr - Vr_old);
            paths(:,3,n+1) = rt - DT.*(2*Vt - Vt_old);      
        end
        Vx_old = Vx; Vr_old = Vr; Vt_old = Vt;
    else
        if forwards == 1
            paths(:,1,n+1) = x + DT.*Vx;
            paths(:,2,n+1) = r + DT.*Vr;
            paths(:,3,n+1) = rt + DT.*Vt;
        else
            paths(:,1,n+1) = x - DT.*Vx;
            paths(:,2,n+1) = r - DT.*Vr;
            paths(:,3,n+1) = rt - DT.*Vt;      
        end
        Vx_old = Vx; Vr_old = Vr; Vt_old = Vt;
    end
    paths(:,4,n+1) = DT;

    % z is vector describing paths still running
    x = paths(:,1,n+1);
    if forwards == 1
        z = (isnan(x) == 0) & (x < x_out);
    else
        z = (isnan(x) == 0) & (x > x_out);
    end
    
    % Transfer finished paths into finished matrix
    m_finish = paths(z == 0,4,1);
    finished_paths(m_finish,:,:) = paths(z == 0,:,:);
    paths(z == 0,:,:) = [];
    
    % Remove old velocities of finished paths
    if order == 2
        Vx_old(z == 0) = [];
        Vr_old(z == 0) = [];
        Vt_old(z == 0) = [];
    end
    
    % Check all paths have finished
    if size(paths,1) == 0
        break
    end
    
end
% Assign in unfinished paths anyway and re order
m = paths(:,4,1);
finished_paths(m,:,:) = paths;

paths = finished_paths;
n = sort(paths(:,4,1));
paths = paths(n,:,:);

end
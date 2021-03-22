function g = ts_add_seal(g,g_seal,bids_pass,bids_seal,plot_stuff)
% Connect two meshes together at the seal and passage interface

% Default to plot the mesh and check the patches
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

% Change blade numbers in the seal
for bid = ts_get_bids(g_seal)
    g_seal{bid+1}.bv.nblade = g{bids_pass(1)+1}.bv.nblade;
    g_seal{bid+1}.bv.fblade = g{bids_pass(1)+1}.bv.fblade;
end

% Get theta distributions from main passage blocks
x = []; r = []; t = [];
for bid = bids_pass
    x = [x ; squeeze(g{bid+1}.x(:,1,:))];
    r = [r ; squeeze(g{bid+1}.r(:,1,:))];
    t = [t ; squeeze(g{bid+1}.rt(:,1,:) ./ g{bid+1}.r(:,1,:))];
end

% Extend theta distributions
nk = g{bid+1}.attribute.nk;
x = [ones(1,nk)*-10 ; x ; ones(1,nk)*10]; r = [r(1,:) ; r ; r(end,:)]; t = [t(1,:) ; t ; t(end,:)];

% Loop over all k indices and write new tangential coordinates into the seal
for bid = ts_get_bids(g_seal)

    % Extract meridional coordinates
    g_seal{bid+1}.x = repmat(g_seal{bid+1}.x(:,:,1),[1 1 nk]);
    g_seal{bid+1}.r = repmat(g_seal{bid+1}.r(:,:,1),[1 1 nk]);
    g_seal{bid+1}.rt = zeros(g_seal{bid+1}.attribute.ni,g_seal{bid+1}.attribute.nj,nk);

    % Loop over all tangential indices and interpolate
    for k = 1:nk
        t_temp = interp1(x(:,k),t(:,k),g_seal{bid+1}.x(:,:,k),'linear');
        g_seal{bid+1}.rt(:,:,k) = t_temp .* g_seal{bid+1}.r(:,:,k);
    end

    % Update kst and ken in patch attributes
    for p = 1:length(g_seal{bid+1}.patch)
        if g_seal{bid+1}.patch{p}.attribute.kst > 0
            g_seal{bid+1}.patch{p}.attribute.kst = nk-1;
        end
        if g_seal{bid+1}.patch{p}.attribute.ken > 1
            g_seal{bid+1}.patch{p}.attribute.ken = nk;
        end
    end

    % Update block attributes
    g_seal{bid+1}.attribute.nk = nk;
end

% Morph radial coordinates of seal blocks to lie on the hub line
for n = 1:length(bids_pass)
    
    % Preallocate shift array
    dr = zeros(size(g_seal{bids_seal(n)+1}.x));
    
    % Loop over all tangential indices
    for k = 1:g_seal{bids_seal(n)+1}.attribute.nk
    
        % Maximum shift at hub line
        dr(:,end,k) = interp1(x(:,k),r(:,k),g_seal{bids_seal(n)+1}.x(:,end,k))...
            - g_seal{bids_seal(n)+1}.r(:,end,k);

        % Interpolate shift down to the bottom of the block
        for i = 1:g_seal{bids_seal(n)+1}.attribute.ni
            dr(i,:,k) = interp1(g_seal{bids_seal(n)+1}.r(i,[1 end],k),dr(i,[1 end],k),...
                g_seal{bids_seal(n)+1}.r(i,:,k));
        end
    end
    
    % Apply the radial shift to the whole block
    t = g_seal{bids_seal(n)+1}.rt ./ g_seal{bids_seal(n)+1}.r;
    g_seal{bids_seal(n)+1}.r = g_seal{bids_seal(n)+1}.r + dr;
    g_seal{bids_seal(n)+1}.rt = t .* g_seal{bids_seal(n)+1}.r;
end

% Get axial coordinates from seals to apply to main passage
for n = 1:length(bids_pass)

    % Get meridional coordinates of seal to join to main passage
    x_seal = g_seal{bids_seal(n)+1}.x(:,end,1);
    r_seal = g_seal{bids_seal(n)+1}.r(:,end,1);
    rt_seal = g_seal{bids_seal(n)+1}.rt(:,end,:);

    % Get index of main passage join to seal
    [~,i_start] = min(abs(g{bids_pass(n)+1}.x(:,1,1) - x_seal(1)));
    [~,i_end] = min(abs(g{bids_pass(n)+1}.x(:,1,1) - x_seal(end)));

    % Get coordinates of main passage block above seal
    x = g{bids_pass(n)+1}.x(i_start:i_end,:,:);
    r = g{bids_pass(n)+1}.r(i_start:i_end,:,:);
    rt = g{bids_pass(n)+1}.rt(i_start:i_end,:,:);
    q = linspace(0,1,i_end-i_start+1);
    qi = linspace(0,1,length(x_seal));

    % Reinterpolate mesh to right number of i-coords
    [~,nj,nk] = size(x); ni = length(qi);
    x_new = zeros(ni,nj,nk); r_new = zeros(ni,nj,nk); rt_new = zeros(ni,nj,nk);
    for j = 1:nj
        for k = 1:nk
            x_new(:,j,k) = interp1(q,x(:,j,k),qi,'linear');
            r_new(:,j,k) = interp1(q,r(:,j,k),qi,'linear');
            rt_new(:,j,k) = interp1(q,rt(:,j,k),qi,'linear');
        end
    end

    % Calculate shift in coordinates of new mesh
    dx = x_new(:,1,:) - repmat(reshape(x_seal,[],1),[1 1 nk]);
    dr = r_new(:,1,:) - repmat(reshape(r_seal,[],1),[1 1 nk]);
    drt = rt_new(:,1,:) - reshape(rt_seal,[],1,nk);

    % Assemble shift matrices with decaying cos shift
    nj_shift = 10;
    q = [cos(pi*(1:nj_shift)/(2*nj_shift)) zeros(1,nj-nj_shift)];
    dx = repmat(dx,[1 nj 1]) .* repmat(q,[ni 1 nk]);
    dr = repmat(dr,[1 nj 1]) .* repmat(q,[ni 1 nk]);
    drt = repmat(drt,[1 nj 1]) .* repmat(q,[ni 1 nk]);

    % Apply shifts
    x_new = x_new - dx;
    r_new = r_new - dr;
    rt_new = rt_new - drt;

    % Reassemble main passage block with right number of points and matching mesh
    if i_end < size(g{bids_pass(n)+1}.x,1)
        i_rear = i_end+1:size(g{bids_pass(n)+1}.x,1);
    else
        i_rear = [];
    end
    if i_start ~= 1
        i_front = 1:i_start-1;
    else
        i_front = [];
    end            
    g{bids_pass(n)+1}.x = [g{bids_pass(n)+1}.x(i_front,:,:) ; x_new ; ...
        g{bids_pass(n)+1}.x(i_rear,:,:)];
    g{bids_pass(n)+1}.r = [g{bids_pass(n)+1}.r(i_front,:,:) ; r_new ; ...
        g{bids_pass(n)+1}.r(i_rear,:,:)];
    g{bids_pass(n)+1}.rt = [g{bids_pass(n)+1}.rt(i_front,:,:) ; rt_new ; ...
        g{bids_pass(n)+1}.rt(i_rear,:,:)];
    
    % Correct rear indices
    if isempty(i_rear) == 0
        i_rear = i_rear + size(x_new,1) - size(x,1);
    end
    
    % Apply shifts to coordinate up and downstream
    if isempty(i_front) == 0
        f = repmat(linspace(0,1,length(i_front))',[1 nj nk]);
        g{bids_pass(n)+1}.x(i_front,:,:) = g{bids_pass(n)+1}.x(i_front,:,:) - ...
            repmat(dx(1,:,:),[length(i_front) 1 1]) .* f;
        g{bids_pass(n)+1}.r(i_front,:,:) = g{bids_pass(n)+1}.r(i_front,:,:) - ...
            repmat(dr(1,:,:),[length(i_front) 1 1]) .* f;
        g{bids_pass(n)+1}.rt(i_front,:,:) = g{bids_pass(n)+1}.rt(i_front,:,:) - ...
            repmat(drt(1,:,:),[length(i_front) 1 1]) .* f;
    end
    if isempty(i_rear) == 0
        f = repmat(linspace(1,0,length(i_rear))',[1 nj nk]);
        g{bids_pass(n)+1}.x(i_rear,:,:) = g{bids_pass(n)+1}.x(i_rear,:,:) - ...
            repmat(dx(end,:,:),[length(i_rear) 1 1]) .* f;
        g{bids_pass(n)+1}.r(i_rear,:,:) = g{bids_pass(n)+1}.r(i_rear,:,:) - ...
            repmat(dr(end,:,:),[length(i_rear) 1 1]) .* f;
        g{bids_pass(n)+1}.rt(i_rear,:,:) = g{bids_pass(n)+1}.rt(i_rear,:,:) - ...
            repmat(drt(end,:,:),[length(i_rear) 1 1]) .* f;
    end
    
    % Update block attributes
    ni = size(g{bids_pass(n)+1}.x,1);
    g{bids_pass(n)+1}.attribute.ni = ni;

    % Update ist & ien in patch attributes
    for pid = ts_get_pids(g,bids_pass(n))
        if g{bids_pass(n)+1}.patch{pid+1}.attribute.ist > 0
            g{bids_pass(n)+1}.patch{pid+1}.attribute.ist = ni-1;
        end
        if g{bids_pass(n)+1}.patch{pid+1}.attribute.ien > 1
            g{bids_pass(n)+1}.patch{pid+1}.attribute.ien = ni;
        end
    end 
end

% Plot meridional mesh
if plot_stuff == 1
    
    % Open the figure window
    figure(); hold on; axis equal; box on; axis equal;
    
    % Plot the passage blocks in black
    for bid = ts_get_bids(g)
        mesh(g{bid+1}.x(:,:,1),g{bid+1}.r(:,:,1),g{bid+1}.rt(:,:,1),...
            zeros(size(g{bid+1}.x(:,:,1:3))));
    end
    
    % Plot the seal blocks in red
    for bid = ts_get_bids(g_seal)
        [ni,nj,~] = size(g_seal{bid+1}.x);
        mesh(g_seal{bid+1}.x(:,:,1),g_seal{bid+1}.r(:,:,1),g_seal{bid+1}.rt(:,:,1),...
            cat(3,ones(ni,nj),zeros(ni,nj,2)));
    end
end

% Update block numbers in seal attributes and patch attributes
bid_st = length(g);
for bid = 0:length(g_seal)-1
    g_seal{bid+1}.attribute.bid = g_seal{bid+1}.attribute.bid + bid_st;
    for pid = 0:length(g_seal{bid+1}.patch)-1
        g_seal{bid+1}.patch{pid+1}.attribute.bid = ...
            g_seal{bid+1}.patch{pid+1}.attribute.bid + bid_st;
        g_seal{bid+1}.patch{pid+1}.attribute.nxbid = ...
            g_seal{bid+1}.patch{pid+1}.attribute.nxbid + bid_st;
    end
end

% Copy flow from the main passage
varnames = {'x' 'ro' 'rorvt' 'roe' 'trans_dyn_vis' 'phi'}; 
guess = zeros(length(bids_pass),length(varnames));
for n = 1:length(bids_pass)
    for v = 1:length(varnames)
        guess(n,v) = mean(mean(g{bids_pass(n)+1}.(varnames{v})(:,1,:)));
    end
end

% Interpolate values in seal blocks
for n = 1:length(g_seal)
    
    % Set radial and axial velocities to zero
    g_seal{n}.rovx = zeros(size(g_seal{n}.x)); g_seal{n}.rovr = zeros(size(g_seal{n}.x)); 
    
    % Copy other values from main passage
    for v = 2:length(varnames)
        g_seal{n}.(varnames{v}) = interp1(guess(:,1),guess(:,v),g_seal{n}.x,'linear','extrap');
    end
end

% Append seal blocks to main passage
g = [g ; g_seal]; bids_seal = bids_seal + bid_st;

% Patch seals onto main passage
for n = 1:length(bids_pass)

    % Find matching indices on passage blocks
    [~, i1] = min(abs(g{bids_seal(n)+1}.x(1,end,1) - g{bids_pass(n)+1}.x(:,1,1)));
    [~, i2] = min(abs(g{bids_seal(n)+1}.x(end,end,1) - g{bids_pass(n)+1}.x(:,1,1)));

    % Define directions and indices for both blocks
    Q0 = [i1-1 i2 ; 0 1 ; 0 g{bids_pass(n)+1}.attribute.nk]; 
    dir0 = [0 -1 2];
    Q1 = [0 g{bids_seal(n)+1}.attribute.ni ; g{bids_seal(n)+1}.attribute.nj-1 ...
        g{bids_seal(n)+1}.attribute.nj ; 0 g{bids_seal(n)+1}.attribute.nk];
    dir1 = [0 -1 2];

    % Add the patches
    g = ts_add_periodic_patch(bids_pass(n), bids_seal(n), Q0, Q1, dir0, dir1, g, 5);
end

% Match periodics and patches
g = ts_match_periodics(g);
g = ts_match_patches(g);

% Set rpms on faces at seal inlet and outlet
for n = 1:length(bids_seal)
    if g{bids_pass(n)+1}.bv.rpm == 0
        g{bids_pass(n)+1}.bv.rpmj1 = 1;
        g{bids_seal(n)+1}.bv.rpmj2 = 1;
    else
        g{bids_pass(n)+1}.bv.rpmj1 = 0;
        g{bids_seal(n)+1}.bv.rpmj2 = 0;        
    end
end

% Plot rpms on faces and check the patches
if plot_stuff == 1
    ts_check_rpms(g); ts_check_patches(g,2);
end


end


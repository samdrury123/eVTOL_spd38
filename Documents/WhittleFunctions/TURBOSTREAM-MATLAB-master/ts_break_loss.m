function [b, h, xrrt] = ts_break_loss(g,inlet,outlet,break_type,use_entropy,plot_stuff)
% TS_BREAK_LOSS  Split loss into 2D and 3D components by tracking streamlines through a blade row
%
%   [b, h] = TS_BREAK_LOSS(g,inlet,outlet,use_entropy,plot_stuff)
%
%   g - cell array of TURBOSTREAM solution
%   inlet - data structure cut plane at inlet
%   outlet - data structure cut plane at outlet
%   break_type - data structure controlling loss breakdown
%   use_entropy - 0 or 1 for using stagnation pressure or entropy loss coefficient
%   plot_stuff - 0 or 1 for showing working
%   b - output data structure with loss breakdown and calculations
%   h - figure window handles
%   xrrt - streamline coordinates and plane intersections
%
%   break_type.type is 'mass', 'radius' or 'vorticity'
%   break_type.val defines threshold for loss breakdown


%% Open figure windows and prepare data

% Set plotting default
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Choose to plot non-uniform average breakdowns
plot_avg = 0;

% Open figure windows
if plot_stuff == 1
    h.domain = figure();
    h.wake = figure();
    h.bl = figure();
    
    figure(h.wake);
    t_mean = ts_plot_wake(inlet, outlet, subplot(1,1,1), use_entropy, 2, 1);
    
    figure(h.bl);
    subplot(1,2,1); hold on; grid on; box on; xlabel('Whirl Angle'); ylabel('Radius');
    subplot(1,2,2); hold on; grid on; box on; xlabel('Stagnation Pressure');
else
    h = [];
end

% Get secondary variables
inlet = ts_secondary(inlet); outlet = ts_secondary(outlet); g = ts_secondary(g);

% Work in double precision
property_names = fieldnames(inlet);
for i = 1:length(property_names)
    if isa(inlet.(property_names{i}), 'single') == 1
        inlet.(property_names{i}) = double(inlet.(property_names{i}));
    end
end

property_names = fieldnames(outlet);
for i = 1:length(property_names)
    if isa(outlet.(property_names{i}), 'single') == 1
        outlet.(property_names{i}) = double(outlet.(property_names{i}));
    end
end

for bid = ts_get_bids(g)
    property_names = fieldnames(g{bid+1});
    for i = 1:length(property_names)
        if isa(g{bid+1}.(property_names{i}), 'single') == 1
            g{bid+1}.(property_names{i}) = double(g{bid+1}.(property_names{i}));
        end
    end
end

% Determine which blocks lie in domain in question
bids = [];
for bid = ts_get_bids(g)
    if min(min(min(g{bid+1}.x))) < max(max(outlet.x)) && ...
           max(max(max(g{bid+1}.x))) > min(min(inlet.x)) 
       bids = [bids bid];
    end
end

% Write domain .hdf5 file
[~,~] = system('rm PV/domain.*');
ts_export_paraview(g,'PV/domain.hdf5','Track',bids,1);

% Plot domain surface
if plot_stuff == 1
    figure(h.domain);
    ts_plot_surface(g,[0 0 0],h.domain,1,1,bids);
end

%% Select radii to split inlet mesh

% Break inlet based on non-dimensional radius
if strcmp(break_type.type,'radius') == 1
    r_hub = min(min(inlet.r)); r_cas = max(max(inlet.r));
    r_split(1) = (r_cas - r_hub) * break_type.val + r_hub;
    r_split(2) = (r_cas - r_hub) * (1-break_type.val) + r_hub;
end

% Break inlet based on mass flow
if strcmp(break_type.type,'mass') == 1
    [r,mass] = ts_mass_average(inlet,'r',2);
    mass = [0 cumsum(mass)]; mass = mass / mass(end);
    r_split(1) = interp1(mass,r,break_type.val,'pchip');
    r_split(2) = interp1(mass,r,(1-break_type.val),'pchip');
end

% Break inlet based upon inlet boundary layer vorticity
if strcmp(break_type.type,'vorticity') == 1
    
    % Get inlet flow conditions
    Vx = ts_mass_average(inlet,'Vx',2);
    r = ts_mass_average(inlet,'r',2);

    % Calculate vorticity
    dVx_dr = diff(Vx) ./ diff(r);
    r = 0.5*(r(1:end-1) + r(2:end));
    
    % Find where vorticity crosses threshold
    j1 = find((dVx_dr(2:end) - break_type.val) .* (dVx_dr(1:end-1) - break_type.val) < 0,1,'first');
    j2 = find((dVx_dr(2:end) + break_type.val) .* (dVx_dr(1:end-1) - break_type.val) < 0,1,'last');

    r_split(1) = interp1(dVx_dr(j1-3:j1+3),r(j1-3:j1+3),break_type.val,'pchip');
    r_split(2) = interp1(dVx_dr(j2-3:j2+3),r(j2-3:j2+3),-break_type.val,'pchip');
end

% Plot split points on inlet profile
if plot_stuff == 1
    figure(h.bl);
    
    % Get flow properties
    r = ts_mass_average(inlet,'r',2); 
    Alpha = ts_mass_average(inlet,'Alpha',2);
    Po = ts_mass_average(inlet,'Po',2);
    
    % Plot inlet profiles
    subplot(1,2,1); plot(Alpha,r,'b-','LineWidth',2); v = axis; axis([v(1:2) min(r) max(r)]);
    subplot(1,2,2); plot(Po,r,'b-','LineWidth',2); v = axis; axis([v(1:2) min(r) max(r)]);
    
    % Plot split points
    subplot(1,2,1); plot(interp1(r,Alpha,r_split),r_split,'ko','MarkerFaceColor','k')
    subplot(1,2,2); plot(interp1(r,Po,r_split), r_split,'ko','MarkerFaceColor','k')
end

% Plot split points on downstream wake
if plot_stuff == 1
    figure(h.wake); v = axis;
    
    % Calculate cartesian coordinates
    t = linspace(-pi/4,pi/4,100); 
    for n = 1:length(r_split)
        [y,z] = pol2cart(t,r_split(n)); plot(z,y,'k-','LineWidth',2);
    end
    
    % Reset axis size
    axis(v);
end

% Record split radii
b.r_split = r_split;

%% Write initial coarse array of points across the passage and track in PARAVIEW

% Find pitchwise extents of inlet plane
rt_split_1 = interp1(inlet.r(:,1),inlet.rt(:,1),b.r_split,'linear');
rt_split_2 = interp1(inlet.r(:,end),inlet.rt(:,end),b.r_split,'linear');

% Axial coordinate at inlet plane
x_split = interp1(inlet.r(:,1),inlet.x(:,1),b.r_split,'linear');

% Linearly spaced points to track
nk = 100; xrrt_in_coarse = cell(length(r_split),1);
for n = 1:length(r_split)
    xrrt_in_coarse{n} = [ones(nk,1) * x_split(n), ones(nk,1) * r_split(n), ...
        linspace(rt_split_1(n),rt_split_2(n),nk)'];
end

% Assemble into cell arrays
g_seeds = cell(length(r_split),1);
for n = 1:length(r_split)
    g_seeds{n}.attribute.bid = n-1;
    g_seeds{n}.attribute.ni = nk; g_seeds{n}.attribute.nj = 1; g_seeds{n}.attribute.nk = 1; 
    g_seeds{n}.x = xrrt_in_coarse{n}(:,1); 
    g_seeds{n}.r = xrrt_in_coarse{n}(:,2); 
    g_seeds{n}.rt = xrrt_in_coarse{n}(:,3); 
end

% Write seed .hdf5 file
[~,~] = system('rm PV/seeds.*');
ts_export_paraview(g_seeds,'PV/seeds.hdf5',{},[],1);

% Plot seed points
if plot_stuff == 1
    figure(h.domain);
    for n = 1:length(g_seeds)
        [y,z] = pol2cart(g_seeds{n}.rt./g_seeds{n}.r,g_seeds{n}.r);
        plot3(g_seeds{n}.x,z,y,'b.');
    end
end

% Track solution in paraview
[~,~] = system('rm PV/lines*.csv');
cd('PV');
[~,~] = system('pvpython pv_track_seeds.py');
cd('..');

% Read line data back to matlab
xyz_coarse = cell(length(r_split),1);
for n = 1:length(r_split)
    s = pv_read_csv(['PV/lines' num2str(n-1) '.csv']); xyz_coarse{n} = s.xyz;
end

% Plot coarse tracked lines
if plot_stuff == 1
    figure(h.domain);
    for n = 1:length(r_split)
        plot3(xyz_coarse{n}(:,:,1),xyz_coarse{n}(:,:,3),xyz_coarse{n}(:,:,2),'b-')
    end
end

%% Find where the blade splits the streamsurfaces, refine points distribution and retrack

% Remove points that do not cross outlet plane
for n = 1:length(r_split)
    q = max(xyz_coarse{n}(:,:,1),[],1) < min(min(outlet.x));
    xyz_coarse{n}(:,q,:) = [];
end

% Find x at midchord where tangential difference is maximum
x = linspace(max(max(inlet.x)),min(min(outlet.x)),20);
n = 1; nk = size(xyz_coarse{1},2); z = zeros(nk,length(x)); 
for k = 1:nk
    q = isnan(xyz_coarse{n}(:,k,1)) == 0;
    x_temp = xyz_coarse{n}(q,k,1); z_temp = xyz_coarse{n}(q,k,3);
    [x_temp,i] = unique(x_temp); z_temp = z_temp(i);
    z(k,:) = interp1(x_temp,z_temp,x,'linear');
end

[~,i] = max(max(diff(z,1,1),[],1),[],2); x_cen = x(i);

% Loop over all sets of streamlines
xrrt_coarse = cell(size(xyz_coarse)); g_seeds = cell(length(r_split),1);
xrrt_in_fine = cell(length(r_split),1);
for n = 1:length(xyz_coarse)
    
    % Calculate polar streamline coordinates
    [t,r] = cart2pol(xyz_coarse{n}(:,:,2),xyz_coarse{n}(:,:,3));
    xrrt_coarse{n} = cat(3,xyz_coarse{n}(:,:,1),r,r.*t);
    
    % Interpolate tangential coordinates of streamlines at x-coordinate at mid-chord
    nk = size(xrrt_coarse{n},2); rt = zeros(nk,1); 
    for k = 1:nk
        q = isnan(xrrt_coarse{n}(:,k,1)) == 0;
        x_temp = xrrt_coarse{n}(q,k,1); rt_temp = xrrt_coarse{n}(q,k,3);
        [x_temp,i] = unique(x_temp); rt_temp = rt_temp(i);
        rt(k) = interp1(x_temp,rt_temp,x_cen,'linear');
    end
    
    % Find discontinuity in tangential coordinate for blade surface
    drt = abs(diff(rt));
    k = find(drt > 5 * mean(drt));
    
    % Select mid-pitch if no discontinuity is present
    if isempty(k) == 1
        k = round(nk/2);
    end
    
    % Define new range of seeds to track
    dk = 3; nk = 98;
    xrrt_in_fine{n} = [ones(nk,1) * x_split(n), ones(nk,1) * r_split(n), ...
        linspace(xrrt_coarse{n}(1,k-dk,3),xrrt_coarse{n}(1,k+dk,3),nk)'];
    
    % Assemble into a cell array
    g_seeds{n}.attribute.bid = n-1;
    g_seeds{n}.attribute.ni = nk; g_seeds{n}.attribute.nj = 1; g_seeds{n}.attribute.nk = 1; 
    g_seeds{n}.x = xrrt_in_fine{n}(:,1); 
    g_seeds{n}.r = xrrt_in_fine{n}(:,2); 
    g_seeds{n}.rt = xrrt_in_fine{n}(:,3); 

end

% Write seed .hdf5 file
[~,~] = system('rm PV/seeds.*');
ts_export_paraview(g_seeds,'PV/seeds.hdf5',{},[],1);

% Plot seed points
if plot_stuff == 1
    figure(h.domain);
    for n = 1:length(g_seeds)
        [y,z] = pol2cart(g_seeds{n}.rt./g_seeds{n}.r,g_seeds{n}.r);
        plot3(g_seeds{n}.x,z,y,'r.');
    end
end

% Track solution in paraview
[~,~] = system('rm PV/lines*.csv');
cd('PV');
[~,~] = system('pvpython pv_track_seeds.py');
cd('..');

% Read line data back to matlab
xyz_fine = cell(length(r_split),1);
for n = 1:length(r_split)
    s = pv_read_csv(['PV/lines' num2str(n-1) '.csv']); xyz_fine{n} = s.xyz;
end

% Plot fine tracked lines
if plot_stuff == 1
    figure(h.domain);
    for n = 1:length(r_split)
        plot3(xyz_fine{n}(:,:,1),xyz_fine{n}(:,:,3),xyz_fine{n}(:,:,2),'r-')
    end
end

% Assemble coarse and fine lines into one sorted array
xyz = cell(length(xyz_fine));
for n = 1:length(xyz_fine)
    
    % Pad shorter array with nans
    ni_coarse = size(xyz_coarse{n},1); ni_fine = size(xyz_fine{n},1);
    if ni_fine > ni_coarse
        xyz_coarse{n} = [xyz_coarse{n} ; nan(ni_fine - ni_coarse,size(xyz_coarse{n},2),3)];
    else
        xyz_fine{n} = [xyz_fine{n} ; nan(ni_coarse - ni_fine,size(xyz_fine{n},2),3)];
    end
    
    % Assemble and sort
    xyz{n} = [xyz_coarse{n} xyz_fine{n}];
    [~,k] = sort(xyz{n}(1,:,3));
    xyz{n} = xyz{n}(:,k,:);
end

% Plot order of lines on domain
if plot_stuff == 1
    figure(h.domain)
    for n = 1:length(xyz)
        plot3(xyz{n}(1,:,1),xyz{n}(1,:,3),xyz{n}(1,:,2),'g-');
    end
end

% Calculate polar streamline coordinates
xrrt_lines = cell(size(xyz));
for n = 1:length(xyz)
    [t,r] = cart2pol(xyz{n}(:,:,2),xyz{n}(:,:,3));
    xrrt_lines{n} = cat(3,xyz{n}(:,:,1),r,r.*t);
end

%% Compute intersections of streamlines with outlet plane

% Extract intersections of streamlines with outlet plane
xrrt_out = cell(length(xrrt_lines),1);
for n = 1:length(xrrt_lines)
    xrrt_out{n} = extract_line(xrrt_lines{n}, outlet);
end

% Remove out of range points
r_cas = max(outlet.r(:)); r_hub = min(outlet.r(:));
for n = 1:length(xrrt_out)
    xrrt_out{n}(xrrt_out{n}(:,2) > r_cas | xrrt_out{n}(:,2) < r_hub,:) = [];
end

% Plot tracked streamlines on wakes
if plot_stuff == 1
    figure(h.wake)
    for n = 1:length(xrrt_out)
        r = xrrt_out{n}(:,2);
        t = xrrt_out{n}(:,3) ./ r;
        plot((t - t_mean) .* r,r,'r-')
    end
end

% Record seed lines
xrrt_in = cell(length(xrrt_lines),1);
for n = 1:length(xrrt_lines)
    xrrt_in{n} = squeeze(xrrt_lines{n}(1,:,:));
end

% Extend seed lines +/- rt
for n = 1:length(xrrt_in)
    xrrt_in{n} = [xrrt_in{n}(1,:) ; xrrt_in{n} ; xrrt_in{n}(end,:)];
    drt = 0.5*abs(rt_split_2(n) - rt_split_1(n));
    xrrt_in{n}(1,3) = xrrt_in{n}(1,3) - drt;
    xrrt_in{n}(end,3) = xrrt_in{n}(end,3) + drt;
end

% Return coordinates of streamsurfaces
xrrt.xrrt = xrrt_lines; xrrt.out = xrrt_out; xrrt.in = xrrt_in;

%% Calculate loss generation in sectors divided by streamsurfaces

% Determine loss varaiable to use
if use_entropy == 1
    loss_varname = {'s','T'};
else
    loss_varname = 'Po';
end

% Preallocate arrays
loss_var_in = zeros(length(xrrt_lines)+1,length(loss_varname)); 
loss_var_out = zeros(length(xrrt_lines)+1,length(loss_varname)); 
m_in = zeros(length(xrrt_lines)+1,1); m_out = zeros(length(xrrt_lines)+1,1); 
A_in = zeros(length(xrrt_lines)+1,1); A_out = zeros(length(xrrt_lines)+1,1); 

% Compute hub sector
[loss_var_in(1,:),m_in(1),A_in(1)] = ...
    ts_mass_average_nonuniform(inlet,loss_varname,xrrt_in{1},1,plot_avg);
[loss_var_out(1,:),m_out(1),A_out(1)] = ...
    ts_mass_average_nonuniform(outlet,loss_varname,xrrt_out{1},1,plot_avg);
% [loss_var_in(1,:),m_in(1),A_in(1)] = ...
%     ts_mass_average_nonuniform_2(inlet,loss_varname,[0 0 0 ; xrrt_in{1} ; 0 0 0],1);
% [loss_var_out(1,:),m_out(1),A_out(1)] = ...
%     ts_mass_average_nonuniform_2(outlet,loss_varname,[0 0 0 ; xrrt_out{1} ; 0 0 0]);

% Compute mid sectors
for n = 1:length(xrrt_out)-1
    line_temp = [xrrt_out{n} ; flipud(xrrt_out{n+1})];
    [loss_var_out(n+1,:),m_out(n+1),A_out(n+1)] = ...
        ts_mass_average_nonuniform(outlet,loss_varname,line_temp,2,plot_avg);
    line_temp = [xrrt_in{n} ; flipud(xrrt_in{n+1})];
    [loss_var_in(n+1,:),m_in(n+1),A_in(n+1)] = ...
        ts_mass_average_nonuniform(inlet,loss_varname,line_temp,2,plot_avg);
%     line_temp = [xrrt_out{n} ; flipud(xrrt_out{n+1})];
%     [loss_var_out(n+1,:),m_out(n+1),A_out(n+1)] = ...
%         ts_mass_average_nonuniform_2(outlet,loss_varname,line_temp);
%     line_temp = [xrrt_in{n} ; flipud(xrrt_in{n+1})];
%     [loss_var_in(n+1,:),m_in(n+1),A_in(n+1)] = ...
%         ts_mass_average_nonuniform_2(inlet,loss_varname,line_temp);
end

% Compute casing sector
[loss_var_in(end,:),m_in(end),A_in(end)] = ...
    ts_mass_average_nonuniform(inlet,loss_varname,xrrt_in{end},3,plot_avg);
[loss_var_out(end,:),m_out(end),A_out(end)] = ...
    ts_mass_average_nonuniform(outlet,loss_varname,xrrt_out{end},3,plot_avg);
% [loss_var_in(end,:),m_in(end),A_in(end)] = ...
%     ts_mass_average_nonuniform_2(inlet,loss_varname,[1 1 1 ; xrrt_in{end} ; 1 1 1]);
% [loss_var_out(end,:),m_out(end),A_out(end)] = ...
%     ts_mass_average_nonuniform_2(outlet,loss_varname,[1 1 1 ; xrrt_out{end} ; 1 1 1]);

%% Calculate loss coefficients

% Entropy generation coefficient
if use_entropy == 1
    V_in = ts_mass_average(inlet,'V',3);
    Loss = 2 * loss_var_out(:,2) .* (loss_var_out(:,1) - loss_var_in(:,1)) / (V_in^2);
end

% Stagnation pressure loss coefficient
if use_entropy == 0
    Po_in = ts_mass_average(inlet,'Po',3);
    P_in = ts_mass_average(inlet,'P',3);
    Loss = (loss_var_in - loss_var_out) / (Po_in - P_in);
    b.Po_in = loss_var_in; b.Po_out = loss_var_out; b.P_ref = Po_in - P_in;
end
    

%% Clean files and return output

% Delete tracking files
[~,~] = system('rm PV/domain.*');
[~,~] = system('rm PV/seeds.*');
[~,~] = system('rm PV/lines*.csv');

% Assign results into data structure
b.Loss = Loss;
b.m_in = m_in; b.m_out = m_out; 
b.A_in = A_in; b.A_out = A_out;


end

function xrrt_line = extract_line(xrrt_paths, outlet)
% Calculate intersections of streamlines and an outlet plane

% Calculate intersections if outlet is not at constant x
if all(outlet.x == mean(outlet.x(:))) == 0

    % Permute paths array
    xrrt_paths = permute(xrrt_paths,[2 3 1]);

    % Re organise paths to intersect outlet plane - definining intersecting line
    xrrt_line = nan(size(xrrt_paths,1),3);
    for i = 1:size(xrrt_paths,1)
        n1 = sum(~isnan(xrrt_paths(i,1,:))); n2 = length(outlet.x(:,1));

        % Line coordinates and differences in #1
        xp1 = repmat(squeeze(xrrt_paths(i,1,1:n1-1)),[1 n2-1]); 
        rp1 = repmat(squeeze(xrrt_paths(i,2,1:n1-1)),[1 n2-1]);
        rtp1 = repmat(squeeze(xrrt_paths(i,3,1:n1-1)),[1 n2-1]);

        xs1 = repmat(squeeze(xrrt_paths(i,1,2:n1)),[1 n2-1]) - xp1; 
        rs1 = repmat(squeeze(xrrt_paths(i,2,2:n1)),[1 n2-1]) - rp1; 
        rts1 = repmat(squeeze(xrrt_paths(i,3,2:n1)),[1 n2-1]) - rtp1; 

        % Outlet coordinates and differences in #2
        xp2 = repmat(outlet.x(1:end-1,1).',[n1-1 1]); 
        rp2 = repmat(outlet.r(1:end-1,1).',[n1-1 1]); 

        xs2 = repmat(outlet.x(2:end,1).',[n1-1 1]) - xp2; 
        rs2 = repmat(outlet.r(2:end,1).',[n1-1 1]) - rp2; 

        % Determine scalar constants
        t1 = ((xp2 - xp1).*rs2 - (rp2 - rp1).*xs2) ./ (xs1.*rs2 - rs1.*xs2);
        t2 = ((xp1 - xp2).*rs1 - (rp1 - rp2).*xs1) ./ (xs2.*rs1 - rs2.*xs1);

        % Remove all those that do not correspond to a mutual intersection
        t1(t1<0) = 0; t1(t1>1) = 0;
        t2(t2<0) = 0; t2(t2>1) = 0;

        % Calcuate x,y and build indexing matrix
        x = xp1 + t1.*xs1; r = rp2 + t2.*rs2; rt = rtp1 + t1.*rts1;

        t1(t1~=0) = 0.5;
        t2(t2~=0) = 0.5;

        t = t1+t2;
        t(t~=1) = 0;

        x = sum(sum(x.*t));
        r = sum(sum(r.*t));
        rt = sum(sum(rt.*t));

        if sum(sum(t)) == 1
            xrrt_line(i,:) = [x r rt];
        end
    end

    % Remove unfinished paths
    xrrt_line(isnan(xrrt_line(:,1)),:) = [];

    % Remove self intersecting regions
    [~,~,i,j] = intersections(xrrt_line(:,3),xrrt_line(:,2));
    i_intersect = [];
    for n = 1:size(i,1)
        if ceil(i(n)) < floor(j(n))
            i_intersect = [i_intersect ceil(i(n)):floor(j(n))];
        end
    end
    xrrt_line(i_intersect,:) = [];
else
    % Interpolate x and r coordinate for each line
    
    % Determine cut coordinate
    x_cut = mean(outlet.x(:));
    
    % Loop over all streamlines and interpolate x and r
    xrrt_line = zeros(size(xrrt_paths,2),3);
    for k = 1:size(xrrt_paths,2)
        q = find(isnan(xrrt_paths(:,k,1)) == 0);
        [~,i] = unique(xrrt_paths(q,k,1));
        xrrt_line(k,2) = interp1(xrrt_paths(q(i),k,1),xrrt_paths(q(i),k,2),x_cut,'pchip');
        xrrt_line(k,3) = interp1(xrrt_paths(q(i),k,1),xrrt_paths(q(i),k,3),x_cut,'pchip');
        xrrt_line(k,1) = x_cut;
    end
    
    % Remove unfinished paths
    xrrt_line(isnan(xrrt_line(:,2)),:) = [];
end

% Repeat +/- a pitch
nblade = double(outlet.bv.nblade);
xrrt_line = [xrrt_line(:,1:2) xrrt_line(:,3) - xrrt_line(:,2) * pi * 2 / nblade ; xrrt_line ; ...
    xrrt_line(:,1:2) xrrt_line(:,3) + xrrt_line(:,2) * pi * 2 / nblade];

end


function [mass_av_property,mass_flux,area] = ...
    ts_mass_average_nonuniform_2(c,property,xrrt_line,plot_stuff)
% Mass average an area divided by xrrt_line

if iscell(property) == 0
    temp = property; clear property; property{1} = temp;
end

% Work in double precision
c = ts_secondary(c);
property_names = fieldnames(c);
for i = 1:length(property_names)
    if isa(c.(property_names{i}), 'single') == 1
        c.(property_names{i}) = double(c.(property_names{i}));
    end
end
    
% Check if plane is structured or unstructured
varnames = {'ro' 'rovx' 'rovr' 'rorvt' 'roe'};
if isfield(c,'faces') == 0

    % Create triangulation of output plane through use of isosurface (marching
    % cubes)
    [F, V] = isosurface(repmat(c.rt,[1 1 2]),repmat(c.r,[1 1 2]),...
        repmat(c.x,[1 1 2]),cat(3,zeros(size(c.rt)), reshape(1:numel(c.rt),size(c.rt))),0.5);

    % Vertex and triangle coordinates
    V0 = V(F(:,1),:); V1 = V(F(:,2),:); V2 = V(F(:,3),:);
    rt_tri = [V0(:,1) V1(:,1) V2(:,1)]; r_tri = [V0(:,2) V1(:,2) V2(:,2)];
    x_tri = [V0(:,3) V1(:,3) V2(:,3)]; 

    % Flow properties
    for v = 1:length(varnames)
        [~,var_temp] = isosurface(repmat(c.rt,[1 1 2]),repmat(c.r,[1 1 2]),...
            repmat(c.(varnames{v}),[1 1 2]),cat(3,zeros(size(c.rt)),...
            reshape(1:numel(c.rt),size(c.rt))),0.5);
        h_tri.(varnames{v}) = [var_temp(F(:,1),3) var_temp(F(:,2),3) var_temp(F(:,3),3)];
    end
else
    % Get triangulated data from unstructured cut
    for v = 1:length(varnames)
        h_tri.(varnames{v}) = c.(varnames{v})(c.faces);
    end
    x_tri = c.x(c.faces); 
    r_tri = c.r(c.faces);
    rt_tri = c.rt(c.faces);
end

% Determine which triangles are in area with inpoly
% q = inpoly(mean(,xrrt_line(:,[3 2]));

% Areas of the triangles
t_tri = rt_tri./r_tri; t_av = sum(t_tri,2) / 3;

Ax = 0.5 * ((r_tri(:,1) .* (t_tri(:,1) - t_av) - r_tri(:,3) .* (t_tri(:,3) - t_av))...
    .* (r_tri(:,2) - r_tri(:,1)) -...
    (r_tri(:,1) .* (t_tri(:,1) - t_av) - r_tri(:,2) .* (t_tri(:,2) - t_av)) .*...
    (r_tri(:,3) - r_tri(:,1)));

Ar = 0.5 * ((r_tri(:,1) .* (t_tri(:,1) - t_av) - r_tri(:,3) .* (t_tri(:,3) - t_av))...
    .* (x_tri(:,2) - x_tri(:,1)) -...
    (r_tri(:,1) .* (t_tri(:,1) - t_av) - r_tri(:,2) .* (t_tri(:,2) - t_av)) .*...
    (x_tri(:,3) - x_tri(:,1)));

At = 0.5 * ((r_tri(:,1) - r_tri(:,3)) .* (x_tri(:,2) - x_tri(:,1)) -...
    (r_tri(:,1) - r_tri(:,2)) .* (x_tri(:,3) - x_tri(:,1)));

% Average flow properties
for v = 1:length(varnames)
    h_av.(varnames{v}) = sum(h_tri.(varnames{v}),2)/3;
end
h_av.rovt = sum(h_tri.rorvt ./ r_tri,2)/3;

% Get mass fluxes and areas
dM = Ax.*h_av.rovx + Ar.*h_av.rovr + At.*h_av.rovt;
mag_dA = (Ax.^2 + Ar.^2 + At.^2).^0.5;

% Get triangles that are in divided region using inpoly
h_av.x = sum(x_tri,2) / 3; h_av.r = sum(r_tri,2) / 3; h_av.rt = sum(rt_tri,2) / 3;
q = inpoly([h_av.rt h_av.r],xrrt_line(:,[3 2]));

% Perform mass averaging on property list
Bottom = sum(dM(q));
area = sum(mag_dA(q));

h_av.av = c.av; h_av.bv = c.bv;
h_av = ts_secondary(h_av);
mass_av_property = zeros(1,length(property));
for v = 1:length(property)
    Top = sum(h_av.(property{v})(q) .* dM(q));
    mass_av_property(1,v) = Top / Bottom;
end

mass_flux = Bottom;
    
end

function [mass_av_property, mass_flux, area] = ...
    ts_mass_average_nonuniform(c,property,xrrt_line, sector_type,plot_stuff)
% Calculate the mass average of a property in an arbitrarily bounded part of the mesh

% Choose cut type
% sector_type = 1 hub area bounded by a line above
% sector_type = 2 mid region bounded by two lines
% sector type = 3 casing area bounded by a line below

% Convert desired properties to cell array
if iscell(property) == 0
    temp = property; clear property; property{1} = temp;
end

% Work in double precision
c = ts_secondary(c);
property_names = fieldnames(c);
for i = 1:length(property_names)
    if isa(c.(property_names{i}), 'single') == 1
        c.(property_names{i}) = double(c.(property_names{i}));
    end
end

xrrt_line = double(xrrt_line);

% Extract intersecting line coordinates and trim length
r_line = xrrt_line(:,2); rt_line = xrrt_line(:,3);
drt_max = max(diff(rt_line)); 
rt_min = min(min(c.rt)) - drt_max; rt_max = max(max(c.rt)) + drt_max;
r_line(rt_line < rt_min | rt_line > rt_max) = [];
rt_line(rt_line < rt_min | rt_line > rt_max) = [];

% Check if plane is structured or unstructured
varnames = {'ro' 'rovx' 'rovr' 'rorvt' 'roe'};
if isfield(c,'faces') == 0

    % Create triangulation of output plane through use of isosurface (marching
    % cubes)
    [F, V] = isosurface(repmat(c.rt,[1 1 2]),repmat(c.r,[1 1 2]),...
        repmat(c.x,[1 1 2]),cat(3,zeros(size(c.rt)), reshape(1:numel(c.rt),size(c.rt))),0.5);

    % Vertex and triangle coordinates
    V0 = V(F(:,1),:); V1 = V(F(:,2),:); V2 = V(F(:,3),:);
    rt_tri = [V0(:,1) V1(:,1) V2(:,1) V0(:,1)]; r_tri = [V0(:,2) V1(:,2) V2(:,2) V0(:,2)];
    x_tri = [V0(:,3) V1(:,3) V2(:,3) V0(:,3)]; 

    % Flow properties
    for v = 1:length(varnames)
        [~,var_temp] = isosurface(repmat(c.rt,[1 1 2]),repmat(c.r,[1 1 2]),...
            repmat(c.(varnames{v}),[1 1 2]),cat(3,zeros(size(c.rt)),...
            reshape(1:numel(c.rt),size(c.rt))),0.5);
        h_tri.(varnames{v}) = [var_temp(F(:,1),3) var_temp(F(:,2),3) ...
            var_temp(F(:,3),3) var_temp(F(:,1),3)];
    end
else
    % Get triangulated data from unstructured cut
    for v = 1:length(varnames)
        h_tri.(varnames{v}) = c.(varnames{v})(c.faces(:,[1 2 3 1]));
    end
    x_tri = c.x(c.faces(:,[1 2 3 1])); 
    r_tri = c.r(c.faces(:,[1 2 3 1]));
    rt_tri = c.rt(c.faces(:,[1 2 3 1]));
end

if plot_stuff == 1
    figure()
    hold on
    z_max = max(c.(property{1})(:));
    contourf(c.rt,c.r,c.(property{1})-z_max-1,20)
    colormap bone;
%     trimesh(F,V(:,1),V(:,2),zeros(size(V,1),1))
    view(2)
    axis equal
    plot(rt_line,r_line,'r-')
end

% Preallocate intersection arrays
n1 = size(rt_tri,1); n2 = length(rt_line);
rt_i = zeros(n1,n2-1,3); r_i = zeros(n1,n2-1,3); x_i = zeros(n1,n2-1,3); 

for v = 1:length(varnames)
    h_i.(varnames{v}) = zeros(n1,n2-1,3);
end

% Loop over all triangle edges
for n = 1:3
    % Grid coordinates and differences in #1
    rp1 = repmat(squeeze(r_tri(:,n)),[1 n2-1]);
    rtp1 = repmat(squeeze(rt_tri(:,n)),[1 n2-1]);

    rs1 = repmat(squeeze(r_tri(:,n+1)),[1 n2-1]) - rp1; 
    rts1 = repmat(squeeze(rt_tri(:,n+1)),[1 n2-1]) - rtp1; 
    
    % Grid properties
    xp1 = repmat(squeeze(x_tri(:,n)),[1 n2-1]);
    xs1 = repmat(squeeze(x_tri(:,n+1)),[1 n2-1]) - xp1; 
    
    % Flow properties
    for v = 1:length(varnames)
        p1.(varnames{v}) = repmat(squeeze(h_tri.(varnames{v})(:,n)),[1 n2-1]);
        s1.(varnames{v}) = repmat(squeeze(h_tri.(varnames{v})(:,n+1)),[1 n2-1]) - p1.(varnames{v}); 
    end
    
    % Line coordinates and differences in #2
    rtp2 = repmat(rt_line(1:end-1,1).',[n1 1]); 
    rp2 = repmat(r_line(1:end-1,1).',[n1 1]); 

    rts2 = repmat(rt_line(2:end,1).',[n1 1]) - rtp2; 
    rs2 = repmat(r_line(2:end,1).',[n1 1]) - rp2; 

    % Determine scalar constants
    t1 = ((rtp2 - rtp1).*rs2 - (rp2 - rp1).*rts2) ./ (rts1.*rs2 - rs1.*rts2);
    t2 = ((rtp1 - rtp2).*rs1 - (rp1 - rp2).*rts1) ./ (rts2.*rs1 - rs2.*rts1);

    % Remove all those that do not correspond to a mutual intersection
    t1(t1<0) = 0; t1(t1>1) = 0;
    t2(t2<0) = 0; t2(t2>1) = 0;

    % Calcuate x,y and build indexing matrix
    rt = rtp1 + t1.*rts1; r = rp2 + t2.*rs2;
    
    % Calculate grid properties at intersections
    x = xp1 + t1.*xs1;
    for v = 1:length(varnames)
        h_temp.(varnames{v}) = p1.(varnames{v}) + t1.*s1.(varnames{v});
    end

    t1(t1~=0) = 0.5;
    t2(t2~=0) = 0.5;

    t = t1+t2;
    t(t~=1) = 0;

    r = r.*t;
    rt = rt.*t;

    rt(rt == 0) = nan; 
    r(r == 0) = nan; 
    
    if plot_stuff == 1
%         plot(rt,r,'b.')
    end

    % Record intersection values
    rt_i(:,:,n) = rt; r_i(:,:,n) = r; x_i(:,:,n) = x; 
    for v = 1:length(varnames)
        h_i.(varnames{v})(:,:,n) = h_temp.(varnames{v});
    end
    
end

% Split into unintersected triangles and intersected ones
q1 = (sum(sum(isnan(rt_i),3) == 3,2) == n2-1) == 0;
q2 = sum(sum(isnan(rt_i),3) == 3,2) == n2-1;

rt_tri_i = rt_tri(q1,:);
r_tri_i = r_tri(q1,:);
rt_tri = rt_tri(q2,1:3);
r_tri = r_tri(q2,1:3);

rt_i(q2,:,:) = [];
r_i(q2,:,:) = [];

% Grid properties 
x_tri_i = x_tri(q1,:); x_tri = x_tri(q2,1:3); x_i(q2,:,:) = [];

for v = 1:length(varnames)
    h_tri_i.(varnames{v}) = h_tri.(varnames{v})(q1,:);
    h_tri.(varnames{v}) = h_tri.(varnames{v})(q2,1:3);
    h_i.(varnames{v})(q2,:,:) = [];
end

% Loop over intersected triangles, divide into polygons and then into new
% triangles
rt_new = []; r_new = []; q_line = (1:length(r_line)-1).';
x_new = []; for v = 1:length(varnames); h_new.(varnames{v}) = []; end;
for n = 1:size(rt_i,1)
    % Pull out all points of interest
    q = isnan(rt_i(n,:,:)) == 0;
    rt = [rt_i(n,q) rt_tri_i(n,1:3)];
    r = [r_i(n,q) r_tri_i(n,1:3)];
    
    x = [x_i(n,q) x_tri_i(n,1:3)];
    for v = 1:length(varnames)
        h.(varnames{v}) = [h_i.(varnames{v})(n,q) h_tri_i.(varnames{v})(n,1:3)];
    end
    
    % Include line points within triangle if present
    q1 = q_line(sum(isnan(rt_i(n,:,:)) == 0,3) == 1);
    if max(diff(q1)) > 0
        q = zeros(length(r_line),length(q1)-1);
        for m = 1:length(q1)-1
            q(q1(m)+1:q1(m+1),m) = 1;
        end
        
        q(:,2:2:end) = [];
        q = sum(q,2) ~= 0;
        
        if sum(isnan(griddata(rt,r,x,rt_line(q),r_line(q)))) == 0
            x = [x griddata(rt,r,x,rt_line(q),r_line(q)).'];
            for v = 1:length(varnames)
                h.(varnames{v}) = [h.(varnames{v}) ...
                    griddata(rt,r,h.(varnames{v}),rt_line(q),r_line(q)).'];
            end

            rt = [rt rt_line(q).'];
            r = [r r_line(q).'];
        end
    end

    % Delaunay triangulation for new triangles
    dt = DelaunayTri(rt.',r.');
    
    if plot_stuff == 1
        triplot(dt,'Color','k')
    end
    
    rt_new = [rt_new ; dt.X(dt.Triangulation(:,1),1) dt.X(dt.Triangulation(:,2),1)...
        dt.X(dt.Triangulation(:,3),1)];
    r_new = [r_new ; dt.X(dt.Triangulation(:,1),2) dt.X(dt.Triangulation(:,2),2)...
        dt.X(dt.Triangulation(:,3),2)];
        
    x_new = [x_new ; x(dt.Triangulation(:,1)).' x(dt.Triangulation(:,2)).'...
        x(dt.Triangulation(:,3)).'];
    for v = 1:length(varnames)
        h_new.(varnames{v}) = [h_new.(varnames{v}) ; h.(varnames{v})(dt.Triangulation(:,1)).' ...
            h.(varnames{v})(dt.Triangulation(:,2)).' h.(varnames{v})(dt.Triangulation(:,3)).'];
    end
end

% Assemble all triangle data and exclude based on line location
rt_tri = [rt_tri ; rt_new]; r_tri = [r_tri ; r_new];
x_tri = [x_tri ; x_new]; 

for v = 1:length(varnames)
    h_tri.(varnames{v}) = [h_tri.(varnames{v}) ; h_new.(varnames{v})];
end

% Construct polygon of interesting area
rt_min = min(min(c.rt)) - drt_max; rt_max = max(max(c.rt)) + drt_max;
r_min = min(min(c.r)); r_max = max(max(c.r));

if sector_type == 1
    S_poly = [rt_min r_min ; rt_line r_line ; rt_max  r_min ; rt_min r_min];
elseif sector_type == 2
    S_poly = [rt_line r_line ; rt_line(1) r_line(1) ];
elseif sector_type == 3
    S_poly = [rt_min r_max ; rt_line r_line ; rt_max  r_max ; rt_min r_max];    
end

if plot_stuff == 1
    plot(S_poly(:,1),S_poly(:,2),'r-')
end

% Exclude triangles with inpoly
q1 = inpoly([rt_tri(:,1) r_tri(:,1)],S_poly);
q2 = inpoly([rt_tri(:,2) r_tri(:,2)],S_poly);
q3 = inpoly([rt_tri(:,3) r_tri(:,3)],S_poly);

q = q1 & q2 & q3;

rt_tri = rt_tri(q,:); r_tri = r_tri(q,:);
x_tri = x_tri(q,:); 
for v = 1:length(varnames)
    h_tri.(varnames{v}) = h_tri.(varnames{v})(q,:);
end

if plot_stuff == 1
    plot([rt_tri rt_tri(:,1)].',[r_tri r_tri(:,1)].','b-')
end

% Areas of the triangles
t_tri = rt_tri./r_tri; t_av = sum(t_tri,2) / 3;

Ax = 0.5 * ((r_tri(:,1) .* (t_tri(:,1) - t_av) - r_tri(:,3) .* (t_tri(:,3) - t_av))...
    .* (r_tri(:,2) - r_tri(:,1)) -...
    (r_tri(:,1) .* (t_tri(:,1) - t_av) - r_tri(:,2) .* (t_tri(:,2) - t_av)) .*...
    (r_tri(:,3) - r_tri(:,1)));

Ar = 0.5 * ((r_tri(:,1) .* (t_tri(:,1) - t_av) - r_tri(:,3) .* (t_tri(:,3) - t_av))...
    .* (x_tri(:,2) - x_tri(:,1)) -...
    (r_tri(:,1) .* (t_tri(:,1) - t_av) - r_tri(:,2) .* (t_tri(:,2) - t_av)) .*...
    (x_tri(:,3) - x_tri(:,1)));

At = 0.5 * ((r_tri(:,1) - r_tri(:,3)) .* (x_tri(:,2) - x_tri(:,1)) -...
    (r_tri(:,1) - r_tri(:,2)) .* (x_tri(:,3) - x_tri(:,1)));

% Average flow properties
for v = 1:length(varnames)
    h_av.(varnames{v}) = sum(h_tri.(varnames{v}),2)/3;
end
h_av.rovt = sum(h_tri.rorvt ./ r_tri,2)/3;

% Get mass fluxes and areas
dM = Ax.*h_av.rovx + Ar.*h_av.rovr + At.*h_av.rovt;
mag_dA = (Ax.^2 + Ar.^2 + At.^2).^0.5;

Bottom = sum(dM);
area = sum(mag_dA);

h_av.x = sum(x_tri,2) / 3; h_av.r = sum(r_tri,2) / 3; h_av.rt = sum(rt_tri,2) / 3;
h_av.av = c.av; h_av.bv = c.bv; h_av = ts_secondary(h_av);
mass_av_property = zeros(1,length(property));
for v = 1:length(property)
    Top = sum(h_av.(property{v}) .* dM);
    mass_av_property(1,v) = Top / Bottom;
end

mass_flux = Bottom;


end



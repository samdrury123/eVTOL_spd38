function [crr,r_corner,h_bl] = ts_shape_factor(c,b,xr_hub,xr_cas,xrrt_interp,plot_stuff)
% TS_SHAPE_FACTOR  Calculate boundary layer parameters and shape factors
%
%   [h_bar,h_corner] = TS_SHAPE_FACTOR(c,b,r_secs,plot_stuff)
%
%   c - structured cut of the blade o-mesh
%   b - assembled block of the o and h block on the suction side
%   xr_hub - matrix of hub line coordinates - [x1 r1 ; x2 r2 ; ...]
%   xr_cas - matrix of casing line coordinates - [x1 r1 ; x2 r2 ; ...]
%   xrrt_interp - matrix of desired coordinates to evaluate shape factors at
%   plot_stuff - 0 or 1 for showing working

% Default coordinates to evaluate at
if exist('xrrt_interp','var') == 0
    xrrt_interp = [];
end

% Default to show plots
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

% Default to calculate crr
if exist('calc_crr','var') == 0
    calc_crr = 1;
end

% Number of points across the boundary layer
nn = 30;


%% Prepare data

% Parameterise midspan section
r_hub = interp1(xr_hub(:,1),xr_hub(:,2),c.x(1,:),'pchip');
r_cas = interp1(xr_cas(:,1),xr_cas(:,2),c.x(1,:),'pchip');
r_nondim = (c.r(1,:) - r_hub) ./ (r_cas - r_hub);
[~,j_mid] = min(abs(r_nondim - 0.5));
t = bl_parameterise_section([c.x(:,j_mid) c.rt(:,j_mid)]);

% Calculate block secondary variables and convert to doubles
b = ts_secondary(b);
varnames = fieldnames(b);
for v = 1:length(varnames)
    if isnumeric(b.(varnames{v})) == 1
        b.(varnames{v}) = double(b.(varnames{v}));
    end
end

% Calculate non-dimensional radius
r_hub = interp1(xr_hub(:,1),xr_hub(:,2),b.x,'pchip');
r_cas = interp1(xr_cas(:,1),xr_cas(:,2),b.x,'pchip');
b.r_nondim = (b.r - r_hub) ./ (r_cas - r_hub);


%% Surface contraction ratio and corner separation extents
    
% Absolute distances based upon chord
d_min = t.chord / 50;
d_le = t.chord / 25;

% Open figure windows
if plot_stuff == 1
    h.cr = figure(); hold on; axis equal; box on;
    xlabel('Chord'); ylabel('Radius');
    
    h.bl = figure(); 
    subplot(2,1,1); hold on; grid on; box on;
    xlabel('Axial Distance'); ylabel('BL Thickness');
    subplot(2,1,2); hold on; grid on; box on;
    xlabel('Chord'); ylabel('Shape Factor');
    
    h.H_te = figure(); 
    subplot(1,2,1); hold on; grid on; box on; axis([0 5 0 1])
    xlabel('Surface Contraction Ratio'); ylabel('Radius');
    subplot(1,2,2); hold on; grid on; box on; axis([1.5 2.5 0 1]);
    xlabel('Trailing Edge Shape Factor'); ylabel('Radius');
end     

% Update attributes and write paraview file
g = cell(1); g{1} = c;
g{1}.attribute.bid = 0; g{1}.attribute.ni = size(g{1}.x,1);
g{1}.attribute.nj = size(g{1}.x,2); g{1}.attribute.nk = size(g{1}.x,3); 
ts_export_paraview(g,'paraview_temp.hdf5',{},[],1);

% Run paraview to track blade streamlines
[~,~] = system('vglrun pvpython ./PV/pv_track_contraction.py "./paraview_temp.xdmf"');
        
% Read csv file and calculate contractions
% xyz_raw = dlmread('./paraview_temp_mid_contraction0.csv',',',1,12);

% Cut into complete streamlines
% dx = diff(xyz_raw(:,1),1);
% q = dx < -d_min;
% i = find(q == 1); i = [0 ; i];
% xyz = nan(length(i)-1,max(diff(i))+2,3); 
% for j = 1:length(i)-1
%     xyz(j,1:(i(j+1)-i(j)),:) = reshape(xyz_raw(i(j)+1:i(j+1),:),1,[],3);
% end
% [t,r] = cart2pol(xyz(:,:,2),xyz(:,:,3));
% xrrt = cat(3,xyz(:,:,1),r,r.*t);

% Read streamline coordinates
s = pv_read_csv('./paraview_temp_mid_contraction0.csv');
xyz = permute(s.xyz,[2 1 3]); xrrt = bl_cart2pol(xyz);

% Get meridional lines near leading and trailing edges
[x,i] = min(c.x,[],1); xr_le = [x.'+ d_le c.r(i,:).'];
[x,i] = max(c.x,[],1); xr_te = [x.'- d_le c.r(i,:).'];

% Get intersection points
xr_1 = []; xr_2 = [];
for j = 1:size(xrrt,1)
    i_end = find(isnan(xrrt(j,:,1)),1)-1;
    if i_end > 1
        [x1_temp,r1_temp] = ...
            intersections(xr_le(:,1),xr_le(:,2),xrrt(j,1:i_end,1),xrrt(j,1:i_end,2),1);
        [x2_temp,r2_temp] = ...
            intersections(xr_te(:,1),xr_te(:,2),xrrt(j,1:i_end,1),xrrt(j,1:i_end,2),1);

        if length(x1_temp) == 1 && length(x2_temp) == 1
            xr_1 = [xr_1 ; x1_temp r1_temp];
            xr_2 = [xr_2 ; x2_temp r2_temp];
        end
    end
end

% Plot tracked streamlines and intersections with leading and trailing edges
if plot_stuff == 1
    figure(h.cr)
    plot(xrrt(:,:,1).',xrrt(:,:,2).','r-')
    plot(xr_le(:,1),xr_le(:,2),'k-');
    plot(xr_te(:,1),xr_te(:,2),'k-');
    plot(xr_1(:,1),xr_1(:,2),'b.')
    plot(xr_2(:,1),xr_2(:,2),'b.') 
end

% Calculate contraction rate
cr = diff(xr_1(:,2)) ./ diff(xr_2(:,2));
x_2 = 0.5 * (xr_2(1:end-1,1) + xr_2(2:end,1));
r_hub = interp1(xr_hub(:,1),xr_hub(:,2),x_2,'pchip');
r_cas = interp1(xr_cas(:,1),xr_cas(:,2),x_2,'pchip');
r = (0.5 * (xr_2(1:end-1,2) + xr_2(2:end,2)) - r_hub) ./ (r_cas - r_hub);
crr = [cr r];

% Plot contraction rate
if plot_stuff == 1
    figure(h.H_te)
    subplot(1,2,1);
    plot(cr,r,'-')
end

% Record spanwise extents of closed corner separations
xr_corner = xr_2([1 end],:);
r_corner = (xr_corner(:,2) - interp1(xr_hub(:,1),xr_hub(:,2),xr_corner(:,1),'pchip')) ./ ...
    (interp1(xr_cas(:,1),xr_cas(:,2),xr_corner(:,1),'pchip') - ...
    interp1(xr_hub(:,1),xr_hub(:,2),xr_corner(:,1),'pchip'));


%% Suction surface boundary layer parameters

% Block sizes
ni = size(b.x,1); nj = size(b.x,2); nk = size(b.x,3);

% Calculate distance from surface
b.d = cat(3,zeros(size(b.x,1),size(b.x,2),1),...
    cumsum((diff(b.x,1,3).^2 + diff(b.r,1,3).^2 + diff(b.rt,1,3).^2).^0.5,3));

% Generate velocity and density interpolants
f.ro = scatteredInterpolant(reshape(b.x,[],1),reshape(b.r,[],1),reshape(b.rt,[],1),reshape(b.ro,[],1));
f.V = scatteredInterpolant(reshape(b.x,[],1),reshape(b.r,[],1),reshape(b.rt,[],1),reshape(b.Vx,[],1));

% Choose relative or absolute velocity
if b.bv.rpm == 0; V = 'V'; else V = 'V_rel'; end;

% Estimate boundary layer thickness from max velocity
[~,k] = max(b.(V),[],3); d_Vmax = zeros(size(k));
for i = 1:ni
    for j = 1:nj
        
        % Interpolate finer velocity at a finer resolution
        k_temp = k(i,j)-2:k(i,j)+2; k_temp(k_temp < 1 | k_temp > nk) = [];
        d_temp = linspace(b.d(i,j,k_temp(1)),b.d(i,j,k_temp(end)),50);
        V_temp = interp1(squeeze(b.d(i,j,k_temp)),squeeze(b.(V)(i,j,k_temp)),d_temp,'spline');
        
        % Take maximum velocity value
        [~,k_temp] = max(V_temp); d_Vmax(i,j) = d_temp(k_temp);
    end
end

% Calculate blade normals
xrrt_blade = cat(3,b.x(:,:,1),b.r(:,:,1),b.rt(:,:,1));
a1 = xrrt_blade(2:end,2:end,:) - xrrt_blade(1:end-1,1:end-1,:);
a2 = xrrt_blade(1:end-1,2:end,:) - xrrt_blade(2:end,1:end-1,:);

n1 = cross(a1,a2); n1 = n1 ./ repmat(sum(n1.^2,3).^0.5,[1 1 3]);
n2 = zeros(size(xrrt_blade));
n2(2:end-1,2:end-1,:) = 0.25*(n1(1:end-1,1:end-1,:) + n1(2:end,1:end-1,:) +...
    n1(1:end-1,2:end,:) + n1(2:end,2:end,:));

% Decide whether to generate seed points or use specified values
if isempty(xrrt_interp) == 1

    % Seed points normal to blade scaled by estimated thickness
    d = repmat(reshape(1.5*hyperbolic_bunch(nn,100/nn,10/nn),[1 1 1 nn]),[ni nj 3]);
    xrrt_interp = repmat(xrrt_blade,[1 1 1 nn]) + d.* repmat(n2 .* repmat(d_Vmax,[1 1 3]),[1 1 1 nn]);

    % Limit calculation to between corner separations
    j_min = find(xrrt_blade(end,:,2) > xr_corner(1,2),1,'first');
    j_max = find(xrrt_blade(end,:,2) < xr_corner(2,2),1,'last');

    % Extract limited coordinates
    I = round(linspace(10,ni-5,10)); J = round(linspace(j_min,j_max,41));
    
else
    
    % Interpolate blade normals at specified points
    f.nx = scatteredInterpolant(reshape(xrrt_blade(:,:,1),[],1),reshape(xrrt_blade(:,:,2),[],1),...
        reshape(xrrt_blade(:,:,3),[],1),reshape(n2(:,:,1),[],1));
    f.nr = scatteredInterpolant(reshape(xrrt_blade(:,:,1),[],1),reshape(xrrt_blade(:,:,2),[],1),...
        reshape(xrrt_blade(:,:,3),[],1),reshape(n2(:,:,2),[],1));    
    f.nrt = scatteredInterpolant(reshape(xrrt_blade(:,:,1),[],1),reshape(xrrt_blade(:,:,2),[],1),...
        reshape(xrrt_blade(:,:,3),[],1),reshape(n2(:,:,3),[],1));   
    n2 = cat(3,f.nx(xrrt_interp(:,:,1),xrrt_interp(:,:,2),xrrt_interp(:,:,3)),...
        f.nr(xrrt_interp(:,:,1),xrrt_interp(:,:,2),xrrt_interp(:,:,3)),...
        f.nrt(xrrt_interp(:,:,1),xrrt_interp(:,:,2),xrrt_interp(:,:,3)));
    
    % Interpolate boundary layer thickness
    f.d_Vmax = scatteredInterpolant(reshape(xrrt_blade(:,:,1),[],1),reshape(xrrt_blade(:,:,2),[],1),...
        reshape(xrrt_blade(:,:,3),[],1),reshape(d_Vmax,[],1));
    d_Vmax = f.d_Vmax(xrrt_interp(:,:,1),xrrt_interp(:,:,2),xrrt_interp(:,:,3));
    
    % Seed points normal to blade scaled by estimated thickness
    d = repmat(reshape(1.5*hyperbolic_bunch(nn,100/nn,10/nn),[1 1 1 nn]),size(xrrt_interp));
    xrrt_interp = repmat(xrrt_interp,[1 1 1 nn]) + d.* repmat(n2 .* repmat(d_Vmax,[1 1 3]),[1 1 1 nn]);
    
    % Use all specifed points
    I = 1:size(xrrt_interp,1); J = 1:size(xrrt_interp,2);
end

% Shape desired coordinates into lists
x_interp = reshape(xrrt_interp(I,J,1,:),[],1);
r_interp = reshape(xrrt_interp(I,J,2,:),[],1);
rt_interp = reshape(xrrt_interp(I,J,3,:),[],1);

% Reshape back into point array
x_bl = reshape(x_interp,[length(I) length(J) nn]);
rt_bl = reshape(rt_interp,[length(I) length(J) nn]);
ro_bl = reshape(f.ro(x_interp,r_interp,rt_interp),[length(I) length(J) nn]);
V_bl = reshape(f.V(x_interp,r_interp,rt_interp),[length(I) length(J) nn]);
r_bl = reshape(r_interp,[length(I) length(J) nn]);

% Normal distance of interpolated points
d_bl = cat(3,zeros(length(I),length(J),1),...
    cumsum((diff(x_bl,1,3).^2 + diff(r_bl,1,3).^2 + diff(rt_bl,1,3).^2).^0.5,3));

% Re-interpolate to edge of boundary layer free stream velocity
[~,k] = max(V_bl,[],3);
for i = 1:size(d_bl,1)
    for j = 1:size(d_bl,2)
        d_temp = hyperbolic_bunch(nn,100/nn,10/nn)*d_bl(i,j,k(i,j));
        ro_bl(i,j,:) = interp1(squeeze(d_bl(i,j,:)),squeeze(ro_bl(i,j,:)),d_temp,'pchip');
        V_bl(i,j,:) = interp1(squeeze(d_bl(i,j,:)),squeeze(V_bl(i,j,:)),d_temp,'pchip');        
        d_bl(i,j,:) = d_temp;
    end
end

% Calculate displacement and momentum thickness
ro_99 = ro_bl(:,:,end); V_99 = V_bl(:,:,end); 
ro = 0.5 * (ro_bl(:,:,1:end-1) + ro_bl(:,:,2:end));
V = 0.5 * (V_bl(:,:,1:end-1) + V_bl(:,:,2:end));

dst = sum( (1 - (ro .* V) ./ repmat(ro_99 .* V_99,[1 1 nn-1])) .* diff(d_bl,1,3) , 3);
th = sum( ((ro .* V) ./ repmat(ro_99 .* V_99,[1 1 nn-1])) .*...
    (1 - V ./ repmat(V_99,[1 1 nn-1])) .* diff(d_bl,1,3) , 3);
de = sum( ((ro .* V) ./ repmat(ro_99 .* V_99,[1 1 nn-1])) .*...
    (1 - (V ./ repmat(V_99,[1 1 nn-1])).^2) .* diff(d_bl,1,3) , 3);

% Shape factor
H = dst ./ th;

% Trailing edge shape factor
[~,i] = max(H,[],1);
H_te = zeros(length(J),1); r_te = zeros(length(J),1);
for j = 1:length(J)
    H_te(j) = H(end-1,j);
    r_hub = interp1(xr_hub(:,1),xr_hub(:,2),x_bl(i(j),j,1),'pchip');
    r_cas = interp1(xr_cas(:,1),xr_cas(:,2),x_bl(i(j),j,1),'pchip');
    r_te(j) = (r_bl(i(j),j,1) - r_hub) / (r_cas - r_hub);
end

% Plot trailing edge shape factor
if plot_stuff == 1
    figure(h.H_te);
    subplot(1,2,2);
    plot(H_te,r_te,'-')
end

% Calculate streamwise coordinate
s = cat(1,zeros(1,nj,1),cumsum(sum(diff(xrrt_blade,1,1).^2,3).^0.5,1));

% Plot momentum thickness growth
j_mid = round(length(J)/2);
if plot_stuff
    figure(h.bl)
    subplot(2,1,1);
    plot(s(I,J(j_mid)),dst(:,j_mid),'r.-');
    plot(s(I,J(j_mid)),th(:,j_mid),'k.-');
    legend('Displacement Thickness','Momentum thickness');
    plot(s(I,J),dst,'r-');
    plot(s(I,J),th,'k-');
    subplot(2,1,2);
    plot(s(I,J),H,'r-');
end

% Record shape factor
r_hub = interp1(xr_hub(:,1),xr_hub(:,2),x_bl(:,:,1),'pchip');
r_cas = interp1(xr_cas(:,1),xr_cas(:,2),x_bl(:,:,1),'pchip');
h_bl.s = s(I,J);
h_bl.r = (r_bl(:,:,1) - r_hub) ./ (r_cas - r_hub);
h_bl.H = H;
h_bl.H_te = H_te;
h_bl.r_te = r_te;
h_bl.dst = dst;
h_bl.th = th;
h_bl.de = de;


end
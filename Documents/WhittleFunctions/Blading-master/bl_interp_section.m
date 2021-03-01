function xrrt = bl_interp_section(xrrt_master,plot_stuff,s_clust,keep_le,flip)
% BL_SECTION  Interpolate a poorly or coarsely defined blade to a high quality and resolution
%
%   xrrt = BL_INTERP_SECTION(xrrt_master,plot_stuff,s_clust)
%
%   xrrt_master - 3D or cell array of original section definition
%   plot_stuff - 0 or 1 for showing working
%   s_clust - specification for concentrating point clustering
%   keep_le - 0 or 1 to keep the same index for leading and trailing edge points
%   flip - 0 or 1 to flip blade nose down or not
%   xrrt - interpolated high-res polar section coordinates
%
%   s_clust is 3D array
%       1-direction is list of all points
%       2-direction is [s_cen s_wid s_hei]
%       s_cen is non-dimensional coordinate to focus on
%       s_wid is width of clustering
%       s_hei is intensity of clustering

if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

if exist('s_cluster','var') == 0
    s_clust = [];
end

if exist('keep_le','var') == 0
    keep_le = 0;
end

if exist('flip','var') == 0
    flip = 0;
end

% Specify blade sizes
ni_edge = 100; ni_surf = 200;

%% Prepare blade for easy manipulation

% Ensure nose points down by checking rt-coordinates
[~,i_le] = min(xrrt_master(:,1)); [~,i_te] = max(xrrt_master(:,1)); 
if xrrt_master(i_le,3) > xrrt_master(i_te,3) && flip == 1
    xrrt_master(:,3) = - xrrt_master(:,3);
end

% Ensure points are specified clockwise by crossing adjacent edges together
a = squeeze(diff(xrrt_master,1));
b = cross(a(2:end,:),a(1:end-1,:));
if sum(sign(b(:,2))) < 0
    xrrt_master = flipdim(xrrt_master,1);
end

% Ensure blade is horizontal
gamma = atand((xrrt_master(i_te,3) - xrrt_master(i_le,3)) / (xrrt_master(i_te,1) - xrrt_master(i_le,1)));
xrrt_master(:,[1 3]) = [xrrt_master(:,1) * cosd(-gamma) - xrrt_master(:,3) * sind(-gamma)...
    xrrt_master(:,1) * sind(-gamma) + xrrt_master(:,3) * cosd(-gamma)];

%% Increase resolution on the suction and pressure side

% Split into two surfaces
if keep_le == 0
    t = bl_guess_camber([xrrt_master(:,1) xrrt_master(:,3)]);
    xrrt_1 = xrrt_master(t.i_1,:); xrrt_2 = xrrt_master(t.i_2,:);
else
    xrrt_2 = xrrt_master(1:round(size(xrrt_master,1)/2),:);
    xrrt_1 = xrrt_master(end:-1:round(size(xrrt_master,1)/2),:);
    t.chord = sum((xrrt_1(1,:) - xrrt_1(end,:)).^2)^0.5;
end

% Check points spacing to find clipped leading and trailing edges
d_1 = sum(diff(xrrt_1,1,1).^2,2).^0.5;
d_2 = sum(diff(xrrt_2,1,1).^2,2).^0.5;
i_a = find(d_1 > d_1(1)*2,1,'first');
i_b = find(d_2 > d_2(1)*2,1,'first');
i_c = find(d_1 > min(d_1(end-5:end))*2,1,'last');
i_d = find(d_2 > min(d_2(end-5:end))*2,1,'last');

% Cubic interpolation in x for to reduce discontinuity in spacing
if isempty(i_a) == 0 && isempty(i_b) == 0 && isempty(i_c) == 0 && isempty(i_d) == 0
    x_temp = linspace(xrrt_1(i_a,1),xrrt_1(i_c+1,1),ni_surf*5).';
    xrrt_temp = [x_temp interp1(xrrt_1(i_a:i_c+1,1),xrrt_1(i_a:i_c+1,2),x_temp,'pchip') ...
        interp1(xrrt_1(i_a:i_c+1,1),xrrt_1(i_a:i_c+1,3),x_temp,'pchip')];
    xrrt_1 = [xrrt_1(1:i_a-1,:) ; xrrt_temp ; xrrt_1(i_c+2:end,:)];

    x_temp = linspace(xrrt_2(i_b,1),xrrt_2(i_d+1,1),ni_surf*5).';
    xrrt_temp = [x_temp interp1(xrrt_2(i_b:i_d+1,1),xrrt_2(i_b:i_d+1,2),x_temp,'pchip') ...
        interp1(xrrt_2(i_b:i_d+1,1),xrrt_2(i_b:i_d+1,3),x_temp,'pchip')];
    xrrt_2 = [xrrt_2(1:i_b-1,:) ; xrrt_temp ; xrrt_2(i_d+2:end,:)];    

    xrrt_master = [xrrt_1 ; flipud(xrrt_2(1:end-1,:))];
end

% Remove points that are too close together
xrrt_coarse = xrrt_master;
d_tol = (max(xrrt_coarse(:,1)) - min(xrrt_coarse(:,1))) / 30000;
d = sum(diff(xrrt_coarse,1,1).^2,2).^0.5;
i = find(d < d_tol);
for j = length(i):-1:1
    xrrt_coarse(i(j),:) = 0.5 * sum(xrrt_coarse(i(j):i(j)+1,:),1);
    xrrt_coarse(i(j)+1,:) = [];
end
xrrt_master = xrrt_coarse;

% Plot master geometry
if plot_stuff == 1
    figure(); hold on; axis equal; grid on; box on;
    axis([min(xrrt_master(:,1))-0.1*t.chord max(xrrt_master(:,1))+0.1*t.chord...
        min(xrrt_master(:,3))-0.1*t.chord max(xrrt_master(:,3))+0.1*t.chord]);
    plot(xrrt_master(:,1),xrrt_master(:,3),'k-');
end

%% Pull out four overlapping curves for LE, sides and TE

% Chordwise coordinates to cut
c_le = [0.01 0.05]; c_te = [0.94 0.99];
c_cen = [0 0.015 0.985 1];

if keep_le == 0
    t = bl_guess_camber([xrrt_master(:,1) xrrt_master(:,3)]);
else
    t.i_2 = (1:round(size(xrrt_master,1)/2))';
    t.i_1 = (size(xrrt_master,1):-1:round(size(xrrt_master,1)/2))';
    xy_1 = xrrt_master(t.i_1,[1 3]); xy_2 = xrrt_master(t.i_2,[1 3]); 
    s_1 = [0 ; cumsum(sum(diff(xy_1,1,1).^2,2).^0.5)]; s_1 = s_1 / max(s_1);
    s_2 = [0 ; cumsum(sum(diff(xy_2,1,1).^2,2).^0.5)]; s_2 = s_2 / max(s_2);

    s_fit = hyperbolic_bunch(12,0.03,0.01)';
    xy_1 = [interp1(s_1,xy_1(:,1),s_fit) interp1(s_1,xy_1(:,2),s_fit)];
    xy_2 = [interp1(s_2,xy_2(:,1),s_fit) interp1(s_2,xy_2(:,2),s_fit)];

    % Fit a polynomial camberline
    xy_cl = 0.5 * (xy_1 + xy_2); x = linspace(xy_1(1,1),xy_1(end,1),100)';
    p = polyfit(xy_cl(2:end-1,1),xy_cl(2:end-1,2),4);
    xy_cl = [x polyval(p,x)];
    s_cl = [0 ; cumsum(sum(diff(xy_cl,1,1).^2,2).^0.5)]; s_cl = s_cl / max(s_cl);
    
    t.s_1 = s_1; t.s_2 = s_2; t.s_cl = s_cl; t.xy_cl = xy_cl;
end

% Find indices closest to cut coordinates
[~,i_1a] = min(abs(t.s_1 - c_le(1)));
[~,i_1b] = min(abs(t.s_1 - c_le(2)));
[~,i_1c] = min(abs(t.s_1 - c_te(1)));
[~,i_1d] = min(abs(t.s_1 - c_te(2)));    

[~,i_2a] = min(abs(t.s_2 - c_le(1)));
[~,i_2b] = min(abs(t.s_2 - c_le(2)));
[~,i_2c] = min(abs(t.s_2 - c_te(1)));
[~,i_2d] = min(abs(t.s_2 - c_te(2)));

% Cut curves out
s_ss_master = t.s_2(i_2a:i_2d);
s_ps_master = t.s_1(i_1a:i_1d);
    
i_le = [t.i_2(i_2b:-1:1) ; t.i_1(2:i_1b)];
i_ss = t.i_2(i_2a:i_2d);
i_ps = t.i_1(i_1a:i_1d);
i_te = [t.i_2(i_2c:end) ; t.i_1(end-1:-1:i_1c)]; 

xrrt_le_master = xrrt_master(i_le,:);
xrrt_te_master = xrrt_master(i_te,:);
xrrt_ss_master = xrrt_master(i_ss,:);
xrrt_ps_master = xrrt_master(i_ps,:);

% Ensure unique values in master arrays
[~,i] = unique(xrrt_le_master,'rows'); xrrt_le_master = xrrt_le_master(sort(i),:);
[~,i] = unique(xrrt_te_master,'rows'); xrrt_te_master = xrrt_te_master(sort(i),:);
[~,i] = unique(xrrt_ss_master,'rows'); xrrt_ss_master = xrrt_ss_master(sort(i),:);
[~,i] = unique(xrrt_ps_master,'rows'); xrrt_ps_master = xrrt_ps_master(sort(i),:);

% Interpolate coordinates for nose, tail and centres
xrt_c_master = cell(4,1);
for c = 1:4
    xrt_c_master{c} = [interp1(t.s_cl,t.xy_cl(:,1),c_cen(c)) ...
        interp1(t.s_cl,t.xy_cl(:,2),c_cen(c))];
end

% Plot centre points
if plot_stuff == 1
    for c = 1:4
        plot(xrt_c_master{c}(1),xrt_c_master{c}(2),'k.');
    end
end

% Plot overlapping curves
if plot_stuff == 1
%     plot(xrrt_le_master(:,1),xrrt_le_master(:,3),'k.','MarkerSize',10)
%     plot(xrrt_te_master(:,1),xrrt_te_master(:,3),'ro','MarkerSize',10)
%     plot(xrrt_ss_master(:,1),xrrt_ss_master(:,3),'b+','MarkerSize',10)
%     plot(xrrt_ps_master(:,1),xrrt_ps_master(:,3),'gx','MarkerSize',10)
end

%% Interpolate leading edge geometry

% Find approximate leading edge angle
phi_nose = atan2(xrt_c_master{2}(2) - xrt_c_master{1}(2),...
        xrt_c_master{2}(1) - xrt_c_master{1}(1));

% Calculate angle coordinate around leading edge from centre
phi_le = atan2(xrt_c_master{2}(2) - xrrt_le_master(:,3), ...
    xrt_c_master{2}(1) - xrrt_le_master(:,1)) - phi_nose;

% Calculate surface angle based on gradient
psi_le = atan2(diff(xrrt_le_master(:,3),1,1),diff(xrrt_le_master(:,1),1,1));
psi_le = 0.5 * (psi_le(2:end) + psi_le(1:end-1)); 
psi_le = [psi_le(1) ; psi_le ; psi_le(end)];
psi_le = [0 ; cumsum(abs(diff(psi_le,1,1)))];

% Calculate distance through surface
s_le = [0 ; cumsum(sum(diff(xrrt_le_master(:,[1 3]),1,1).^2,2).^0.5)];

% Determine cutoff centre based angle for interpolation
phi_min = -1.5; phi_max = 1.5;
q = phi_le > phi_min & phi_le < phi_max;

% Non-dimensionalise surface angle and distance
psi_le = psi_le(q); psi_le = (psi_le - min(psi_le)) / (max(psi_le) - min(psi_le));
s_le = s_le(q); s_le = (s_le - min(s_le)) / (max(s_le) - min(s_le));

% Phi based on minimum surface angle change
[psi_le,i] = unique(psi_le); phi_le_temp = phi_le(q); 
phi_1 = interp1(psi_le,phi_le_temp(i),linspace(0,1,round(ni_edge/2))');

% Phi based on surface distance change
phi_2 = interp1(s_le,phi_le(q),linspace(0,1,round(ni_edge/2))');

% Combine both distributions
phi = unique([phi_1 ; phi_2]);

% Smooth phi distribution and increase to desired resolution
phi = interp1(linspace(0,1,length(phi))',smooth(smooth(phi)),linspace(0,1,ni_edge)','pchip');

% Fit leading edge
xrrt_LE_fit = zeros(ni_edge,3);
f = fit(phi_le, xrrt_le_master(:,1), 'smoothingspline');
xrrt_LE_fit(:,1) = f(phi);
f = fit(phi_le, xrrt_le_master(:,2), 'smoothingspline');
xrrt_LE_fit(:,2) = f(phi);
f = fit(phi_le, xrrt_le_master(:,3), 'smoothingspline');
xrrt_LE_fit(:,3) = f(phi);

% Plot fitted leading edge
if plot_stuff == 1
    plot(xrrt_LE_fit(:,1),xrrt_LE_fit(:,3),'r.-')
end

%% Interpolate trailing edge geometry

% Find trailing edge angle
phi_tail = atan2(xrt_c_master{4}(:,2) - xrt_c_master{3}(:,2),...
    xrt_c_master{4}(:,1) - xrt_c_master{3}(:,1));

% Angle coordinate around trailing edge from centre
phi_te = atan2(xrrt_te_master(:,3) - xrt_c_master{3}(2),...
    xrrt_te_master(:,1) - xrt_c_master{3}(1)) - phi_tail;

% Calculate surface angle based on gradient
psi_te = atan2(diff(xrrt_te_master(:,3),1,1),-diff(xrrt_te_master(:,1),1,1));
psi_te = 0.5 * (psi_te(2:end) + psi_te(1:end-1)); 
psi_te = [psi_te(1) ; psi_te ; psi_te(end)];
psi_te = [0 ; cumsum(abs(diff(psi_te,1,1)))];

% Calculate distance through surface
s_te = [0 ; cumsum(sum(diff(xrrt_te_master(:,[1 3]),1,1).^2,2).^0.5)];

% Determine cutoff centre based angle for interpolation
phi_min = -1.5; phi_max = 1.5;
q = phi_te > phi_min & phi_te < phi_max;

% Non-dimensionalise surface angle and distance
psi_te = psi_te(q); psi_te = (psi_te - min(psi_te)) / (max(psi_te) - min(psi_te));
s_te = s_te(q); s_te = (s_te - min(s_te)) / (max(s_te) - min(s_te));

% Phi based on minimum surface angle change
[psi_te,i] = unique(psi_te); phi_te_temp = phi_te(q); 
phi_1 = interp1(psi_te,phi_te_temp(i),linspace(0,1,round(ni_edge/2))');

% Phi based on surface distance change
phi_2 = interp1(s_te,phi_te(q),linspace(0,1,round(ni_edge/2))');

% Combine both distributions
phi = unique([phi_1 ; phi_2]);

% Smooth phi distribution and increase to desired resolution
phi = interp1(linspace(0,1,length(phi))',smooth(smooth(phi)),linspace(0,1,ni_edge)','pchip');

% Interpolate trailing edge
xrrt_TE_fit = zeros(ni_edge,3);

f = fit(phi_te(:), xrrt_te_master(:,1), 'smoothingspline');
xrrt_TE_fit(:,1) = f(phi);
f = fit(phi_te(:), xrrt_te_master(:,2), 'smoothingspline');
xrrt_TE_fit(:,2) = f(phi);
f = fit(phi_te(:), xrrt_te_master(:,3), 'smoothingspline');
xrrt_TE_fit(:,3) = f(phi);

% Plot fitted trailing edge
if plot_stuff == 1
    plot(xrrt_TE_fit(:,1),xrrt_TE_fit(:,3),'r.-')
end

%% Interpolate suction and pressure side geometry

ds_le_ss = abs(diff(interp1(xrrt_ss_master(:,1),s_ss_master,xrrt_LE_fit(1:2,1))));
ds_le_ps = abs(diff(interp1(xrrt_ps_master(:,1),s_ps_master,xrrt_LE_fit(end-1:end,1))));
ds_te_ss = abs(diff(interp1(xrrt_ss_master(:,1),s_ss_master,xrrt_TE_fit(end-1:end,1))));
ds_te_ps = abs(diff(interp1(xrrt_ps_master(:,1),s_ps_master,xrrt_TE_fit(1:2,1))));

% Interpolate suction and pressure surface
xrrt_SS_fit = zeros(ni_surf,3); xrrt_PS_fit = zeros(ni_surf,3);

% Smooth suction surface geometric progression
s_ss_le = interp1(xrrt_ss_master(:,1),s_ss_master,xrrt_LE_fit(1,1),'pchip');
s_ss_te = interp1(xrrt_ss_master(:,1),s_ss_master,xrrt_TE_fit(end,1),'pchip');
s_tot = s_ss_te - s_ss_le;
s = hyperbolic_bunch(ni_surf+2,ds_te_ss/s_tot,ds_le_ss/s_tot) * s_tot + s_ss_le;

% Cluster points closer to specified s locations
if isempty(s_clust) == 0
    for j = 1:size(s_clust,1)
        s_cen = s_clust(j,1); s_wid = s_clust(j,2); s_hei = s_clust(j,3);
        f = zeros(size(s)); 
        s_f = s(s > s_cen - s_wid & s < s_cen + s_wid);
        f(s > s_cen - s_wid & s < s_cen + s_wid) = -s_hei * sin(pi * (s_f - s_cen) / s_wid);
        s = s + f;
    end
end

xrrt_SS_fit(:,1) = interp1(s_ss_master,xrrt_ss_master(:,1),s(2:end-1),'pchip');
xrrt_SS_fit(:,2) = interp1(s_ss_master,xrrt_ss_master(:,2),s(2:end-1),'pchip');
xrrt_SS_fit(:,3) = interp1(s_ss_master,xrrt_ss_master(:,3),s(2:end-1),'pchip');

if plot_stuff == 1
    plot3(xrrt_SS_fit(:,1),xrrt_SS_fit(:,3),xrrt_SS_fit(:,2),'b.-')
end

% Smooth pressure surface geometric progression
s_ps_le = interp1(xrrt_ps_master(:,1),s_ps_master,xrrt_LE_fit(end,1),'pchip');
s_ps_te = interp1(xrrt_ps_master(:,1),s_ps_master,xrrt_TE_fit(1,1),'pchip');
s_tot = s_ps_te - s_ps_le;
s = hyperbolic_bunch(ni_surf+2,ds_te_ps/s_tot,ds_le_ps/s_tot) * s_tot + s_ps_le;

% Cluster points closer to specified s locations
if isempty(s_clust) == 0
    for j = 1:size(s_clust,1)
        s_cen = s_clust(j,1); s_wid = s_clust(j,2); s_hei = s_clust(j,3);
        f = zeros(size(s)); 
        s_f = s(s > s_cen - s_wid & s < s_cen + s_wid);
        f(s > s_cen - s_wid & s < s_cen + s_wid) = -s_hei * sin(pi * (s_f - s_cen) / s_wid);
        s = s + f;
    end
end

xrrt_PS_fit(:,1) = interp1(s_ps_master,xrrt_ps_master(:,1),s(2:end-1),'pchip');
xrrt_PS_fit(:,2) = interp1(s_ps_master,xrrt_ps_master(:,2),s(2:end-1),'pchip');
xrrt_PS_fit(:,3) = interp1(s_ps_master,xrrt_ps_master(:,3),s(2:end-1),'pchip');

if plot_stuff == 1
    plot3(xrrt_PS_fit(:,1),xrrt_PS_fit(:,3),xrrt_PS_fit(:,2),'b.-')
end

%% Return interpolated geometry

% Assemble all sections into array
xrrt = [xrrt_LE_fit(ni_edge/2:-1:1,:,:) ; xrrt_SS_fit ; flipud(xrrt_TE_fit) ;...
    flipud(xrrt_PS_fit) ; xrrt_LE_fit(end:-1:ni_edge/2,:,:)];

if plot_stuff == 1
%     plot3(xrrt(:,1),xrrt(:,3),xrrt(:,2),'g-')
end

% Rotate blade back
xrrt(:,[1 3]) = cat(3,xrrt(:,1) * cosd(gamma) - xrrt(:,3) * sind(gamma),...
    xrrt(:,1) * sind(gamma) + xrrt(:,3) * cosd(gamma));

end
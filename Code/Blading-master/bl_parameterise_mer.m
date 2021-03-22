function m = bl_parameterise_mer(xrrt,xr_hub,xr_cas,plot_stuff)
% BL_PARAMETERISE_MER  Parameterise meridional lines in the vicinity of a blade geometry
%
%   b = BL_PARAMETERISE_MER(xrrt,xr_hub,xr_cas,plot_stuff)
%   
%   xrrt - 3D float array of blade polar coordinates
%   xr_hub - 2D float list of hub coordinates
%   xr_cas - 2D float list of casing coordinates
%   plot_stuff - 0 or 1 for showing working
%   m - meridional definition output

if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Take unique hub and casing points
xr_hub = unique(xr_hub,'stable','rows'); xr_cas = unique(xr_cas,'stable','rows');

% Open figure window
if plot_stuff == 1
    figure(); hold on; grid on; box on; axis equal; view(2);
end

% Plot blade mesh
if plot_stuff == 1
    mesh(xrrt(:,:,1),xrrt(:,:,2),xrrt(:,:,3),zeros(size(xrrt)));
end

% Cut down hub and casing lines to half a chord up and downstream of blade
j_mid = round(size(xrrt,2)/2);
t = bl_guess_camber(xrrt(:,j_mid,[1 3])); dx = 0.5 * t.chord; i_te = t.i_te;
xr_hub = xr_hub(xr_hub(:,1) > xrrt(1,1,1) - dx & xr_hub(:,1) < xrrt(i_te,1,1) + dx,:);
xr_cas = xr_cas(xr_cas(:,1) > xrrt(1,end,1) - dx & xr_cas(:,1) < xrrt(i_te,end,1) + dx,:);

% Plot shortened meridional lines
if plot_stuff == 1
    plot(xr_hub(:,1),xr_hub(:,2),'k-')
    plot(xr_cas(:,1),xr_cas(:,2),'k-')
end

% Find platform and bleed gaps by discontinuities in angle in meridional lines
phi = atan2(diff(xr_hub(:,2)),diff(xr_hub(:,1))) * 360 / (2*pi);
i_dis_hub = find(abs(diff(phi)) > 2) + 1;

phi = atan2(diff(xr_cas(:,2)),diff(xr_cas(:,1))) * 360 / (2*pi);
i_dis_cas = find(abs(diff(phi)) > 2) + 1;

% Plot discontinuous points
if plot_stuff == 1
    plot(xr_hub(i_dis_hub,1),xr_hub(i_dis_hub,2),'r.')
    plot(xr_hub(i_dis_cas,1),xr_hub(i_dis_cas,2),'r.')
end    

% Find intersections of leading and trailing edges
xr_le = squeeze(xrrt(1,:,1:2)); xr_te = squeeze(xrrt(i_te,:,1:2)); 
[x,r] = intersections(xr_hub(:,1),xr_hub(:,2),xr_le(:,1),xr_le(:,2),1); xr_hub_le = [x(1) r(1)];
[x,r] = intersections(xr_cas(:,1),xr_cas(:,2),xr_le(:,1),xr_le(:,2),1); xr_cas_le = [x(1) r(1)];
[x,r] = intersections(xr_hub(:,1),xr_hub(:,2),xr_te(:,1),xr_te(:,2),1); xr_hub_te = [x(1) r(1)];
[x,r] = intersections(xr_cas(:,1),xr_cas(:,2),xr_te(:,1),xr_te(:,2),1); xr_cas_te = [x(1) r(1)];

% Find axial gaps between discontinuities
if isempty(i_dis_hub) == 0
    m.dx_gap_le_hub = abs(diff(xr_hub(i_dis_hub(xr_hub(i_dis_hub,1) < xr_hub_le(1)),1)));
    m.dx_gap_te_hub = abs(diff(xr_hub(i_dis_hub(xr_hub(i_dis_hub,1) > xr_hub_te(1),1))));
    m.x_gap_le_hub = xr_hub(find(i_dis_hub(xr_hub(i_dis_hub,1) < xr_hub_le(1)),1,'last'),1);
    m.x_gap_te_hub = xr_hub(find(i_dis_hub(xr_hub(i_dis_hub,1) > xr_hub_te(1)),1,'first'),1);
end
if isempty(i_dis_cas) == 0
    m.dx_gap_le_cas = abs(diff(xr_cas(i_dis_cas(xr_cas(i_dis_cas,1) < xr_cas_le(1),1))));
    m.dx_gap_te_cas = abs(diff(xr_cas(i_dis_cas(xr_cas(i_dis_cas,1) > xr_cas_te(1),1))));
    m.x_gap_le_cas = xr_cas(find(i_dis_cas(xr_cas(i_dis_cas,1) < xr_cas_le(1)),1,'last'),1);
    m.x_gap_te_cas = xr_cas(find(i_dis_cas(xr_cas(i_dis_cas,1) > xr_cas_te(1)),1,'first'),1);
end

% Find platform leading edge radii
% m.r_le_hub_plat = []; 
% if length(i_dis_hub) >= 2
%     if xr_hub(i_dis_hub(1),1) < xr_hub_le(1) && xr_hub(i_dis_hub(2),1) < xr_hub_le(1)
%         
%     end
% end

% Find leading and trailing edge areas
A_le = pi * (xr_cas_le(2)^2 - xr_hub_le(2)^2);
A_te = pi * (xr_cas_te(2)^2 - xr_hub_te(2)^2);

% Find leading and trailing edge mean radii
r_le = 0.5 * (xr_hub_le(2) + xr_cas_le(2));
r_te = 0.5 * (xr_hub_te(2) + xr_cas_te(2));

% Find gradients of lines at leading and trailing edges
drdx_hub = grad_mg(xr_hub(:,1),xr_hub(:,2));
drdx_cas = grad_mg(xr_cas(:,1),xr_cas(:,2));
drdx_hub_le = interp1(xr_hub(:,1),drdx_hub,xr_hub_le(1),'pchip');
drdx_hub_te = interp1(xr_hub(:,1),drdx_hub,xr_hub_te(1),'pchip');
drdx_cas_le = interp1(xr_cas(:,1),drdx_cas,xr_cas_le(1),'pchip');
drdx_cas_te = interp1(xr_cas(:,1),drdx_cas,xr_cas_te(1),'pchip');

% Return meridional parameters
m.A_le = A_le; m.A_te = A_te;
m.r_le = r_le; m.r_te = r_te;
m.drdx_hub_le = drdx_hub_le; m.drdx_hub_te = drdx_hub_te;
m.drdx_cas_le = drdx_cas_le; m.drdx_cas_te = drdx_cas_te;


end
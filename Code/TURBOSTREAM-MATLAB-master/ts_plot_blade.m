function [p,f,h,l] = ts_plot_blade(c,r_secs,h,col,plot_stuff,plot_cp)
% Plot interesting results on a TS blade row, peak suction, exit angle,
% loss and incidence

% Set defaults and open figure window if necessary
if exist('col','var') == 0 || isempty(col) == 1
    col = [0 0 0];
end
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end
if exist('h','var') == 0 || isempty(h) == 1 && plot_stuff == 1
    h = figure();
elseif plot_stuff == 1 && exist('h','var') ~= 0
    figure(h);
end

if exist('plot_cp','var') == 0
    plot_cp = [];
end

% Constants for use in shape factor calculation
s_bl_min = 0.45; s_bl_mid = 0.7;

%% Get interpolated blade section data

% Pre-allocate interpolated arrays
nj = length(r_secs); ni = c.blade.attribute.ni;
x_sects = zeros(ni,nj); rt_sects = zeros(ni,nj); P_sects = zeros(ni,nj); r_sects = zeros(ni,nj);
T_sects = cell(nj,1); s_sects = zeros(ni,nj); 

% Non-dimensional r-coordinate of blade
c.blade = ts_secondary(c.blade);
r_blade_hub = interp1(c.xr_hub(:,1),c.xr_hub(:,2),c.blade.x,'pchip','extrap');
r_blade_cas = interp1(c.xr_cas(:,1),c.xr_cas(:,2),c.blade.x,'pchip','extrap');
c.blade.r_nondim = (c.blade.r - r_blade_hub) ./ (r_blade_cas - r_blade_hub);

% Delete blade data above and below desired sections
varnames = fieldnames(c.blade); r_nondim = c.blade.r_nondim(1,:);
q = r_nondim > min(r_secs)-0.1*mean(diff(r_secs)) & r_nondim < max(r_secs)+0.1*mean(diff(r_secs));
for v = 1:length(varnames)
    if isnumeric(c.blade.(varnames{v})) == 1
        c.blade.(varnames{v}) = c.blade.(varnames{v})(:,q);
    end
end

% Re-interpolate sections at desired coordinates
for i = 1:ni
    [~,j] = unique(c.blade.r_nondim(i,:));
    x_sects(i,:) = interp1(c.blade.r_nondim(i,j),c.blade.x(i,j),r_secs,'pchip');
    rt_sects(i,:) = interp1(c.blade.r_nondim(i,j),c.blade.rt(i,j),r_secs,'pchip');
    r_sects(i,:) = interp1(c.blade.r_nondim(i,j),c.blade.r(i,j),r_secs,'pchip');
    P_sects(i,:) = interp1(c.blade.r_nondim(i,j),c.blade.P(i,j),r_secs,'pchip');
end

% Analyse blade surface
for j = 1:nj
    T_sects{j} = bl_parameterise_section([x_sects(:,j) rt_sects(:,j)]);
    s_sects(T_sects{j}.i_1,j) = T_sects{j}.s_cl_1_raw;
    s_sects(T_sects{j}.i_2,j) = T_sects{j}.s_cl_2_raw;
end
s_sects(1,:) = s_sects(end,:);

% Determine whether rotating
if c.blade.bv.rpm ~= 0
    use_rel = 1;
else
    use_rel = 0;
end


%% Incidence onto the leading edge

% Find stagnation point from maximum pressure and calculate incidence
psi_stag = zeros(1,nj); 
for j = 1:nj
        
    % Centre of leading edge at 3% chord
    xrt_le_cen = [interp1(T_sects{j}.s_cl,T_sects{j}.xy_cam(:,1),0.03,'pchip') ...
        interp1(T_sects{j}.s_cl,T_sects{j}.xy_cam(:,2),0.03,'pchip')];
    
    % Angle of surface from leading edge centre
    phi = atan2(xrt_le_cen(2) - rt_sects(:,j),xrt_le_cen(1) - x_sects(:,j)) * 360 / (2*pi);
    phi_le = atan2(xrt_le_cen(2) - T_sects{j}.xy_le(2),xrt_le_cen(1) - ...
        T_sects{j}.xy_le(1)) * 360 / (2*pi);
    phi = phi - phi_le;
    
    % Angle of blade surface
    drtdx = grad_mg(x_sects(:,j),rt_sects(:,j)); psi = atand(drtdx);
    
    % Increase resolution and interpolate to find peak pressure location
    [~, i] = max(P_sects(:,j));
    phi_temp = linspace(phi(i-2),phi(i+2),1000);
    P_temp = interp1(phi(i-2:i+2),P_sects(i-2:i+2,j),phi_temp,'spline');
    [~, i_max] = max(P_temp);    
    phi_stag = phi_temp(i_max);
    
    % Angle stagnation point makes with aerofoil relative to leading edge
    % angle
    if use_rel == 0
        psi_stag(j) = 90 - T_sects{j}.chi_le + ...
        interp1(phi(i - 2:i + 2),psi(i - 2:i + 2),phi_stag,'pchip');
    else
        psi_stag(j) = T_sects{j}.chi_le + 90 - ...
            interp1(phi(i - 2:i + 2),psi(i - 2:i + 2),phi_stag,'pchip');
    end
end


%% Blade pressure distributions

% Inlet pressures
r_nondim_in = ts_mass_average(c.inlet,'r_nondim',2).';

P_in = ts_mass_average(c.inlet,'P',2).';
P_in_secs = interp1(r_nondim_in,P_in,r_secs,'pchip');

if use_rel == 0
    Po_in = ts_mass_average(c.inlet,'Po',2).';
else
    Po_in = ts_mass_average(c.inlet,'Po_rel',2).';
end
Po_in_secs = interp1(r_nondim_in,Po_in,r_secs,'pchip');

% Pressure coefficient
Cp_sects = (P_sects - repmat(P_in_secs,[ni 1])) ./ repmat(Po_in_secs - P_in_secs,[ni 1]);

% Integrated pressure loading
Cp_tot = zeros(1,nj); for j = 1:nj; Cp_tot(j) = trapz(s_sects(:,j),-Cp_sects(:,j)); end;

% Find peak suction point and fit a polynomial to increase resolution
s_Cp_max = zeros(1,nj); i_ss = cell(nj,1); s_ss = cell(nj,1);
for j = 1:nj
    
    % Extract suction surface pressure and coordwise coordinate
    if mean(Cp_sects(T_sects{j}.i_2,j)) < mean(Cp_sects(T_sects{j}.i_1,j))
        i_ss{j} = T_sects{j}.i_2; s_ss{j} = T_sects{j}.s_cl_2_raw; n_ss = 2;
    else
        i_ss{j} = T_sects{j}.i_1; s_ss{j} = T_sects{j}.s_cl_1_raw; n_ss = 1;
    end
    
    % Bound extracted pressures with numbers close by for robustness
    Cp = Cp_sects(i_ss{j},j); s = s_ss{j};
    s = [s(1)-2e-9 ; s(1)-1e-9 ; s ; s(end)+1e-9 ; s(end)+2e-9];
    Cp = [Cp(1)+1e-9 ; Cp(1)+1e-9 ; Cp ; Cp(end)+1e-9 ; Cp(end)+1e-9];
    
    % Find peak suction
    [~, i] = min(Cp);
    
    % Increase resolution
    s_temp = linspace(s(i - 2),s(i + 2),1000);
    Q = polyval(polyfit(s(i - 2:i + 2),Cp(i - 2:i + 2),3),s_temp);
    [~, i_max] = min(Q);    
    s_Cp_max(j) = s_temp(i_max);
    
end


%% Exit whirl angles

% Get mass averaged exit angle
r_nondim_out = ts_mass_average(c.outlet,'r_nondim',2);
if use_rel == 0
    Alpha_out = ts_mass_average(c.outlet,'Alpha',2);
else
    Alpha_out = ts_mass_average(c.outlet,'Alpha_rel',2);
end
Alpha_out_secs = interp1(r_nondim_out,Alpha_out,r_secs,'pchip');


%% Radial loss coefficient 

% Extract inlet and outlet paramters
[s_in, mass_in] = ts_mass_average(c.inlet,'s',2);
if use_rel == 0
    V_in = ts_mass_average(c.inlet,'V',2);
else
    V_in = ts_mass_average(c.inlet,'V_rel',2);
end
[s_out, mass_out] = ts_mass_average(c.outlet,'s',2);
T_out = ts_mass_average(c.outlet,'T',2);

% Calculate mass flow fractions at inlet and exit
mass_in = [0 cumsum(mass_in)]; mass_out = [0 cumsum(mass_out)];
mass_in = mass_in / max(mass_in); mass_out = mass_out / max(mass_out);

% Calculate entropy loss coefficient
ent_coef = 2*T_out.*(s_out - interp1(mass_in,s_in,mass_out,'pchip')) ./...
    (interp1(mass_in,V_in,mass_out,'pchip').^2);
ent_coef_secs = interp1(r_nondim_out,ent_coef,r_secs,'pchip');
loss = 2 * ts_mass_average(c.outlet,'T',3) * (ts_mass_average(c.outlet,'s',3) -...
    ts_mass_average(c.inlet,'s',3)) ./ ts_mass_average(c.inlet,'V',3)^2;
        

%% Static pressure rise coefficient

% Inlet and outlet pressures
P_in = ts_mass_average(c.inlet,'P',2); P_out = ts_mass_average(c.outlet,'P',2);

% Midspan dynamic head
V_in = ts_mass_average(c.inlet,'V_rel',2); ro = ts_mass_average(c.inlet,'ro',2);
P_dyn = 0.5 * interp1(r_nondim_in,ro,0.5,'pchip') * interp1(r_nondim_in,V_in,0.5,'pchip')^2;

% Static pressure rise coefficient
C_static = (P_out - interp1(mass_in,P_in,mass_out,'pchip')) / P_dyn;
C_static_secs = interp1(r_nondim_out,C_static,r_secs,'pchip');


%% Exit mass flow function

% Axial Mach number
Mx = ts_mass_average(c.outlet,'M',2) .* cosd(ts_mass_average(c.outlet,'Alpha',2));

% Exit flow function
ga = c.outlet.av.ga;
mfunc = ga * (ga - 1)^(-0.5) * Mx .* (1 + 0.5 * (ga-1) * Mx.^2) .^ (-0.5 * (ga + 1) / (ga-1));
mfunc_secs = interp1(r_nondim_out,mfunc,r_secs,'pchip');


%% Shape factor profile and value at trailing edge

% Determine coordinates to evaluate shape factor at
ni_bl = 10; xrrt_surf = zeros(ni_bl,nj,3); s_bl = zeros(ni_bl,nj);
for j = 1:nj
    
    % Value of camberline coordinate at the trailing edge circle
    s_circ = s_ss{j}(i_ss{j} == T_sects{j}.i_circ(1) | i_ss{j} == T_sects{j}.i_circ(2));
    
    % Streamwise positions for boundary layer calculation
    s_bl(:,j) = linspace(s_bl_min,s_circ-0.01,ni_bl);
    
    % Interpolate surface coordinates 
    xrrt_surf(:,j,1) = interp1(s_ss{j},x_sects(i_ss{j},j),s_bl(:,j),'pchip');
    xrrt_surf(:,j,2) = interp1(s_ss{j},r_sects(i_ss{j},j),s_bl(:,j),'pchip');
    xrrt_surf(:,j,3) = interp1(s_ss{j},rt_sects(i_ss{j},j),s_bl(:,j),'pchip');
    
end

% Evaluate boundary layer shape factor
if isfield(c,'ss') == 1
    [~,~,bl] = ts_shape_factor(c.blade,c.ss,c.xr_hub,c.xr_cas,xrrt_surf,0);
    Hb = bl.H; Hb_te = Hb(end,:);
else
    Hb = zeros(size(s_bl)); Hb_te = zeros(1,size(s_bl,2));
end

% Integrate shape factor profile to determine whether the blade is front or rear loaded
aHb_1 = zeros(1,nj); aHb_2 = zeros(1,nj); 
for j = 1:nj
    
    % Interpolate shape factor at increased resolution
    s_temp = linspace(s_bl(1,j),s_bl(end,j),100);
    hb_temp = interp1(s_bl(:,j),Hb(:,j),s_temp,'pchip');
    
    % Interpolate differences from linear shape factor distribution
    dHb = hb_temp - interp1(s_bl([1 end],j),Hb([1 end],j),s_temp);

    % Split into two sections 
    q_1 = s_temp < s_bl_mid; q_2 = s_temp >= s_bl_mid;
    
    % Calculate integrated areas
    aHb_1(j) = trapz(s_temp(q_1),dHb(q_1)); aHb_2(j) = trapz(s_temp(q_2),dHb(q_2));
    
end


%% Plot results
if plot_stuff == 1
    
    % Switch to figure window
    figure(h); set(h,'Position',[721 924 1800 900])
    
    % Incidence
    subplot_tight(1,8,1,[0.08,0.01]); hold on; grid on; box on; 
    xlabel('Incidence'); set(gca,'yticklabel','');
    plot(psi_stag, r_secs, '.-', 'Color', col)
    axis auto; v = axis; axis([v(1:2) 0 1]);
    
    % Peak suction
    subplot_tight(1,8,2,[0.08,0.01]); hold on; grid on; box on; 
    xlabel('Peak Suction'); set(gca,'yticklabel','');
    plot(s_Cp_max, r_secs, '.-', 'Color', col)
    axis auto; v = axis; axis([v(1:2) 0 1]);
    
    % Shape factor areas
    subplot_tight(1,8,3,[0.08,0.01]); hold on; grid on; box on; 
    xlabel('H-bar Area'); set(gca,'yticklabel','');
    plot(aHb_1, r_secs, '.-', 'Color', col)
    plot(aHb_2, r_secs, '.-', 'Color', col)
    axis auto; v = axis; axis([v(1:2) 0 1]);
    
    % Trailing edge shape factor
    subplot_tight(1,8,4,[0.08,0.01]); hold on; grid on; box on; 
    xlabel('Exit H-bar'); set(gca,'yticklabel','');
    plot(Hb_te, r_secs, '.-', 'Color', col)
    axis auto; v = axis; axis([v(1:2) 0 1]);    
    
    % Exit whirl angle
    subplot_tight(1,8,5,[0.08,0.01]); hold on; grid on; box on; 
    xlabel('Exit Whirl'); set(gca,'yticklabel','');
    plot(Alpha_out, r_nondim_out, '-', 'Color', col)
    plot(Alpha_out_secs, r_secs, '.', 'Color', col)
    
    % Static pressure rise coefficient
    subplot_tight(1,8,6,[0.08,0.01]); hold on; grid on; box on; 
    xlabel('Pressure Rise'); set(gca,'yticklabel','');
    plot(C_static, r_nondim_out, '-', 'Color', col)
    plot(C_static_secs, r_secs, '.', 'Color', col)   
    
    % Exit flow function
    subplot_tight(1,8,7,[0.08,0.01]); hold on; grid on; box on; 
    xlabel('Exit Flow'); set(gca,'yticklabel','');
    plot(mfunc, r_nondim_out, '-', 'Color', col)
    plot(mfunc_secs, r_secs, '.', 'Color', col)   
    
    % Loss coefficient
    subplot_tight(1,8,8,[0.08,0.01]); hold on; grid on; box on; 
    xlabel('Entropy Coeff'); set(gca,'yticklabel','');
    l = plot(ent_coef, r_nondim_out, '-', 'Color', col);
    plot(ent_coef_secs, r_secs, '.', 'Color', col)
    axis([-0.02 0.3 0 1]);
end

% Section pressure coefficients
if isempty(plot_cp) == 0 
    figure(plot_cp); set(plot_cp,'position',[3 30 578 891]);
    for j = 1:nj
        a = subplot_tight(nj,1,j); hold on; grid on; box on; set(a,'xticklabel',[],'yticklabel',[]);
        plot(s_sects(:,j),-Cp_sects(:,j),'-','Color',col)
        axis auto; v = axis; axis([v(1:2) -1 v(4)]);
    end
end


%% Return design output arguments

% Design parameters evaluated on sections
p.psi_stag = psi_stag;
p.s_Cp_max = s_Cp_max; 
p.aHb_1 = aHb_1;
p.aHb_2 = aHb_2;
p.Alpha = Alpha_out_secs;
p.Hb_te = Hb_te;
p.loss = loss;
p.Cp_tot = Cp_tot;

% Radial variation of secondary design parameters
p.radial.ent_coef = ent_coef;
p.radial.Alpha = Alpha_out;
p.radial.C_static = C_static;
p.radial.mfunc = mfunc;
p.radial.r_nondim = r_nondim_out;


%% Find flow parameters for MISES

% Mach number
if use_rel == 0
    f.M = interp1(r_nondim_in,ts_mass_average(c.inlet,'M',2),r_secs);
else
    f.M = interp1(r_nondim_in,ts_mass_average(c.inlet,'M_rel',2),r_secs);
end

% Inlet angle
if use_rel == 0
    f.Alpha = interp1(r_nondim_in,ts_mass_average(c.inlet,'Alpha',2),r_secs);
else
    f.Alpha = interp1(r_nondim_in,ts_mass_average(c.inlet,'Alpha_rel',2),r_secs);
end

% AVDR
f.AVDR = interp1(r_nondim_in,ts_mass_average(c.inlet,'rovx',2)...
    ./ts_mass_average(c.outlet,'rovx',2),r_secs);

% Reynolds number
if use_rel == 0
    V = interp1(r_nondim_in,ts_mass_average(c.inlet,'V',2),r_secs);
else
    V = interp1(r_nondim_in,ts_mass_average(c.inlet,'V_rel',2),r_secs);
end
ro = interp1(r_nondim_in,ts_mass_average(c.inlet,'ro',2),r_secs);
f.Re = zeros(1,nj);
for j = 1:nj
    f.Re(j) = ro(j) * V(j) * T_sects{j}.tchord / c.inlet.av.viscosity;
end


end

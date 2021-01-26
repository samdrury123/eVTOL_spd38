function g = exp_traverse_calculate(e,c,x,r,N,des,coeffs,calculate_error,extrap_walls)

% Default to hodson coefficients
if exist('coeffs','var') == 0 || isempty(coeffs) == 1
    coeffs = 'hodson';
end

% Default to calculate probe error
if exist('calculate_error','var') == 0
    calculate_error = 1;
end

% Default to extrapolate traverse to walls
if exist('extrap_walls','var') == 0
    extrap_walls = 1;
end

% Load traverse data if filename is specified
if ischar(e) == 1
    load(e)
end

% Recall air properties
exp_air

% Compatible with old traverses
% if isfield(e,'P_raw') == 1
%     e.P = e.P_raw;
%     e.T = e.T_raw;
%     e.time = e.time_raw;
% end

% Simulate rig logging for old data
if isfield(e,'rig') == 0
    e.rig.P = reshape(e.P,[],size(e.P,3));
    e.rig.T = reshape(e.T,[],size(e.T,3));
    e.rig.time = reshape(e.time,[],size(e.time,3));
    
    % Sort in time order
    time_temp = datenum(e.rig.time);
    [~,i] = sort(time_temp);
    e.rig.P = e.rig.P(i,:);
    e.rig.T = e.rig.T(i,:);
    e.rig.time = e.rig.time(i,:);
    
    % Save other variables as defaults
    e.rig.Pa = ones(size(e.rig.P,1),1) * e.Pa;
    e.rig.Ha = ones(size(e.rig.P,1),1) * 30;
    e.rig.R = ones(size(e.rig.P,1),1) * air.R;
    e.rig.cp = ones(size(e.rig.P,1),1) * air.cp;
    e.rig.ga = ones(size(e.rig.P,1),1) * air.ga;
    e.rig.w = ones(size(e.rig.P,1),1) * des.U_mid/r.mid;
end

% Reduce old, multiple Re, calibration maps in size
if isfield(c,'P_raw') == 1 && size(c.P_raw,3) > 1
    c.P_raw = c.P_raw(:,:,1,:);
    c.C_Po = c.C_Po(:,:,1);
    c.C_P = c.C_P(:,:,1);
    c.C_Beta = c.C_Beta(:,:,1);
    c.C_Alpha = c.C_Alpha(:,:,1);
    c.Re = c.Re(:,:,1);
end

% Array sizes
np = size(e.rig.P,2); nt = size(e.rig.T,2); [~,nk] = size(e.r);

% Use stored calibration map if desired
if isempty(c) == 1
    c = {e.c};
end


%% Non-dimensionalise pressures and temperatures to average operating conditions at sea level

% Interpolate rig operating point for every traverse point
g = exp_rig_interp(e,N,r);

% Single out channels for reference
T_ref = mean(g.T_raw(:,:,N.T_ref),3); P_ref = g.P_raw(:,:,N.P_ref); 

% Calculate reference density
U_ref = g.U_mid; R_ref = g.R; ro_ref = P_ref ./ (R_ref .* T_ref);

% Normalise into coefficients with reference values
Cp = (g.P_raw - repmat(P_ref,[1 1 np])) ./ repmat(0.5 * ro_ref .* U_ref.^2,[1 1 np]);
Ct = g.T_raw ./ repmat(T_ref,[1 1 nt]);

% Determine whether to scale to mean operating point or outside reference
if isfield(des,'T') == 0
    T_des = mean(mean(T_ref)); P_des = mean(mean(P_ref)); U_des = mean(mean(U_ref));
else
    T_des = des.T; P_des = des.P; U_des = des.U_mid;
end

% Recalculate temperatures and pressures from scaled coefficients
ro_des = P_des / (air.R * T_des);
g.P_corr = Cp * 0.5 * ro_des * U_des^2 + P_des;
g.T_corr = Ct * T_des;
g.U_mid = U_des;

% Recalculate dynamic head
Po = g.P_corr(:,:,N.rig.Po); P = g.P_corr(:,:,N.rig.P); ga = air.ga; R = air.R;
To = mean(g.T_corr(:,:,N.rig.To),3);
M = ( ((Po ./ P).^((ga-1)./ga) - 1) .* (2./(ga - 1)) ).^0.5;
T = To ./ (1 + M.^2 .* (air.ga -1) / 2);
ro = P ./ (R .* T);
g.P_dyn = 0.5 * ro .* U_des.^2;

% Coordinates
% e.t = e.t - mean(mean(e.t));
g.x = ones(size(e.r)) * x;
g.r = e.r; g.rt = e.r .* (2*pi*e.t/360);
[g.y, g.z] = pol2cart(g.rt./g.r, g.r);

% figure(); hold on;
% plot(reshape(permute(g.P_raw,[2 1 3]),[],size(g.P_raw,3)),'-')
% plot(reshape(permute(g.P_corr,[2 1 3]),[],size(g.P_corr,3)),'--')
% figure(); hold on;
% plot(reshape(permute(g.T_raw,[2 1 3]),[],size(g.T_raw,3)),'-')
% plot(reshape(permute(g.T_corr,[2 1 3]),[],size(g.T_corr,3)5),'--')


%% Apply five hole calibration map to traverse pressures

% Get pressures and temperatures
T_mid = mean(g.T_corr(:,:,N.T_fhp),3);
P_FHP = g.P_corr(:,:,N.FHP);
% T_mid = mean(g.T_raw(:,:,N.T_fhp),3);
% P_FHP = g.P_raw(:,:,N.FHP);

% Use single calibration map if only one specified
if iscell(c) == 0
    c_temp{1} = c;
    c = c_temp;
end

% Correct rig pressure rise drift
if isfield(des,'psi') == 1
    dpsi = mean(g.psi) - des.psi;
    dP = dpsi * mean(g.P_dyn);
    P_FHP = P_FHP + dP;
end

% Correct pressure for finite hole location offsets
P_FHP = exp_fhp_shift(P_FHP,g.r,g.rt,e.probe);

% Apply calibration map
d = cell(size(c)); A = zeros(size(d));
for n = 1:length(c)
    c{n} = exp_fhp_calibrate(c{n}.Iota,c{n}.Tau,c{n}.P_raw,e.probe,c{n}.Ta,coeffs,0);
    try
        if calculate_error == 0; error; end;
        d{n} = exp_fhp_squares(c{n},P_FHP,e.probe,mean(e.Iota(:)),0,coeffs,0);
        A(n) = mean(mean(d{n}.A));
    catch
        d{n} = exp_fhp_calculate(c{n},P_FHP,e.probe,mean(e.Iota(:)),0,coeffs);
        A(n) = -1; d{n}.A = -1;
    end
end
[~,n] = min(A);
d = d{n};

% Assign flow parameters
g.P = d.P; g.Po = d.Po; 
g.Alpha = d.Alpha; g.Beta = d.Beta;


%% Calculate flow properties

% Fill in out of range parts by interpolating in tangential direction
varnames = {'Po','P','Alpha','Beta'};
for j = 1:size(g.Po,1)
    
    % Find out of range points
    q = isnan(g.Po(j,:));

    % Interpolate flow properties
    if sum(q) > 0
        for v = 1:length(varnames)
            g.(varnames{v})(j,q) = interp1(g.rt(j,q == 0),g.(varnames{v})(j,q == 0),...
                g.rt(j,q),'pchip');
        end
    end
end

% % Ensure flow direction is forwards
% g.P(g.Po < g.P) = g.Po(g.Po < g.P)-1;

% % Apply pressure correction
% g.Po = g.Po - dP;

% Calculate density and assign midspan temperature
g.ro = g.P ./ (air.R * T_mid);
g.T_mid = T_mid;

% Record FHP error
g.A = d.A; g.A_av = A;
    
% Velocity
M_temp = ((2/(air.ga-1)) * ((g.P ./ g.Po) .^ ((1-air.ga)/air.ga) - 1)).^0.5;
To_temp = (1 + (M_temp.^2)*(air.ga-1)/2) .* g.T_mid;
g.V = ((air.cp * To_temp).^0.5) .* ((air.ga-1).^0.5) .* M_temp .* ...
    ((1 + (M_temp.^2)*(air.ga-1)/2).^(-0.5));

% Velocity components
g.Vt = g.V .* sind(g.Alpha);
g.Vx = g.V .* cosd(g.Alpha) .* cosd(g.Beta);
g.Vr = g.V .* cosd(g.Alpha) .* sind(g.Beta);

% Rotational speed
g.bv.rpm = U_des * 60 / (2 * pi * r.mid);

% Primary flow variables
g.rovx = g.ro .* g.Vx;
g.rovr = g.ro .* g.Vr;
g.rorvt = g.ro .* g.r .* g.Vt;
g.roe = g.ro .* (air.cv * g.T_mid + 0.5 * g.V.^2);

% Check whether to extrapolate traverse
if extrap_walls > 0
    
    % Determine whether the hub is spinning or not
    if mean(g.Vt(1,:)) > mean(g.Vt(:)) && extrap_walls == 1
        Vt_hub = r.hub * U_des / r.mid;
    else
        Vt_hub = 0;
    end

    % Extrapolate measurement to the walls
    g.ro = g.ro([1 1:end end],:); g.rovx = [zeros(1,nk) ; g.rovx ; zeros(1,nk)];
    g.rovr = [zeros(1,nk) ; g.rovr ; zeros(1,nk)];
    g.rorvt = [g.ro(1,:) * r.hub * Vt_hub ; g.rorvt ; zeros(1,nk)];
    g.roe = g.roe([1 1:end end],:);

    % Define wall coorrdinates
    g.x = g.x([1 1:end end],:);
    g.rt = [r.hub * g.rt(1,:) ./ g.r(1,:) ; g.rt ; r.cas * g.rt(end,:) ./ g.r(end,:)];
    g.r = [ones(1,nk) * r.hub ; g.r ; ones(1,nk) * r.cas];

    % Copy reference values to the walls
    varnames = {'psi' 'phi' 'mrtp' 'nrt' 'Rp' 'P_dyn' 'w' 'P_raw' 'T_raw' 'R' 'cp' 'ga' ...
        'T_mid' 'P_corr' 'T_corr'};
    for v = 1:length(varnames)
        g.(varnames{v}) = g.(varnames{v})([1 1:end end],:,:);
    end
    
    % Delete casing extrapolation for half span traverses
    if g.r(end-1,1) == r.mid
        varnames = [varnames 'ro' 'rovx' 'rovr' 'rorvt' 'roe' 'x' 'r' 'rt'];
        for v = 1:length(varnames); g.(varnames{v})(end,:,:) = []; end;        
    end
    
    % Cubic interpolation of entire grid
    
end

% Remove all secondary flow variables
varnames = {'P','Po','Alpha','Beta','V','Vt','Vx','Vr'};
for v = 1:length(varnames)
    g = rmfield(g,varnames{v});
end

% Secondary variables and nondimensional r coordinate
g = ts_secondary(g);
g.r_nondim = (g.r - r.hub) / (r.cas - r.hub);


%% Pressure coefficients based on reference pressure and temperatures

% Pressure coefficients
g.Cpo = (g.Po - g.P_corr(:,:,N.P_ref)) ./ g.P_dyn;
g.Cp = (g.P - g.P_corr(:,:,N.P_ref)) ./ g.P_dyn;

% Probe coefficients
g.C_Alpha = d.C_Alpha; g.C_Beta = d.C_Beta;
g.C_Po = d.C_Po; g.C_P = d.C_P;

% Pressure and temperature plumbing
g.N = e.N;


%% Remove nearly overlapping points

% Identify overlapping point
k_near = find(abs(diff(g.t(1,:))) < 1e-5);

% Loop over all variable names and average near points
varnames = fieldnames(g);
for n = 1:length(k_near)
    for v = 1:length(varnames)
        if isnumeric(g.(varnames{v})) && size(g.(varnames{v}),2) == nk
            g.(varnames{v}) = [g.(varnames{v})(:,1:k_near(n)-1,:) ...
                mean(g.(varnames{v})(:,k_near(n):k_near(n)+1,:),2) ...
                g.(varnames{v})(:,k_near(n)+2:end,:)];
        end
    end
    nk = size(g.r,2);
end


end
function g = ts_secondary(g,use_double)
% Calculate and append all secondary flow variables for a matlab data cell array

if exist('use_double','var') == 0
    use_double = 0;
end

% Determine if g is a block (cut) or domain of blocks or array of domains
if iscell(g) == 1
    for n = 1:length(g);
        if iscell(g{n}) == 1
            for m = 1:length(g{n})
                if isfield(g{n}{1},'av') == 1
                    g{n}{m}.av = g{n}{m}.av;
                end
                g{n}{m} = calculate_block(g{n}{m},use_double);
            end
        else
            if isfield(g{1},'av') == 1
                g{n}.av = g{1}.av;
            end
            g{n} = calculate_block(g{n},use_double);
        end
    end
else
    g = calculate_block(g,use_double);
end
                  

end

function b = calculate_block(b,use_double)

% Get global variables and calculate gas constants and reference values
global Pref Tref
if isempty(Pref) == 1; Pref = 101325.8; Tref = 288.15; end;
    
% Get gas properties and block speed
if isfield(b,'cp') == 1
    cp = b.cp; ga = b.ga; mu = 1.8000e-05; use_air = 0;
elseif isfield(b, 'av') == 1
    cp = b.av.cp;
    ga = b.av.ga;
    if isfield(b.av,'viscosity') == 1; mu = b.av.viscosity; else mu = 1.8000e-05; end;
    use_air = 0;
else
    % Use air gas properties
    cp = 1003.308; ga = 1.40047; mu = 1.8000e-05;
    use_air = 1;
end

cv = cp ./ ga;
R = cp - cv;

if isfield(b, 'bv') == 1 && isfield(b.bv, 'rpm') == 1
    if length(b.bv) > 1
        b.bv = b.bv(1);
    end
    rpm = b.bv.rpm;
    omega = rpm * 2 * pi / 60;
else
    omega = 0;
end

% Singles
b.x = single(b.x); b.r = single(b.r); b.rt = single(b.rt);
b.ro = single(b.ro); b.roe = single(b.roe);
b.rovx = single(b.rovx); b.rorvt = single(b.rorvt); b.rovr = single(b.rovr);

% Velocities
b.Vx = b.rovx ./ b.ro;
b.Vr = b.rovr ./ b.ro;
b.Vt = (b.rorvt ./ (b.ro .* b.r));
b.V = (b.Vx.^2 + b.Vr.^2 + b.Vt.^2).^0.5;
b.Vm = sign(b.Vx) .* (b.Vx.^2 + b.Vr.^2).^0.5;
b.Vt_rel = (b.Vt - b.r*omega);
b.V_rel = (b.Vt_rel.^2 + b.Vm.^2).^0.5;

% Temperatures & Mach Numbers
E = b.roe ./ b.ro;
b.T = (E - 0.5 * b.V.^2) ./ cv;
b.M = b.V ./ (ga .* R .* b.T).^0.5;
b.M_rel = b.V_rel ./ (ga .* R .* b.T).^0.5;
b.To = b.T .* (1+ (ga-1)/2 .* (b.M).^2);
b.To_rel = b.T .* (1+ (ga-1)/2 .* (b.M_rel).^2);

% Pressures
b.P = b.ro .* R .* b.T;
b.Po = b.P .* (1+ (ga-1)/2 .* (b.M).^2).^(ga ./ (ga-1));
b.Po_rel = b.P .* (1+ (ga-1)/2 .* (b.M_rel).^2).^(ga ./ (ga-1));

% Angles
b.Alpha = atan2(b.Vt,b.Vm) * 360 / (2 * pi);
b.Alpha_rel = atan2(b.Vt_rel,b.Vm) * 360 / (2 * pi);
b.Beta = atan2(b.Vr,b.Vx) * 360 / (2 * pi);

% Entropy & Enthalpy
b.s = cp .* log(b.T/Tref) - R .* log(b.P/Pref);
b.sfunc = exp(-b.s ./ R);
b.Ts = b.s .* b.T;
b.h = cp .* (b.T-Tref);
b.ho = b.h + 0.5*b.V.^2;
b.ho_rel = b.h + 0.5*b.V_rel.^2;

% Cartesian coords
T = b.rt./b.r;
[Y, Z] = pol2cart(T, b.r);
b.y = Y; b.z = Z; b.t = T;

b.Vy = b.Vr .* cos(T) - b.Vt .* sin(T);
b.Vz = b.Vt .* cos(T) + b.Vr .* sin(T);

b.Vy_rel = b.Vr .* cos(T) - b.Vt_rel .* sin(T);
b.Vz_rel = b.Vt_rel .* cos(T) + b.Vr .* sin(T);

% Add non-dimensional r coords for cuts
if size(b.r,3) == 1
    b.r_nondim = (b.r - min(min(b.r))) / (max(max(b.r)) - min(min(b.r)));
end

% Add new block variables
roref = Pref ./ (R .* Tref);
sref = cp .* log(Tref) - R .* log(Pref);
b.bv.Pref = Pref; b.bv.Tref = Tref; 
b.bv.roref = roref; b.bv.sref = sref; 
b.bv.omega = omega;
b.bv.R = R; b.bv.cv = cv; b.bv.cp = cp; b.bv.ga = ga; b.bv.mu = mu;

% Add new application variables
b.av.cv = cv;
b.av.R = R;    
if use_air == 1
    b.av.cp = cp;
    b.av.ga = ga;
    b.av.viscosity = mu;
end

% Convert to doubles if required
if use_double == 1
    varnames = fieldnames(b);
    for v = 1:length(varnames)
        if isnumeric(b.(varnames{v})) == 1
            b.(varnames{v}) = double(b.(varnames{v}));
        end
    end
end

end
function g = exp_rig_interp(e,N,r)
% EXP_RIG_INTERP  Interpolate asynchronous rig logging onto the traverse time points

% Array sizes
np = size(e.rig.P,2); nt = size(e.rig.T,2); [nj,nk] = size(e.r);

% Construct time vectors for other data
time_rig = datenum(e.rig.time);
time_trav = datenum(reshape(e.time,[],6));

% Reinterpolate errors in atomospherics
Pa = e.rig.Pa; Pa_av = median(Pa);
q = abs((Pa - Pa_av) / Pa_av) > 0.01;
e.rig.Pa(q) = interp1(time_rig(q == 0),e.rig.Pa(q == 0),time_rig(q == 1),'linear','extrap');

% Reinterpolate errors in shaft speed
w = e.rig.w; w_av = median(w);
q = abs((w - w_av) / w_av) > 0.01;
e.rig.w(q) = interp1(time_rig(q == 0),e.rig.w(q == 0),time_rig(q == 1),'linear','extrap');

% Interpolate all rig properties at traverse time points
varnames = fieldnames(e.rig);
for v = 1:length(varnames)
    e_rig.(varnames{v}) = interp1(time_rig,e.rig.(varnames{v}),time_trav,'linear','extrap');
end

% Calculate missing flow properties
e_rig.N = N.rig;
e_rig = exp_rig_calculate(e_rig,r.mid);

% Calculate pressure rise coefficient and other 1D numbers
e_rig.P = e_rig.P - repmat(e_rig.P(:,N.rig.Pa) - e_rig.Pa,[1 np]);
g.P_dyn = 0.5 * e_rig.ro .* e_rig.U_mid.^2;
g.psi = (mean(e_rig.P(:,N.rig.P_out),2) - mean(e_rig.P(:,N.rig.Po),2)) ./ g.P_dyn;
g.phi = e_rig.phi; g.mrtp = e_rig.mrtp; g.nrt = e_rig.nrt; 
g.Rp = mean(e_rig.P(:,N.rig.P_out),2) ./ mean(e_rig.P(:,N.rig.Po),2);
g.w = e_rig.w;

% Correct pressures for chained DSAs
P_off = zeros(size(e.P));
if isfield(N,'P_off') == 1
    q = N.P_off ~= 0; P_off(:,:,q) = e.P(:,:,N.P_off(q));
end

% Get raw temperatures and pressures
g.P_raw = e.P - repmat(e.P(:,:,N.rig.Pa) - reshape(e_rig.Pa,[nj nk]),[1 1 np]) + P_off; 
g.T_raw = reshape(e_rig.T,[nj nk nt]);

% Blade speed and density
g.U_mid = r.mid * reshape(e_rig.w,[nj nk]);
g.ro = reshape(e_rig.ro,[nj nk]);

% Gas constants
g.R = e_rig.R; g.cp = e_rig.cp; g.ga = e_rig.ga;

% Reshape into arrays
varnames = fieldnames(g);
for v = 1:length(varnames)
    if size(g.(varnames{v}),2) == 1
        g.(varnames{v}) = reshape(g.(varnames{v}),size(e.r));
    end
end


end
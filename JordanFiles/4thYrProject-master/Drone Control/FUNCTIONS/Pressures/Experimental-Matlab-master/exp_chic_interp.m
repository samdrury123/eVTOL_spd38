function e = exp_chic_interp(e)
% EXP_CHIC_INTERP  Interpolate all data to the same rate

% Find array sizes
P_chan = size(e.P,2); T_chan = size(e.T,2);
if isfield(e,'V') == 1
    V_chan = size(e.V,2);
end

% Smooth torque measurments
if isfield(e.N,'V_torque') == 1
    e.V(:,2) = smooth(e.V(:,2),20000);
end

% Find time window covered by all data arrays
n_interp = size(e.P,1);
if isfield(e,'V') == 1
    time_min = (max([datenum(e.time_t(1,:)) datenum(e.time_p(1,:)) ...
        datenum(e.time_v(1,:))]));
    time_max = (min([datenum(e.time_t(end,:)) datenum(e.time_p(end,:)) ...
        datenum(e.time_v(end,:))]));
    time_interp = linspace(time_min,time_max,n_interp).';
else
    time_interp = datenum(e.time_p);
end

% Interpolate pressures onto desired time range
e.P_interp = zeros(n_interp,P_chan);
[~,i] = unique(datenum(e.time_p));
for n = 1:P_chan
    e.P_interp(:,n) = interp1(datenum(e.time_p(i,:)),e.P(i,n),time_interp,'linear');
end

% Interpolate temperatures
e.T_interp = zeros(n_interp,T_chan);
[~,i] = unique(datenum(e.time_t));
for n = 1:T_chan
    e.T_interp(:,n) = interp1(datenum(e.time_t(i,:)),e.T(i,n),time_interp,'linear');
end

% Interpolate pxie
if isfield(e,'V') == 1
    e.V_interp = zeros(n_interp,V_chan);
    [~,i] = unique(datenum(e.time_v));
    for n = 1:V_chan
        e.V_interp(:,n) = interp1(datenum(e.time_v(i,:)),e.V(i,n),time_interp,'linear');
    end
end

% Record time vector
e.time_interp = time_interp;

end
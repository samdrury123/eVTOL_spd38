function c = ts_shift_cut(c,t_shift,n_rep,straight_edges)
% TS_SHIFT_CUT  Move and repeat a cut plane for plotting
%
%   c = TS_SHIFT_CUT(c,t_shift,n_rep)
%
%   c - struct definining cut plane
%   t_shift - positive theta shift in radians before repeating and re-interpolating
%   n_rep - total number of passages to return
%   straight_edges - optional argument 0 or 1 to re-interpolate for constant theta edges

% Default to natural edges
if exist('straight_edges','var') == 0
    straight_edges = 0;
end

% Calculate secondary variables
c = ts_secondary(c);

% All fields in the plane
varnames = fieldnames(c);
for v = length(varnames):-1:1
    if isfloat(c.(varnames{v})) == 0 || ismatrix(c.(varnames{v})) == 0
        varnames(v) = [];
    end
end

% Size of the plane
[nj,nk] = size(c.x);

% Pitch of the cut plane
dt = abs(c.t(1,1) - c.t(1,end));

% Find mesh index at desired tangential shift
t_av = mean(c.t,1); t_rel = [0 cumsum(diff(t_av))];
[~,k_shift] = min(abs(t_rel - t_shift));

% Shift the coordinates
for v = 1:length(varnames)
    if strcmp(varnames{v},'t') == 0
        c.(varnames{v}) = [c.(varnames{v})(:,k_shift:end) c.(varnames{v})(:,2:k_shift)];
    else
        c.(varnames{v}) = [c.t(:,k_shift:end) c.t(:,2:k_shift) + dt];
    end
end

% Repeat the cut plane
for n = 1:n_rep-1
    for v = 1:length(varnames)
        if strcmp(varnames{v},'t') == 0
            c.(varnames{v}) = [c.(varnames{v}) c.(varnames{v})(:,2:end)];
        else
            c.(varnames{v}) = [c.t c.t(:,2:end) + dt * n];
        end
    end
end

% Centre the plane
t_av = mean(c.t,1); t_shift = 0.5 * (t_av(1) + t_av(end));
c.t = c.t - t_shift;

% Interpolate straight edges
if straight_edges == 1
    
    % Determine tangential coordinates to interpolate to
    t_min = t_av(1) - t_shift; t_max = t_av(end) - t_shift;
    t_interp = linspace(t_min,t_max,nk * 5);

    % Repeat tangential coordinate by +/- a pitch
    t = [c.t - (t_max - t_min) c.t(:,2:end-1) c.t + (t_max - t_min)];
    
    % Interpolate all variables
    for v = 1:length(varnames)
        var_temp = [c.(varnames{v}) c.(varnames{v})(:,2:end-1) c.(varnames{v})];
        var_interp = zeros(nj,nk*5);
        for j = 1:nj; var_interp(j,:) = interp1(t(j,:),var_temp(j,:),t_interp,'pchip'); end;
        c.(varnames{v}) = var_interp;
    end
    
    % Record tangential coordinates
    c.t = repmat(t_interp,[nj 1]);
end

% Recalculate primary and secondary variables
c.rt = c.t .* c.r; c = ts_secondary(c);


end
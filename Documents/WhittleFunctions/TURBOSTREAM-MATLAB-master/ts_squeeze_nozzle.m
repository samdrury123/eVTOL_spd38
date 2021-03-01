function g = ts_squeeze_nozzle(g,x_lim,xr_hub_old,xr_cas_old,xr_hub_new,xr_cas_new)
% TS_SQUEEZE_NOZZLE  Morph a mesh to a new nozzle shape

% Maximum and minimum axial coordinate
x_max = -inf; x_min = inf;
for bid = ts_get_bids(g)
    x_max = max(max(g{bid+1}.x(:)),x_max); 
    x_min = min(min(g{bid+1}.x(:)),x_min); 
end

% Interpolate radial shift
x = linspace(x_lim,x_max,100)';
dr_hub = interp1(xr_hub_new(:,1),xr_hub_new(:,2),x,'linear','extrap') - ...
    interp1(xr_hub_old(:,1),xr_hub_old(:,2),x,'linear','extrap');
dr_cas = interp1(xr_cas_new(:,1),xr_cas_new(:,2),x,'linear','extrap') - ...
    interp1(xr_cas_old(:,1),xr_cas_old(:,2),x,'linear','extrap');

% Extend shifts upstream
x = [x_min ; x]; dr_hub = [0 ; dr_hub]; dr_cas = [0 ; dr_cas];

% Loop over all blocks and check downstream of limit
for bid = ts_get_bids(g)
    if max(g{bid+1}.x(:)) > x_lim

        % Interpolate non-dimensional radius
        r_hub = interp1(xr_hub_old(:,1),xr_hub_old(:,2),g{bid+1}.x,'linear','extrap');
        r_cas = interp1(xr_cas_old(:,1),xr_cas_old(:,2),g{bid+1}.x,'linear','extrap');
        r_nondim = (g{bid+1}.r - r_hub) ./ (r_cas - r_hub);

        % Calculate radial shift
        dr_hub_temp = interp1(x,dr_hub,g{bid+1}.x,'linear','extrap'); 
        dr_cas_temp = interp1(x,dr_cas,g{bid+1}.x,'linear','extrap'); 
        dr = r_nondim .* (dr_cas_temp - dr_hub_temp) + dr_hub_temp;

        % Apply radial shift
        t = g{bid+1}.rt ./ g{bid+1}.r;
        g{bid+1}.r = g{bid+1}.r + dr;    
        g{bid+1}.rt = g{bid+1}.r .* t;

    end
end

end


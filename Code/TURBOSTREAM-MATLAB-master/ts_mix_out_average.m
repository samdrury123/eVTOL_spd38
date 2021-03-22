function mix = ts_mix_out_average(c,plot_stuff)
% Use a mix out average to average a flow solution tangentially

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Get secondary variables
c = ts_secondary(c);

% Cell relative tangential velocity
omega = c.bv.rpm * 2 * pi / 60;
c.rovt = c.rorvt ./ c.r - (c.r .* c.ro * omega);

% Determine whether cut is structured or unstructured
if isfield(c,'faces') == 1
    
    % Resolution of bins
%     nj = 10; dj = 0.1;
    nj = 21; dj = 0.02;
    
    % Get cartesian coordinates and face indices
    y = c.y; z = c.z; 
    faces = c.faces; f1 = faces(:,1); f2 = faces(:,2); f3 = faces(:,3);
    
    % Calculate the area of the face triangles using Heron's formula
    a1 = ((y(f2) - y(f1)).^2 + (z(f2) - z(f1)).^2).^0.5;
    a2 = ((y(f3) - y(f2)).^2 + (z(f3) - z(f2)).^2).^0.5;
    a3 = ((y(f1) - y(f3)).^2 + (z(f1) - z(f3)).^2).^0.5;
    d = sort([a1 a2 a3],2,'descend');
    a1 = d(:,1); a2 = d(:,2); a3 = d(:,3);
    Ax_list = 0.25*sqrt(4*a1.*a1.*a2.*a2 - (a1.*a1 + a2.*a2 - a3.*a3).^2);

    % Correct imaginary areas
    if any(imag(Ax_list)~=0)
        Ax_list = real(Ax_list);
    end
    
    % Average all properties
    varnames = fieldnames(c);
    for v = 1:length(varnames)
        if isfloat(c.(varnames{v})) == 1
            avg_list.(varnames{v}) = (1/3) * (c.(varnames{v})(f1) + c.(varnames{v})(f2) + ...
                c.(varnames{v})(f3));
        end
    end
    
    % Define radii to bin triangles in
    r_bin = (max(c.r) - min(c.r)) * hyperbolic_bunch(nj,dj,dj) + min(c.r);
    
    % Loop over all radii and collect triangles in range
    for j = 1:length(r_bin)-1
        
        % Look for average radii in range
        q = avg_list.r > r_bin(j) & avg_list.r < r_bin(j+1);
        nk = sum(q);
        
        % Loop over all variables and collect in matrix padded with zeros
        for v = 1:length(varnames)
            if isfloat(c.(varnames{v})) == 1
                avg.(varnames{v})(j,1:nk) = avg_list.(varnames{v})(q);
            end
        end
        
        % Collect area into matrix
        Ax(j,1:nk) = Ax_list(q);
    end
    
else

    % Calculate cell centred averages of all properties
    varnames = fieldnames(c);
    for v = 1:length(varnames)
        if isfloat(c.(varnames{v})) == 1
            avg.(varnames{v}) = 0.25 * (c.(varnames{v})(1:end-1,1:end-1) + ...
                c.(varnames{v})(2:end,1:end-1) + ...
                c.(varnames{v})(1:end-1,2:end) + ...
                c.(varnames{v})(2:end,2:end));
        end
    end

    % Calculate area components by cross product of diagonal vectors
    c.t = c.rt ./ c.r;
    t_av = 0.25 * (c.t(2:end,2:end) + c.t(1:end-1,1:end-1) + c.t(1:end-1,2:end) + c.t(2:end,1:end-1));
%     x1 = c.x(2:end,2:end) - c.x(1:end-1,1:end-1);
%     x2 = c.x(1:end-1,2:end) - c.x(2:end,1:end-1);
    r1 = c.r(2:end,2:end) - c.r(1:end-1,1:end-1);
    r2 = c.r(1:end-1,2:end) - c.r(2:end,1:end-1);
    rt1 = c.r(2:end,2:end) .* (c.t(2:end,2:end) - t_av)...
        - c.r(1:end-1,1:end-1) .* (c.t(1:end-1,1:end-1) - t_av);
    rt2 = c.r(1:end-1,2:end) .* (c.t(1:end-1,2:end) - t_av)...
        - c.r(2:end,1:end-1) .* (c.t(2:end,1:end-1) - t_av);

    Ax = 0.5 * (r1.*rt2 - r2.*rt1);
    % Ar = 0.5 * (x2.*rt1 - x1.*rt2);
    % At = 0.5 * (x1.*r2 - x2.*r1);
end

% Assume constant area of mixing and constant x coordinate

% Guess a density distribution
nk = sum(avg.ro ~= 0,2);
mix.ro = sum(avg.ro,2) ./ nk;

% Set up Mach number functions
ga = c.av.ga; cp = c.av.cp; R = c.av.R; cv = c.av.cv; 
M_interp = 0:0.001:2;
V_func_interp = (ga-1)^0.5 * M_interp .* (1 + 0.5 * (ga-1) * M_interp.^2) .^ (-0.5);
T_func_interp = (1 + 0.5 * (ga-1) * M_interp.^2) .^ (-1);

% Plot the density
if plot_stuff == 1
    figure(); hold on; grid on; box on; xlabel('Density'); ylabel('Radius');
    r = avg.r; r(r == 0) = nan; r = nanmean(r,2); plot(r,mix.ro,'-');
end

% Iterate on the density
for n = 1:5

    % Conservation of mass to get mixed out axial velocity
    mix.Vx = sum(avg.rovx .* Ax,2) ./ (sum(Ax,2) .* mix.ro);

    % Conservation of axial momentum to get mixed out static pressure
    mix.P = ( sum(Ax .* (avg.ro .* avg.Vx.^2 + avg.P),2) - mix.ro .* sum(Ax,2) .* mix.Vx.^2 ) ./ ...
        sum(Ax,2);

    % Conservation of tangential momentum to get mixed out tangential velocity
    mix.Vt = sum(Ax .* (avg.ro .* avg.r .* avg.Vx .* avg.Vt),2) ./ ...
        (sum(Ax .* avg.r,2) .* mix.ro .* mix.Vx);

    % Destruction of radial momentum
    mix.Vr = zeros(size(mix.Vx));

    % Total temperature from first law of thermodynamics
    mix.To = sum(Ax .* avg.rovx .* avg.To,2) ./ (sum(Ax,2) .* mix.ro .* mix.Vx);
    
    % Mach number from total temperature and velocity relation
    mix.V = (mix.Vx.^2 + mix.Vr.^2 + mix.Vt.^2).^0.5;
    V_func = mix.V .* (cp * mix.To) .^ (-0.5);
    mix.M = interp1(V_func_interp,M_interp,real(V_func));
    
    % Static temperature from Mach number relation
    T_func = interp1(M_interp,T_func_interp,mix.M);
    mix.T = T_func .* mix.To;
    
    % Density from ideal gas law
    mix.ro = mix.P ./ (R * mix.T);
    
    % Record flow condition in primary variables
    mix.rovx = mix.ro .* mix.Vx;
    mix.rovr = mix.ro .* mix.Vr;
    mix.rorvt = mix.ro .* mix.Vt .* sum(avg.r,2) ./ nk;
    mix.roe = mix.ro .* (cv * mix.T + 0.5 * mix.V.^2);    

    % Plot the density update
    if plot_stuff == 1
        r = avg.r; r(r == 0) = nan; r = nanmean(r,2); plot(r,mix.ro,'-'); 
    end
    
end

% Record mixed out coordinates in bins for unstructured cuts
if isfield(c,'faces') == 1
    avg.x(avg.x == 0) = nan; avg.rt(avg.rt == 0) = nan;
    mix.x = nanmean(avg.x,2); mix.rt = nanmean(avg.rt,2);
end

% Reassign properties back to nodes
varnames = fieldnames(mix);
for v = 1:length(varnames)
    if numel(mix.(varnames{v})) > 1
        var_new = zeros(length(mix.(varnames{v}))+1,1);
        var_new(2:end-1) = 0.5 * (mix.(varnames{v})(1:end-1) + ...
            mix.(varnames{v})(2:end));

        var_new(1) = mix.(varnames{v})(1);
        var_new(end) = mix.(varnames{v})(end);
        mix.(varnames{v}) = var_new;
    end
end

% Record mixed out coordinates
if isfield(c,'faces') == 1
    mix.r = r_bin';
else
    mix.r = mean(c.r,2);
    mix.x = mean(c.x,2);
    mix.rt = mean(c.rt,2);
end

% Gas properties
mix.av = c.av; mix.bv = c.bv;

% Get secondary mixed out variables
mix = ts_secondary(mix);

% Plot the final densit
if plot_stuff == 1
    plot(mix.r,mix.ro,'k-'); 
end

end
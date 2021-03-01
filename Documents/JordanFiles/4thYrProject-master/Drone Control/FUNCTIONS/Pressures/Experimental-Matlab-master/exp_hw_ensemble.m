function g = exp_hw_ensemble(e,frac,dt,plot_stuff,fit_order)
% EXP_HW_ENSEMBLE  Phase lock average hot-wire data with or without time interpolation

% Default to no shift
if exist('dt','var') == 0
    dt = 0;
end

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Default to quartic fit through all available points
if exist('fit_order','var') == 0
    fit_order = 4;
end

% Find falling edges of once per rev signals
% a_log = e.V_shaft(:,:,:,2:end) == 0 & e.V_shaft(:,:,:,1:end-1) == 1;

% Find rising and falling edges of shaft signal
a_fall = e.V_shaft(:,:,:,2:end) == 0 & e.V_shaft(:,:,:,1:end-1) == 1;
a_rise = e.V_shaft(:,:,:,2:end) == 1 & e.V_shaft(:,:,:,1:end-1) == 0;

% Loop over all hot-wire positions and ensemble average
rate = 1 / mean(diff(e.dtime)); f = mean(e.rig.w) / (2*pi);
[nj,nk,ni,~] = size(e.V_shaft); nt = round(frac * rate / f); nn = length(e.dtime);
V_av = zeros(nj,nk,ni,nt); V_p = zeros(nj,nk,ni,nt); w = zeros(nj,nk,ni); 
q = (1:nt) + round(nt * dt / (2 * pi * frac));
for i = 1:ni
    for j = 1:nj
        for k = 1:nk    
            
            % Find indices of rising and falling edges
            a_ind_fall = 1:nn-1; a_ind_fall = a_ind_fall(a_fall(j,k,i,:));
            a_ind_rise = 1:nn-1; a_ind_rise = a_ind_rise(a_rise(j,k,i,:));
            if e.V_shaft(j,k,i,1) == 0; a_ind_rise(1) = []; end
            if e.V_shaft(j,k,i,end) == 0; a_ind_fall(end) = []; end
            
            % Offset indices by probe location
            a_ind = 0.5 * (a_ind_fall + a_ind_rise);
            a_ind = a_ind + round(rate * e.t(j,k,i) / (f * 360));
            a_ind(a_ind < 0 | a_ind > nn) = [];

            % Smooth shaft triggering by polynomial fit of speed
            a = 1:length(a_ind); a_ind = round(polyval(polyfit(a,a_ind,2),a));
            a_ind(a_ind < 1 | a_ind > size(e.V,4)) = [];

            % Trim signal
            a = a_ind(1):a_ind(end);
            V = squeeze(e.V(j,k,i,a)); a_ind = a_ind - a_ind(1) + 1;

            % Re-interpolate hot-wire readings to constant shaft position
            n_points = round(mean(diff(a_ind))); n_revs = length(a_ind)-1;
            ph = interp1(a_ind,0:n_revs,1:a_ind(end),'pchip')';
            ph_const = linspace(0,ph(end),n_points*n_revs+1)';
            V_hw_interp = interp1(ph,V,ph_const,'pchip');

            % Compute ensemble average
            V_ens = mean(reshape(V_hw_interp(1:end-1),[n_points n_revs]),2);
            
            % Compute velocity fluctuation
            V_fluc = V_hw_interp(1:end-1) - repmat(V_ens,[n_revs 1]);
            
            % Reduce to desired fraction of annulus
            V_av(j,k,i,:) = V_ens(q); V_p(j,k,i,:) = V_fluc(q);

%             % Find indices of falling edges
%             a_ind = 1:length(a_log(j,k,i,:)); a_ind = a_ind(a_log(j,k,i,:));
% 
%             % Offset indices by probe location
%             a_ind = a_ind + round(rate * e.t(j,k,i) / (f * 360));
%             a_ind(a_ind < 0 | a_ind > nn) = [];
%             
%             % Remove spurious triggerings of once per rev signals
%             dt = diff(a_ind); a_ind(dt < 0.5 * max(dt)) = [];
% 
%             % Average first nt points
%             for m = 1:length(a_ind)-1
%                 V_av(j,k,i,:) = V_av(j,k,i,:) + e.V(j,k,i,a_ind(m):a_ind(m)+nt-1);
%             end
%             V_av(j,k,i,:) = V_av(j,k,i,:)/(length(a_ind)-1);

            % Calculate disk angular velocity
            w(j,k,i) = rate * 2 * pi / mean(diff(a_ind));

        end
    end
end

% Calculate rotor relative tangential coordinate
t_rel = e.dtime(1:nt) * mean(w(isnan(w(:)) == 0));

% Plot all voltage readings against probe angle to check the fit
if plot_stuff == 1
    figure(); hold on; grid on; box on; xlabel('Probe Angle / ^\circ'); ylabel('Velocity / ms^{-1}');
    plot(reshape(permute(e.Iota,[3 1 2]),[ni nj*nk]),...
        reshape(permute(mean(V_av,4),[3 1 2]),[ni nj*nk]),'k-')
end

% Remove points if angles are skewed
if fit_order == 2
    
end

% Construct angle vector and calculate polynomial regression coefficients
if fit_order == 2
    x = linspace(0,1,ni)'; X = [ones(length(x),1) x x.^2];
elseif fit_order == 4
    x = linspace(0,1,ni)'; X = [ones(length(x),1) x x.^2 x.^3 x.^4];
end
X = (X' * X) \ X'; X = reshape(X,[1 1 size(X,1) 1 size(X,2)]);

% Calculate quadratic coefficients by multiplying regression coefficicients by voltage readings
p = flip(mmat(X,V_av,[3 5]),3);

% Plot polynomial fits
if plot_stuff == 1
    p_av = mean(p,4); ni = 100; x_ch = repmat(reshape(linspace(0,1,ni),[1 1 ni]),[nj nk 1]);
    if fit_order == 2
        V_ch = repmat(p_av(:,:,1),[1 1 ni]) .* x_ch.^2 + repmat(p_av(:,:,2),[1 1 ni]) .* x_ch + ...
            repmat(p_av(:,:,3),[1 1 ni]);
    elseif fit_order == 4
        V_ch = repmat(p_av(:,:,1),[1 1 ni]) .* x_ch.^4 + repmat(p_av(:,:,2),[1 1 ni]) .* x_ch.^3 + ...
            repmat(p_av(:,:,3),[1 1 ni]) .* x_ch.^2 + repmat(p_av(:,:,4),[1 1 ni]) .* x_ch + ...
            repmat(p_av(:,:,5),[1 1 ni]);
    end
    Iota_ch = x_ch .* repmat(e.Iota(:,:,end) - e.Iota(:,:,1),[1 1 ni]) + repmat(e.Iota(:,:,1),[1 1 ni]);
    plot(reshape(permute(Iota_ch,[3 1 2]),[ni nj*nk]),reshape(permute(V_ch,[3 1 2]),[ni nj*nk]),'r-')
end

% Calculate flow angle at turning point of quadratic
if fit_order == 2
    x_max = -p(:,:,2,:) ./ (2 * p(:,:,1,:));
    V = p(:,:,1,:) .* x_max.^2 + p(:,:,2,:) .* x_max + p(:,:,3,:);
end

% Calculate flow angle at turning point of quartic using Newton-Raphson method
if fit_order == 4
    [~,i] = max(V_av,[],3); x_max = x(i);
    for n = 1:10
        y = 4 * p(:,:,1,:) .* x_max.^3 + 3 * p(:,:,2,:) .* x_max.^2 + 2 * p(:,:,3,:) .* x_max + ...
            p(:,:,4,:);
        dydx = 12 * p(:,:,1,:) .* x_max.^2 + 6 * p(:,:,2,:) .* x_max + 2 * p(:,:,3,:);
        x_max = x_max - y ./ dydx;
    end
    V = p(:,:,1,:) .* x_max.^4 + p(:,:,2,:) .* x_max.^3 + p(:,:,3,:) .* x_max.^2 + ...
        p(:,:,4,:) .* x_max + p(:,:,5,:);
end

% Recalculate flow angle from non-dimensional maximum, only works on equally linearly spaced matrices
Alpha = x_max .* repmat(e.Iota(:,:,end) - e.Iota(:,:,1),[1 1 1 nt]) + repmat(e.Iota(:,:,1),[1 1 1 nt]);

% Squeeze dimensions
Alpha = real(reshape(Alpha,[nj nk nt])); V = real(reshape(V,[nj nk nt])); w = repmat(mean(w,3),[1 1 nt]);

% Calculate meridional velocity
Vm = V .* cosd(Alpha);

% Calculate cartesian plotting coordinates
rt = (t_rel * squeeze(e.r(:,1,1))')'; r = repmat(e.r(:,1),[1 nt]); t = rt./r;
[y,z] = pol2cart(t - mean(t(:)),r); x = repmat(e.x(:,1),[1 nt]);

% Record processed data
g.Vm = Vm; g.V = V; g.Alpha = Alpha; g.w = w; g.V_p = V_p;
g.x = x; g.r = r; g.rt = rt; g.t = t; g.z = z; g.y = y;


end
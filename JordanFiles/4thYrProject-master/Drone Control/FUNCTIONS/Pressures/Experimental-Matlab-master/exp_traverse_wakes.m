function tr_wake = exp_traverse_wakes(s,search,park,probe,N,des,err,r_mid)
% EXP_TRAVERSE_WAKES  Survey downstream of stators to find wakes
%
%   tr_wake = EXP_TRAVERSE_WAKES(s,search,park,probe)
%
%   s - data structure with open instrument connections
%   search - data structure with limits of movement for wake search
%   park - data structure of stepper offsets to convert from rig coordinates to steps
%   probe - current probe for offsets


%% Create coarse search grid to move probe over

% Matrices of linearly spaced r and theta coordinates
r = repmat(reshape(linspace(search.r_min,search.r_max,search.nj),[search.nj 1]),[1 search.nk]);
t = repmat(reshape(linspace(search.t_min,search.t_max,search.nk),[1 search.nk]),[search.nj 1]);

% Probe angle
Iota = search.Iota .* ones(search.nj,search.nk);

% Check type of traverse gear
if strcmp(park.type,'XY') == 1
    
    % Calculate cartesian coordinates
    [y,z] = pol2cart(deg2rad(t),r);
    
    % Calculate traverse grid from cartesian coordinates
    q = round(cat(3,(y + probe.r_off) * park.d(1) + park.c(1),z * park.d(2) + park.c(2),...
        (Iota + probe.Iota_off) * park.d(3) + park.c(3)));    
    
elseif strcmp(park.type,'RT') == 1
    
    % Calculate traverse grid from polars directly
    q = round(cat(3,(r + probe.r_off) * park.d(1) + park.c(1),t * park.d(2) + park.c(2),...
        (Iota + probe.Iota_off) * park.d(3) + park.c(3)));

end

% Calculate traverse grid in steps
% q = round(cat(3,(r + probe.r_off) * park.d(1) + park.c(1),t*park.d(2) + park.c(2),...
%     Iota*park.d(3) + park.c(3)));

% Create index vectors for tangential or radial stepping
J = repmat(reshape(search.nj:-1:1,[search.nj 1]),[1 search.nk]); 
K = repmat(reshape(1:search.nk,[1 search.nk]),[search.nj 1]); 
if strcmp(park.type,'XY') == 1
    J = reshape(J,[],1); K = reshape(K,[],1);    
else
    J = reshape(J',[],1); K = reshape(K',[],1);
end
 
%% Initialise monitoring window

% Open figure window
figure(); hold on; box on; axis equal;

% Plot mesh of traverse grid
[x,y] = pol2cart(2*pi*(t+90)/360,r);
mesh(x,y,zeros(search.nj,search.nk),zeros(search.nj,search.nk,3));

% Calculate offset grid
x_off = 0.25 * (x(1:end-1,1:end-1) + x(1:end-1,2:end) + x(2:end,1:end-1) + x(2:end,2:end));
y_off = 0.25 * (y(1:end-1,1:end-1) + y(1:end-1,2:end) + y(2:end,1:end-1) + y(2:end,2:end));
x_off = [x(1,1) 0.5 * (x(1,1:end-1) + x(1,2:end)) x(1,end) ; ...
    0.5 * (x(1:end-1,1) + x(2:end,1)) x_off 0.5 * (x(1:end-1,end) + x(2:end,end)) ;...
    x(end,1) 0.5 * (x(end,1:end-1) + x(end,2:end)) x(end,end)];
y_off = [y(1,1) 0.5 * (y(1,1:end-1) + y(1,2:end)) y(1,end) ; ...
    0.5 * (y(1:end-1,1) + y(2:end,1)) y_off 0.5 * (y(1:end-1,end) + y(2:end,end)) ;...
    y(end,1) 0.5 * (y(end,1:end-1) + y(end,2:end)) y(end,end)];


%% Survey exit and record centre hole pressures

% Set rig operating point for the first time
uiwait(msgbox('Start Rig & Set Traverse Support'));
exp_rig_set(s,N,des,err,r_mid);

% Move to first location
j = J(1); k = K(1);
exp_stepper_coord(s.park,[0 0 q(j,k,3)],park.limit);
exp_stepper_coord(s.park,[0 q(j,k,2) q(j,k,3)],park.limit);
exp_stepper_coord(s.park,q(j,k,:),park.limit);

% Initialise pressure array
P = zeros(search.nj,search.nk);

% Loop over all locations and record centre hole pressure
for n = 1:search.nj*search.nk

    % Get current indices
    j = J(n); k = K(n);
    
    % Move to next location
    if strcmp(park.type,'XY') == 1
        
        % Radial motion for cartesian traverse gears
        exp_stepper_move(s.park,1,q(j,k,1),park.limit)
        exp_stepper_position(s.park,1,q(j,k,1))
        exp_stepper_move(s.park,2,q(j,k,2),park.limit)
        exp_stepper_position(s.park,2,q(j,k,2))

    else
        
        % Tangential motion for polar traverse gears
        if k == 1
            exp_stepper_coord(s.park,[q(j,k,1) q(j,k,2)-park.q_comp q(j,k,3)],park.limit);
        end
        exp_stepper_move(s.park,2,q(j,k,2),park.limit)
        exp_stepper_position(s.park,2,q(j,k,2))
    end
    
    % Settle probe
    if isfield(probe,'time_cen') == 1
        pause(probe.time_cen);
    else
        pause(probe.time);
    end
    
    % Record centre hole pressure
    P_temp = exp_dsa_read(s.dsa);
    P(j,k) = P_temp(N.P_FHP(1));
    
    % Plot patch of pressure
    patch([x_off(J(n),K(n)) ; x_off(J(n)+1,K(n)) ; x_off(J(n),K(n)+1)],...
        [y_off(J(n),K(n)) ; y_off(J(n)+1,K(n)) ; y_off(J(n),K(n)+1)],...
        ones(3,1)*P(j,k),'EdgeColor','none','FaceColor','flat')
    patch([x_off(J(n)+1,K(n)+1) ; x_off(J(n)+1,K(n)) ; x_off(J(n),K(n)+1)],...
        [y_off(J(n)+1,K(n)+1) ; y_off(J(n)+1,K(n)) ; y_off(J(n),K(n)+1)],...
        ones(3,1)*P(j,k),'EdgeColor','none','FaceColor','flat')
end

% Drive motors back home
exp_stepper_coord(s.park,[0 0 0],park.limit);

%% Fit gaussians to measured wakes to find centres

figure(); hold on; grid on; box on;

% Loop over each slice, and find centre of wake
tr_wake = [];
for j = 1:size(r,1);
    t_temp = t(j,:).'; r_temp = r(j,:).';
    P_temp = P(j,:).'; P_temp = max(P_temp) - P_temp;
    [~,k_max] = max(P_temp);
    
    % Fit a gaussian pulse to the wake
    ft = fittype( 'gauss1' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Normalize = 'on';
    k = k_max-3:k_max+3; k(k > length(P_temp) | k < 1) = [];
    f = fit( t_temp(k), P_temp(k), ft, opts );    
    
    % Plot fit
    plot(t_temp,P_temp,'k.-')
    t_fit = linspace(min(t_temp),max(t_temp),1000);
    P_fit = feval(f,t_fit);
    plot(t_fit,P_fit,'r-')
    
    % Return wake locations
    [~,i] = max(P_fit);
    tr_wake = [tr_wake ; t_fit(i) mean(r_temp)]; 
end


end
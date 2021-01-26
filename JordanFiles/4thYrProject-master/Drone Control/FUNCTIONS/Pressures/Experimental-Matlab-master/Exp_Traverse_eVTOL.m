%   Exp_Traverse                               
%                               
%   A script to run a series of traverses on any rig
%
%   Written by James Taylor                
%   October 2016
%   Edited by Megan McEveley for eVTOL traverse
%   August 2019
%
%   Choose 3 types of calls with run_type - 'grid', 'trav', 'off' or 'full'

% Default to a complete traverse
if exist('run_type','var') == 0
    run_type = 'full';
end 

%% Determine traverse type and load calibration map

if strcmp(probe_name(1:3),'FHP') == 1 || strcmp(probe_name(1),'P') == 1
    
    % Five hole probe traverse
%     trav_type = 'FH'; load([dr.exp calib_name]);
    trav_type = 'FH'; load('FHP_M_40_fine_01-Aug-2019_Bingo_1.mat');
    pxie.nchan = 1; pxie.chan = N.V_shaft; pxie.time = 0.2;
    
elseif strcmp(probe_name(1:2),'HW') == 1
    
    % Hot wire traverse
    trav_type = 'HW'; load([dr.exp calib_name]);
    pxie.nchan = 2; pxie.chan = [N.V_shaft N.V_hw N.V_hw_av]; 
    pxie.time = probe.(probe_name).time;
    
elseif strcmp(probe_name(1:2),'To') == 1
    
    % Total temperature traverse
    trav_type = 'To'; tc08.nchan = tc08.nchan + 1;
    
end

%% Initialise all instruments

% Check whether current call is for a dry run
if strcmp(run_type,'grid') == 0

    % Connect and power up stepper motors
   s.park = exp_stepper_initialise(park,0);
    
    % Power up other instruments if full traverse
    if strcmp(run_type,'off') == 1 || strcmp(run_type,'full') == 1
    
        % Initialise and zero DSA
        s.dsa = exp_dsa_initialise(dsa);

        % Initialise thermocouple logger
%         s.tc08 = exp_tc08_initialise(tc08);

        % Initialise voltage logger
        s.pxie = exp_pxie_initialise(pxie);
        
        % Initialise RTD logger
%         s.rtd = exp_pxie_initialise(rtd);

        % Point to low speed lab PTH sensor
%         if strcmp(N_rig.PTH,'Vaisala') == 1
% %             if isempty(timerfind) == 1; PTUstart(ptu.loc); end
%             s.ptu = ptu;
%         end
    end

    % Home all axes simultaneously
    exp_stepper_home(s.park,park.home,1)
end

% Number of pressure and temperature channels
if exist('s','var') ~= 0 && isfield(s,'rtd') == 0
    T_chan = length(tc08.chan);
else
    T_chan = length(tc08.chan) + length(rtd.chan);
end
P_chan = dsa.nchan;

% Global voltage and temperature NI-DAQ variables
global V time_v T time_t
V = []; time_v = []; T = []; time_t = []; 

% Start a listeners with global NI-DAQ variables
if strcmp(run_type,'off') == 1 || strcmp(run_type,'full') == 1
    s.pxie.NotifyWhenDataAvailableExceeds = s.pxie.NumberOfScans;
%     s.rtd.NotifyWhenDataAvailableExceeds = s.rtd.NumberOfScans;
    addlistener(s.pxie,'DataAvailable',@exp_pxie_background);
%     addlistener(s.rtd,'DataAvailable',@exp_pxie_background);
end

% Narrow limit range based on traverse location
% if isstruct(park.limit) == 1
%     park.limit = park.limit.(['T' num2str(st)]);
% end

% % Alter limit range based upon probe offset
% park.limit(1,:) = park.limit(1,:) + park.d(1) * probe.(probe_name).r_off; 
% park.limit(3,:) = park.limit(3,:) + park.d(3) * (0 + probe.(probe_name).Iota_off); 
% park.limit = round(park.limit);


%% Set up traverse grid and point order

% Calculate nearest flow coefficients from flow function
if isfield(des,'phi') == 0
    T_sl = 293.15; P_sl = 101325; ro_sl = P_sl / (air.R * T_sl);
    Vx = des.mrtp * P_sl / (ro_sl * (air.cp * T_sl)^0.5);
    U = r_mid * des.nrt * (air.cp * T_sl)^0.5;
    phi_loop = round(100*Vx./U)/100;
else
    phi_loop = des.phi;
end

% Loop over all desired flow coefficients
nk_bl = nk;
for a = 1:length(phi_loop)

% Clear the previous experiment
if exist('e','var') ~= 0; clear e; end;
    
% Get current design point
if isfield(des,'phi') == 0
    des_current.mrtp = des.mrtp(a);
    des_current.nrt = des.nrt(a);
else
    des_current.phi = des.phi(a);
end
    
% Assemble r coordinates for traverse clustering on endwalls if necessary
nk = nk_bl;
if isempty(r_fac) == 0 && exist('r','var') == 0
    if numel(r_fac) == 1; r_fac = [r_fac r_fac]; end
    r = hyperbolic_bunch(nj,r_fac(2)/nj,r_fac(1)/nj) * (r_max - r_min) + r_min;
elseif exist('r','var') == 0
    r = linspace(r_min,r_max,nj);
end
r = repmat(reshape(r,[nj 1]),[1 nk]);

% Determine theta clustering as a function of r
if exist('t_fac','var') ~= 0 && length(t_fac) > 1
%     t_facs = interp1([r_hub r_mid r_cas],t_fac([1 2 1]),r(:,1),'pchip');
    t_facs = interp1([0 0.25 0.75 1] * (r_cas - r_hub) + r_hub,t_fac([1 2 2 1]),r(:,1),'pchip');
elseif exist('t_fac','var') ~= 0
    t_facs = t_fac * ones(nj,1);
end

% Determine which distribution of tangential coordinates to use
if isfloat(t_clust) == 1
    
    % Hyperbolically clustered on one tangential location
    t = []; n_pitch = length(t_clust);
    t_pitch = (t_max - t_min) / n_pitch;
    for n = 1:n_pitch
        t_1 = t_min + (n-1) * t_pitch;
        t_2 = t_min + n * t_pitch;

        n_1 = round((t_clust(n) - t_1) * nk / ((t_2 - t_1) * n_pitch));
        n_2 = round(nk / n_pitch) - n_1 + 1;

        t_temp_1 = hyperbolic_bunch(n_1,t_fac/n_1,1/n_1) * (t_clust(n) - t_1) + t_1;
        t_temp_2 = hyperbolic_bunch(n_2,1/n_2,t_fac/n_2) * (t_2 - t_clust(n)) + t_clust(n);

        t = [t t_temp_1 t_temp_2];
    end
    t(diff(t) < 1e-6) = [];
    t = repmat(reshape(t,[1 nk]),[nj 1]);
    
elseif strcmp(t_clust,'geom') == 1
    
    % Check if wakes have already been found for the same geometry
    tr_wake = [];
    A = dir([dr.exp geom_name '*Ph0' strrep(num2str(phi_loop(a)*100),'.','')...
        '*T' num2str(st) '_*.mat']);
    if isempty(A) == 0 && search.override == 0
        i = 0; time_temp = 0; 
        for n = 1:length(A)
            if A(n).datenum > time_temp
                time_temp = A(n).datenum;
                i = n;
            end
        end
        load([dr.exp A(i).name]); 
        tr_wake = e.tr_wake; clear('e');
    end
    
    % Check if wake location is specified explicitly in a cell array
    if exist('tr_wakes','var') ~= 0; tr_wake = tr_wakes{a}; end
    
    % Search for wakes with quick probe traverse if required
    if isempty(tr_wake) == 1 || search.override == 1
        tr_wake = exp_traverse_wakes(s,search,park,probe.(probe_name),N_rig,des_current,err,r_mid);
    end
    e.tr_wake = tr_wake; t_pitch = search.t_pitch;
    
    % Fit polynomial through found wakes and calculate tangential spacing
    p = polyfit(tr_wake(:,2),tr_wake(:,1),2); clear tr_wake;
    t = zeros(nj,nk);
    for j = 1:nj
        t_mid = polyval(p,r(j,1)); 
        t_max = t_mid + t_pitch/2;
        t_min = t_mid - t_pitch/2;
        t_temp = [hyperbolic_bunch(round((nk+1)/2),t_facs(j),1) * (t_mid - t_min) + t_min ...
            hyperbolic_bunch(round((nk+1)/2),1,t_facs(j)) * (t_max - t_mid) + t_mid];
        t_temp(diff(t_temp) < 1e-5) = [];
        t(j,:) = t_temp;
    end
    
elseif strcmp(t_clust,'linear') == 1
    
    % Uniform distribution
    t = repmat(reshape(linspace(t_min,t_max,nk),[1 nk]),[nj 1]);
    
end

% Repeat t and r coordinates for multiple passage traverses
rs = r; ts = t; nk = length(bl) * (nk - 1) + 1;
for n = 1:length(bl)-1
    t = [t ts(:,2:end) - (bl(n)-2) * t_pitch];
    r = [r rs(:,2:end)];
end

% % Define angle matrix
% if ni == 1
%     if isnumeric(Iota) == 0
%         Iota = Iota(r,t);
%     elseif numel(Iota) == 1
%         Iota = Iota .* ones(nj,nk);
%     end
%     dim_coord = 3;
% else
%     if numel(Iota_min) == 1
%         Iota = repmat(reshape(linspace(Iota_min,Iota_max,ni),[1 1 ni]),[nj nk 1]);
%     else
%         si = repmat(reshape(linspace(0,1,ni),[1 1 ni]),[nj nk 1]); 
%         Iota = si .* repmat(Iota_max - Iota_min,[1 nk ni]) + repmat(Iota_min,[1 nk ni]);
%     end
%     r = repmat(r,[1 1 ni]); t = repmat(t,[1 1 ni]);
%     dim_coord = 4;
% end

% Move tangential coordinates if limit is triggered
t_lim = (park.limit(2,2) - park.c(2)) / park.d(2);
if max(t(:)) > t_lim
    k_lim = find(max(t,[],1) > t_lim,1) - 1;
    t = [t(:,k_lim:end) - length(bl) * t_pitch t(:,2:k_lim)];
end

% Check type of traverse gear
if strcmp(park.type,'XY') == 1
    
%     % Calculate cartesian coordinates
%     [y,z] = pol2cart(deg2rad(t),r);
%     
%     % Calculate traverse grid from cartesian coordinates
%     q = round(cat(dim_coord,(y + probe.(probe_name).r_off) * park.d(1) + park.c(1),...
%         z * park.d(2) + park.c(2),...
%         (Iota + probe.(probe_name).Iota_off) * park.d(3) + park.c(3)));    
    
    % Calculate cartesian coordinates
    x_min = 0, x_max = 0.275;
    y_min = 0, y_max = 0.300;
    nx = 20, ny = 20;
    x = linspace(x_min, x_max, n);
    y = linspace(y_min, y_max, n);
   
    % Calulcate traverse grid from cartesian coordinates
    dim_coord = 2;
    q = round(cat(dim_coord, x * park.d(3), y * park.d(4)));

elseif strcmp(park.type,'RT') == 1
    
    % Calculate traverse grid from polars directly
    q = round(cat(dim_coord,(r + probe.(probe_name).r_off) * park.d(1) + park.c(1),...
        t * park.d(2) + park.c(2),...
        (Iota + probe.(probe_name).Iota_off) * park.d(3) + park.c(3)));

end

% Create index vectors for either circumferential or radial stepping
if strcmp(dir_type,'radial') == 1
    J = repmat(reshape(1:nj,[nj 1]),[1 nk]); 
    K = repmat(reshape(1:nk,[1 nk]),[nj 1]); 
    J = reshape(J,[],1); K = reshape(K,[],1);
elseif strcmp(dir_type,'tang') == 1
    J = repmat(reshape(nj:-1:1,[nj 1]),[1 nk]); 
    K = repmat(reshape(1:nk,[1 nk]),[nj 1]); 
    J = reshape(permute(J,[2 1 3]),[],1); K = reshape(permute(K,[2 1 3]),[],1);
end

% Display minimum and maximum coordinates in each dimension
if ni == 1
    disp(['Max Steps ' num2str(reshape(max(max(q,[],1),[],2),[1 3]))])
    disp(['Min Steps ' num2str(reshape(min(min(q,[],1),[],2),[1 3]))])
else
    disp(['Max Steps ' num2str(reshape(max(max(max(q,[],1),[],2),[],3),[1 3]))])
    disp(['Min Steps ' num2str(reshape(min(min(min(q,[],1),[],2),[],3),[1 3]))])
end

% Open a figure window to display the grid
if size(q,1) > 1 && size(q,2) > 1 
    figure('position',[471 85 788 820]); subplot(2,1,1); hold on; grid on; box on; axis equal;
    xlabel('Z / m'); ylabel('Y / m')

    % Plot the mesh of the traverse grid
    if strcmp(park.type,'RT') == 1
        [y,z] = pol2cart(deg2rad(t),r);
    end
    mesh(z(:,:,1), y(:,:,1), zeros(size(y(:,:,1))), zeros([size(y,1) size(y,2) 3]))

    % Display the step counts
    subplot(2,1,2); hold on; grid on; box on; axis equal;
    xlabel('Axis 2 / steps'); ylabel('Axis 1 / steps');
    if size(q,4) == 1
        mesh(q(:,:,2), q(:,:,1), zeros(size(r)), zeros([size(r) 3]))
    else
        mesh(q(:,:,1,2), q(:,:,1,1), zeros(size(r(:,:,1))), zeros([size(r,1) size(r,2) 3]))
    end
end

% Check if steps exceed the limits
q_lim = reshape(q,[],3); q_lim = [min(q_lim,[],1)' max(q_lim,[],1)'];
if any(q_lim(:,1) < park.limit(:,1)) || any(q_lim(:,2) > park.limit(:,2))
    disp(['Min Limits ' num2str(park.limit(:,1)')]);
    disp(['Max Limits ' num2str(park.limit(:,2)')]);
    error('Stepper coordinates greater than traverse limits');
end

% Stop traverse execution if current call is for a dry run
if strcmp(run_type,'grid') == 1
    return
end


%% Initialise data structure

% Initialise raw pressure matrices
if strcmp(trav_type,'FH') == 1
    e.P = nan(nj,nk,P_chan);
    e.time = nan(nj,nk,6);
elseif strcmp(trav_type,'HW') == 1
    e.V_av = nan(nj,nk,ni);
    e.time = nan(nj,nk,ni,6);
end

% Initialise raw rig operating point information
rignames = [{'P' 'T' 'Pa' 'Ta' 'Ha' 'w' 'time'} fieldnames(des)'];
for v = 1:length(rignames)
    e.rig.(rignames{v}) = [];
end

% Save grid steps and real coordinates
e.q = q; e.Iota = Iota; e.J = J; e.K = K; e.r = r; e.t = t;

% Save configuation information
if isfield(des_current,'phi') == 1
    e.phi = des_current.phi;
else
    e.mrtp = des_current.mrtp;
    e.nrt = des_current.nrt;
end
e.N = N; e.geom_name = geom_name;

% Save probe information and calibration map
e.probe_name = probe_name; e.probe = probe.(probe_name);

if strcmp(trav_type,'FH') == 1
    e.calib_name = calib_name; e.c = c; 
end


%% Open traverse monitoring matlab session

% Traverse name
if ni == 1
    trav_suffix = [num2str(nk) 'x' num2str(nj)];
else
    trav_suffix = [num2str(nk) 'x' num2str(nj) 'x' num2str(ni)];
end

% Determine new filename - no overwriting
n = 1;
file_prefix = [geom_name '_' trav_suffix '_Ph0' strrep(num2str(phi_loop(a)*100),'.','') '_' date...
    '_T' num2str(st) '_' trav_type];
file_prefix = strrep(file_prefix,'.','p');
while exist([dr.exp file_prefix '_' num2str(n) '.mat'],'file') ~= 0
    n = n+1;
end
filename = [dr.exp file_prefix '_' num2str(n)];

% Save initial traverse data
save([filename '.mat'],'e')

% Open new matlab instance with monitoring function
save('current_traverse.mat')
system(['matlab -nosplash -nodesktop -r exp_traverse_monitor(''' filename ''') &'])

% Save hot wire raw voltage files
if strcmp(trav_type,'HW') == 1
    nt = pxie.rate*pxie.time+1;
    V_hw = nan(nj, nk, ni, nt, 'single');
    V_shaft = true(nj, nk, ni, nt);
end


%% Start rig and set flow coefficient

% Record zero flow voltage
if strcmp(trav_type,'HW') == 1
    uiwait(msgbox('Rig off zero flow voltage'));
    [V_temp,t_temp] = exp_pxie_read(s.pxie);
    disp(mean(V_temp));
    e.V0_start = single(V_temp(:,3));
end

% Set rig operating point for the first time
if strcmp(run_type,'full') == 1
    uiwait(msgbox('Start Rig & Set Traverse Support'));
    [e_rig,on_point] = exp_rig_set(s,N_rig,des_current,err,r_mid);
else
    on_point = 1;
end


%% Execute traverse

% Move probe to starting position in reverse axis order
j = J(1); k = K(1);
if ni == 1
    q_temp = squeeze(q(j,k,:));
else
    i = 1;
    q_temp = squeeze(q(j,k,i,:));
end
exp_stepper_coord(s.park,[0 0 q_temp(3)],park.limit);
exp_stepper_coord(s.park,[0 q_temp(2) q_temp(3)],park.limit);
exp_stepper_coord(s.park,q_temp,park.limit);

% Loop over all points in grid, moving in either tangential or radial
% direction, record and process traverse data
for n = 1:nj*nk
    
    % Get current coordinates
    j = J(n); k = K(n); q_old = q_temp;
    if ni == 1
        q_temp = squeeze(q(j,k,:));
    else
        i = 1;
        q_temp = squeeze(q(j,k,i,:));
    end
    
    % Display current indices
    disp(['Measuring Point J = ' num2str(j) ' K = ' num2str(k)]);
    
    % Tangential motion
    if strcmp(dir_type,'tang') == 1
        
        % Move to new radius and extend tangential direction to remove backlash
        if k == 1
            
            % If using a hotwire drive to mid radius before tangential positioning
            if strcmp(trav_type,'HW') == 1
                exp_stepper_coord(s.park,[0 q_temp(2)-park.q_comp q_temp(3)],park.limit);
            end
            
            exp_stepper_coord(s.park,[q_temp(1) q_temp(2)-park.q_comp q_temp(3)],park.limit);
        end
        
        % Position HW at mid radius first
        if strcmp(trav_type,'HW') == 1
            exp_stepper_move(s.park,1,0,park.limit)
            exp_stepper_position(s.park,1,0)
        end
        
        % Move to new tangential coordinates
        exp_stepper_move(s.park,2,q_temp(2),park.limit)
        
        % Query tangential position
        exp_stepper_position(s.park,2,q_temp(2))
        
        % Reposition HW at intended radius
        if strcmp(trav_type,'HW') == 1
            exp_stepper_move(s.park,1,q_temp(1),park.limit)
            exp_stepper_position(s.park,1,q_temp(1))
        end
    end
    
    % Radial motion
    if strcmp(dir_type,'radial') == 1
        
        % Move to new tangential location, repositioning at mid radius
        if j == 1
            exp_stepper_coord(s.park,[0 q_temp(2) q_temp(3)],park.limit);
        end
        
        % Move to new radial coordinates
        exp_stepper_move(s.park,1,q_temp(1),park.limit)
        
        % Move tangential or rotary directions if required
        if ni == 1 && q_temp(3) ~= q_old(3); exp_stepper_move(s.park,3,q_temp(3),park.limit); end
        if q_temp(2) ~= q_old(2); exp_stepper_move(s.park,2,q_temp(2),park.limit); end
        
        % Query radial position
        exp_stepper_position(s.park,1,q_temp(1))
        
        % Query tangential and rotary positions
        if ni == 1 && q_temp(3) ~= q_old(3); exp_stepper_position(s.park,3,q_temp(3)); end
        if q_temp(2) ~= q_old(2); exp_stepper_position(s.park,2,q_temp(2)); end
        
    end
    
    % Pause double for settling at first point in cycle
    if (strcmp(dir_type,'radial') == 1 && j == 1) || (strcmp(dir_type,'tang') == 1 && k == 1) && ...
            strcmp(trav_type,'HW') == 0
        pause(probe.(probe_name).time);
    end
    
    % Loop over probe angles if specified
    for i = 1:ni
        
        % Move rotary table
        if ni > 1
            exp_stepper_move(s.park,3,q(j,k,i,3),park.limit)
            exp_stepper_position(s.park,3,q(j,k,i,3))
        end
        
        % Check if full traverse with instrument logging is specified
        if strcmp(run_type,'off') == 1 || strcmp(run_type,'full') == 1
        
            % Check if rig is on point while probe is settling
            tic; 
            if on_point == 0 || mod(n,n_check) == 0 || n == 1 || n == length(J) ...
                    || strcmp(trav_type,'HW') == 1
                [e.rig,on_point,e_temp] = exp_rig_monitor(s,e.rig,N_rig,des_current,err,r_mid); 
            end
            
            % Reset rig operating point if not on
            if on_point == 0 && i == 1 && strcmp(run_type,'full') == 1
                exp_rig_set(s,N_rig,des_current,err,r_mid);
                [e.rig,on_point,e_temp] = exp_rig_monitor(s,e.rig,N_rig,des_current,err,r_mid); 
            end

            % Pause for remaining probe settling
            settled_time = toc;
            pause(max(probe.(probe_name).time - settled_time,0));

            % Check the type of instrument used for the traverse
            if strcmp(trav_type,'FH') == 1

                % Log DSA pressures
                e.P(j,k,:) = exp_dsa_read(s.dsa);

            elseif strcmp(trav_type,'HW') == 1

                % Record hotwire voltages from rig monitoring data logging
                V_hw(j,k,i,:) = e_temp.V(:,N.V_hw);
                e.V_av(j,k,i) = mean(e_temp.V(:,N.V_hw_av));
                V_shaft(j,k,i,e_temp.V(:,N.V_shaft) < pxie.V_thresh) = 0;
                e.dtime = e_temp.dtime;
                disp(['Vav = ' num2str(e.V_av(j,k,i))]);
                
                % Check ADC is not clipping because gain is too high
                V_max = max(abs(V_hw(j,k,i,:)));
                if V_max > 10.9
                    disp(['ADC is clipping, V_max = ' num2str(V_max)]);
                end

            end

            % Log time
            if ni == 1
                e.time(j,k,:) = clock;
            else
                e.time(j,k,i,:) = clock;
            end

            % Save data to disk
            save([filename '.mat'],'e')
        end
    end
    
    % Save a row of hotwire voltages to file in case of a crash
    if strcmp(trav_type,'HW') == 1
        if k == nk && strcmp(dir_type,'tang') == 1
            e_temp.V_hw = V_hw(j,1:nk,1:ni,1:nt); 
            e_temp.V_shaft = V_shaft(j,1:nk,1:ni,1:nt);
            save(['V_hw_' num2str(j) '.mat'],'e_temp')
        elseif j == nj && strcmp(dir_type,'radial') == 1
            e_temp.V_hw = V_hw(1:nj,k,1:ni,1:nt); 
            e_temp.V_shaft = V_shaft(1:nj,k,1:ni,1:nt);
            save(['V_hw_' num2str(k) '.mat'],'e_temp')
        end
    end
end

% Return motors to zeros and re-home
exp_stepper_coord(s.park,[0 0 0],park.limit);
if length(phi_loop) > 1
    exp_stepper_home(s.park,park.home,1)
end

% Turn off motors and quit the script if grid test is complete
if strcmp(run_type,'trav') == 1
    for n = 1:park.nchan
        exp_stepper_command(s.park,n,'OFF');
        exp_stepper_command(s.park,n,'E0');
    end
    fclose(s.park); delete(s.park); s = rmfield(s,'park');
    return;
end

% Save all data
if strcmp(trav_type,'HW') == 1
    e.V_shaft = V_shaft; clear V_shaft;
    e.V_hw = V_hw; clear V_hw;
    
    uiwait(msgbox('Rig off zero flow voltage'));
    [V_temp,t_temp] = exp_pxie_read(s.pxie);
    e.V0_end = single(V_temp(:,3));
    disp(mean(V_temp));
end
save([filename '.mat'],'e','-v7.3')

% Finish looping over flow coefficients
end


%% Stop rig and close all instruments down

% Turn off the rig if it is running
if strcmp(run_type,'full') == 1
    uiwait(msgbox('Stop Rig'));
end

% Disconnect dsa
if iscell(s.dsa) == 1
    for n = 1:length(s.dsa); fclose(s.dsa{n}); end
    s = rmfield(s,'dsa');
else
    fclose(s.dsa); delete(s.dsa); s = rmfield(s,'dsa');
end

% Close thermocouple sensor
calllib('usbtc08','usb_tc08_close_unit',s.tc08.h);
s = rmfield(s,'tc08');

% Disconnect motors
for n = 1:park.nchan
    exp_stepper_command(s.park,n,'OFF');
    exp_stepper_command(s.park,n,'E0');
end
fclose(s.park); delete(s.park); s = rmfield(s,'park');

% Close pxie connection
delete(s.pxie);

% Stop PTU logger
if strcmp(N_rig.PTH,'Vaisala') == 1; delete(timerfind); end

% Clear up connection structure
clear('s'); instrreset; delete(instrfindall);


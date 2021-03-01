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
    
% Five hole probe traverse
trav_type = 'FH'; load('FHP_M_40_fine_01-Aug-2019_Bingo_1.mat');
pxie.nchan = 1; pxie.chan = N.V_shaft; pxie.time = 0.2;


%% Initialise all instruments

% Check whether current call is for a dry run
if strcmp(run_type,'grid') == 0

    % Connect and power up stepper motors
   s.park = exp_stepper_initialise(park,0);
    
    % Power up other instruments if full traverse
    if strcmp(run_type,'off') == 1 || strcmp(run_type,'full') == 1
    
        % Initialise and zero DSA
        s.dsa = exp_dsa_initialise(dsa);

%         % Initialise thermocouple logger
%         s.tc08 = exp_tc08_initialise(tc08);
%         
        % Initialise voltage logger
        s.pxie = exp_pxie_initialise(pxie);
        
    end

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


%% Set up traverse grid and point order

% % Calculate cartesian coordinates
% x_min = 0; x_max = ;
% y_min = 0; y_max = 0.275;
% ni = 20; nj = 20;
% x = repmat(linspace(x_min, x_max, ni),[nj 1]);
% y = repmat(linspace(y_min, y_max, nj)',[1 ni]);

% Assemble r coordinates for traverse clustering on endwalls if necessary
if isempty(r_fac) == 0 %&& exist('r','var') == 0
    if numel(r_fac) == 1; r_fac = [r_fac r_fac]; end
    r = hyperbolic_bunch(nj,r_fac(2)/nj,r_fac(1)/nj) * (r_max - r_min) + r_min;
else %if exist('r','var') == 0
    r = linspace(r_min,r_max,nj);
end
r = repmat(reshape(r,[nj 1]),[1 nk]);

t = repmat(reshape(linspace(t_min,t_max,nk),[1 nk]),[nj 1]);

% Define angle matrix
if ni == 1
    if isnumeric(Iota) == 0
        Iota = Iota(r,t);
    elseif numel(Iota) == 1
        Iota = Iota .* ones(nj,nk);
    end
    dim_coord = 3;
else
    if numel(Iota_min) == 1
        Iota = repmat(reshape(linspace(Iota_min,Iota_max,ni),[1 1 ni]),[nj nk 1]);
    else
        si = repmat(reshape(linspace(0,1,ni),[1 1 ni]),[nj nk 1]); 
        Iota = si .* repmat(Iota_max - Iota_min,[1 nk ni]) + repmat(Iota_min,[1 nk ni]);
    end
    r = repmat(r,[1 1 ni]); t = repmat(t,[1 1 ni]);
    dim_coord = 4;
end


% Check type of traverse gear
if strcmp(park.type,'XY') == 1
    
    % Calculate cartesian coordinates
    [y,z] = pol2cart(deg2rad(t),r);
    
%     Calculate traverse grid from cartesian coordinates
    q = round(cat(dim_coord,(y + probe.(probe_name).r_off) * park.d(1) + park.c(1),...
        z * park.d(2) + park.c(2),...
        (Iota + probe.(probe_name).Iota_off) * park.d(3) + park.c(3)));    
    
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
    
end
% Check if steps exceed the limits
q_lim = reshape(q,[],4); q_lim = [min(q_lim,[],1)' max(q_lim,[],1)'];
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
    e.P = nan(ni,nj,P_chan);
    e.time = nan(ni,nj,6);
    e.rpm = nan(ni, nj, 2);
%     e.temp = nan(ni , nj, T_chan);
elseif strcmp(trav_type,'HW') == 1
    e.V_av = nan(ni,nj,ni);
    e.time = nan(ni,nj,ni,6);
end

% Initialise raw rig operating point information
rignames = [{'P' 'T' 'Pa' 'Ta' 'Ha' 'w' 'time'} fieldnames(des)'];
for v = 1:length(rignames)
    e.rig.(rignames{v}) = [];
end
% Save grid steps and real coordinates
e.q = q; e.Iota = Iota; e.J = J; e.K = K; e.r = r; e.t = t;

% Save configuation information
e.nrt = des.nrt; e.mrtp = des.mrtp; e.N = N; e.geom_name = geom_name;

% Save probe information and calibration map
e.probe_name = probe_name; e.probe = probe.(probe_name);
e.calib_name = calib_name; e.c = c; 


%% Open traverse monitoring matlab session

% Traverse name
trav_suffix = [num2str(ni) 'x' num2str(nj)];

% Determine new filename - no overwriting
n = 1;
file_prefix = [geom_name '_' trav_suffix '_NRT' strrep(num2str(round(des.nrt*100)),'.','') '_' date];
file_prefix = strrep(file_prefix,'.','p');
while exist([file_prefix '_' num2str(n) '.mat'],'file') ~= 0
    n = n+1;
end
filename = [file_prefix '_' num2str(n)];

% Save initial traverse data
save([filename '.mat'],'e')

% Open new matlab instance with monitoring function
save('current_traverse.mat')
system(['matlab -nosplash -nodesktop -r exp_traverse_monitor(''' filename ''') &'])
% exp_traverse_monitor(''' filename ''')

%% Start rig, set flow coefficient and log pressure and temperature

% Start rig and set the shaft speed for the first time
uiwait(msgbox('Start Rig'));

rpm = 0;
while abs(rpm - rpm_target) / rpm_target > 0.05
    rpm = exp_shaft_freq(s.pxie,1,0) * 30 / pi;
    disp(['RPM = ' num2str(rpm)]);
end

% Log starting temperature and atmospheric pressure
  prompt = {'Record rig temperature /K', 'Record atmospheric pressure /kPa'};
  dlgtitle = 'Start readings';
  T1 = inputdlg(prompt(1));
  P1 = inputdlg(prompt(2));


%% Execute traverse

% % Move probe to starting position in reverse axis order
% i = I(1); j = J(1);
% q_temp = squeeze(q(i,j,:));
% exp_stepper_coord(s.park,[0 ; 0 ; q_temp],park.limit);

% Move probe to starting position in reverse axis order
j = J(1); k = K(1);
if ni == 1
    q_temp = squeeze(q(j,k,:));
else
    i = 1;
    q_temp = squeeze(q(j,k,i,:));
end
% exp_stepper_coord(s.park,[0 0 q_temp(3)],limit);
% exp_stepper_coord(s.park,[0 q_temp(2) q_temp(3)],limit);
% exp_stepper_coord(s.park,q_temp,limit);
exp_stepper_coord(s.park,[0 ; 0 ; q_temp(1:2)],park.limit);


% Loop over all points in grid, moving in either tangential or radial
% direction, record and process traverse data
% for n = 1:ni*nj
%     
%     % Get current coordinates
%     i = I(n); j = J(n); q_old = q_temp;
%     q_temp = squeeze(q(i,j,:));
%     
    
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
    disp(['Measuring Point J = ' num2str(i) ' K = ' num2str(j)]);
    
    % Move to new location
    exp_stepper_coord(s.park,[0 ; 0 ; q_temp(1:2)],park.limit);
    
    % Check if full traverse with instrument logging is specified
    if strcmp(run_type,'full') == 1

%         Get and log the shaft speed
        tic;
        rpm = exp_shaft_freq(s.pxie,1,0) * 30 / pi;
        disp(['RPM = ' num2str(rpm)]);
        e.rpm(i,j,:) = rpm;
        
%         % Log tc08 temperatures
%         M = exp_tc08_read(s.tc08);
%         e.temp(i,j,1) = M(1);
%         e.temp(i,j,2) = M(2);
%         disp(['Inlet temp = ' num2str(e.temp(i,j,1))]);
%         disp(['Motor temp = ' num2str(e.temp(i,j,2))]);
%         
        
        % Pause for remaining probe settling
        settled_time = toc;
        pause(max(probe.(probe_name).time - settled_time,0));

        % Log DSA pressures
        e.P(i,j,:) = exp_dsa_read(s.dsa);

        % Log time
        e.time(i,j,:) = clock;

        % Save data to disk
        save([filename '.mat'],'e')
    end    
end

% Return motors to zeros and re-home
exp_stepper_coord(s.park,[0 0 0 0],park.limit);

% Turn off motors and quit the script if grid test is complete
if strcmp(run_type,'off') == 1
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


%% Record end rig temperature and atmospheric pressure

%Log end temperature and atmospheric pressure
prompt = {'Record rig temperature /K', 'Record atmospheric pressure /kPa'};
dlgtitle = 'End readings';
T2 = inputdlg(prompt(1));
P2 = inputdlg(prompt(2));

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

% Disconnect motors
for n = 1:park.nchan
    exp_stepper_command(s.park,n,'OFF');
    exp_stepper_command(s.park,n,'E0');
end
fclose(s.park); delete(s.park); s = rmfield(s,'park');

% Close pxie connection
delete(s.pxie);

% Close thermocouple sensor
% calllib('usbtc08','usb_tc08_close_unit',s.tc08.h);
% s = rmfield(s,'tc08');

% Clear up connection structure
clear('s'); instrreset; delete(instrfindall);


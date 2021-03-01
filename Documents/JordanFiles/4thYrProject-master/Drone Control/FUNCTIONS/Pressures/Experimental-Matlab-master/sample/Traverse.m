% Run traverses on Gibbons rig

% Stepper motors
% Channel 1 is radial direction
% Channel 2 is tangential direction
% Channel 3 is probe yaw or axial direction

% DSA
% Channel 1 is probe centre
% Channel 2 is probe upper left         ] 
% Channel 3 is probe upper right        ] (When looking at probe head with casing above)
% Channel 4 is probe lower right        ]
% Channel 5 is probe lower left         ]
% Channel 6 is bell total (1)
% Channel 7 is bell casing static (1)
% Channel 8 is IGV inlet casing static (2)
% Channel 9 is rotor inlet casing static (3)
% Channel 10 is stator inlet mid total in sector (4)
% Channel 11 is stator inlet casing static in sector(4)
% Channel 12 is stator inlet casing static (4)
% Channel 13 is stator exit casing static (5)
% Channel 14 is downstream exit casing static (6)
% Channel 15 & 16 is open to atmosphere

% Thermocouples
% Channel 1 is IGV inlet SMURF side (2)
% Channel 2 is IGV inlet drive side (2)
% Channel 3 is IGV inlet crane side (2)
% Channel 4 is stator inlet in traverse sector (4)
% Channel 5 is in rotor casing on drive side (-)

% Clear workspace
close all; clear all; clc; instrreset; delete(instrfindall);

% Load Gibbons Data
Gibbons_Exp

% Record channels for useful data
N.Po_1 = 6;
N.P_1_cas = 7;
N.P_2_cas = 8;
N.P_3_cas = 9;
N.Po_4_sec = 10;
N.P_4_sec = 11;
N.P_4_cas = 12;
N.P_5_cas = 13;
N.P_6_cas = 14;
N.Pa = 16;
N.P_ref = 12;

N.To_2_mid = 1:3;
N.To_4_sec = 4;
% N.T_cas = 5;
% N.T_hub = 6;

N.V_shaft = 1; % AI 0

N_rig.P = N.P_2_cas; N_rig.T = N.To_2_mid(2); N_rig.Po = N.Po_1; N_rig.Pa = N.Pa; N_rig.V = 1;

% Select probe
N.FHP = 1:5; probe_name = 'FHP_G'; calib_name = 'FHP_G_61_fine_18-Nov-2015_Bingo_1.mat'; 
% N.ToP = 5; probe_name = 'ToP_A'; nk = round(round(nk/2)/2)*2+1; nj = round(nj/2);
% N.V_hw = 2; probe_name = 'HW_B'; calib_name = 'HW_B_velocity_06-Nov-2015_1.mat'; V_thresh = 4;
% N.FHP = 1:5; probe_name = 'FHP_HL'; calib_name = 'FHP_HL_61_fine_17-Nov-2015_Bingo_1.mat'; 

% Select traverse configuration

% Midheight at IGV exit for uniformity of calibration area
% trav_suffix = '11x1'; st = 3; bl = 2; dir_type = 'tang'; Iota = 10; 
% nk = 11; t_min = 2; t_max = 6.8; t_clust = 'linear';
% nj = 1; r_min = 0.5 * (r_hub + r_cas); r_max = 0.5 * (r_hub + r_cas); r_fac = [];

% Midheight at IGV exit for yaw calibrations
% trav_suffix = '1x1x80'; st = 3; bl = 2; dir_type = 'tang'; 
% Iota = [linspace(-20,40,41) linspace(39,-19,39)]; 
% nk = 80; t_min = 2+0.01; t_max = 2-0.01; t_clust = 'linear';
% nj = 1; r_min = 0.5 * (r_hub + r_cas); r_max = 0.5 * (r_hub + r_cas); r_fac = [];

% Single passage at IGV exit
% trav_suffix = '19x25'; st = 3; bl = 2; dir_type = 'radial'; Iota = 10; ni = 1;
% nk = 19; t_min = 3.9; t_max = 8.7; t_clust = 6.3; t_fac = 0.3;
% nj = 25; r_min = r_hub + 0.5e-3 + probe.(probe_name).D/2; 
% r_max = r_cas - 0.5e-3 - probe.(probe_name).D/2; r_fac = 0.15;

% Fine single passage at stator inlet
% trav_suffix = '15x27'; st = 4; bl = 2; dir_type = 'radial'; Iota = 50; ni = 1;
% nk = 15; t_clust = 'linear'; 
% % t_min = -2.5; t_max = 2.3; t_pitch = 4.8; % 75 Blades
% % t_min = -2; t_max = 3.14; t_pitch = 5.14; % 70 Blades
% % t_min = -2.5; t_max = 3.04; t_pitch = 5.54; % 65 Blades
% t_min = -2.9; t_max = 5.5; t_pitch = 8.4; % Largest span
% nj = 27; r_min = r_hub + 0.7e-3 + probe.(probe_name).D/2; 
% r_max = r_cas - 0.4e-3 - probe.(probe_name).D/2; r_fac = 0.15;

% Coarse single passage at stator inlet
% trav_suffix = '7x15'; st = 4; bl = 2; dir_type = 'radial'; Iota = 50; ni = 1;
% nk = 7; t_clust = 'linear'; 
% % t_min = -2.5; t_max = 2.3; t_pitch = 4.8; % 75 Blades
% % t_min = -2; t_max = 3.14; t_pitch = 5.14; % 70 Blades
% % t_min = -2.5; t_max = 3.04; t_pitch = 5.54; % 65 Blades
% t_min = -2.9; t_max = 5.5; t_pitch = 8.4; % Largest span
% nj = 15; r_min = r_hub + 0.7e-3 + probe.(probe_name).D/2; 
% r_max = r_cas - 0.4e-3 - probe.(probe_name).D/2; r_fac = 0.15;

% Fine double passage at stator exit
trav_suffix = '54x27'; st = 5; bl = [1 2]; dir_type = 'tang'; Iota = 13; ni = 1;
nk = 27; t_fac = [0.22 0.15]; t_clust = 'geom';
nj = 27;  r_min = r_hub + 0.3e-3 + probe.(probe_name).D/2; 
r_max = r_cas - 0.65e-3 - probe.(probe_name).D/2; r_fac = 0.15;

% Coarse single passage at stator exit
% trav_suffix = '15x15'; st = 5; bl = 2; dir_type = 'tang'; Iota = 13; ni = 1;
% nk = 15; t_fac = [0.22 0.15]; t_clust = 'geom';
% nj = 15;  r_min = r_hub + 0.3e-3 + probe.(probe_name).D/2; 
% r_max = r_cas - 0.65e-3 - probe.(probe_name).D/2; r_fac = 0.15;

% Fine single passage at stator exit
% trav_suffix = '27x27'; st = 5; bl = 2; dir_type = 'tang'; Iota = 13; ni = 1;
% nk = 27; t_fac = [0.22 0.15]; t_clust = 'geom';
% nj = 27;  r_min = r_hub + 0.3e-3 + probe.(probe_name).D/2; 
% r_max = r_cas - 0.65e-3 - probe.(probe_name).D/2; r_fac = 0.15;

% Fine single passage midpsan traverse stator exit 
% trav_suffix = '107x1'; st = 5; bl = 2; dir_type = 'tang'; Iota = 13; ni = 1;
% nk = 107; t_fac = [0.22 0.15]; t_clust = 'geom';
% nj = 1;  r_min = r_mid; r_max = r_mid; r_fac = [];


% Hot wire stator exit
% trav_suffix = '17x17'; st = 5; bl = 2; dir_type = 'tang';
% ni = 6; Iota_min = -20; Iota_max = 50;
% nk = 17; t_fac = [0.22 0.15]; t_clust = 'geom';
% nj = 17;  r_min = r_hub + 1e-3; 
% r_max = r_cas - 1e-3; r_fac = 0.15;

% Hot wire stator inlet
% trav_suffix = '3x27x6'; st = 4; bl = 2; dir_type = 'radial';
% ni = 6; Iota_min = 20; Iota_max = 90;
% nk = 3; t_clust = 'linear';
% % t_min = -2.5; t_max = 2.3; t_pitch = 4.8; % 75 Blades
% t_min = 1; t_max = 5.8; t_pitch = 4.8; % 75 Blades - Ivor's Ring
% % t_min = -2; t_max = 3.14; t_pitch = 5.14; % 70 Blades
% % t_min = -2.5; t_max = 3.04; t_pitch = 5.54; % 65 Blades
% t_max = t_max - t_pitch/nk;
% nj = 27;  r_min = r_hub + 1.3e-3; 
% r_max = r_cas - 1.1e-3; r_fac = 0.15;



% Limits for wake searching
% search.t_min = 3; search.t_max = 7.8; search.t_pitch = 4.8; % 75 Blades - Ivor's Ring
search.t_min = 1; search.t_max = 5.8; search.t_pitch = 4.8; % 75 Blades
% search.t_min = 1.4; search.t_max = 6.77; search.t_pitch = 5.37; % 67 Blades - Bronwyn
% search.t_min = 3; search.t_max = 8.14; search.t_pitch = 5.14; % 70 Blades
% search.t_min = 0; search.t_max = 5.54; search.t_pitch = 5.54; % 65 Blades
search.r_min = r_hub + 8e-3; search.r_max = r_cas - 8e-3;
search.nk = 12; search.nj = 3; search.Iota = 13;
search.override = 1;

% Set tolerances in operating point to 0.25%
err.phi = 0.0025;
err.nrt = 0.0025;
err.mrtp = 0.0025; 

% Set desired operating points
des.mrtp = 0.4746; des.nrt = 0.6823;                               % Design Point
% des.phi = 0.5891;                                                  % Design Point
% des.mrtp = [0.4746 0.385 0.36]; des.nrt = ones(1,3) * 0.6823;      % Winglet hot wires
% des.mrtp = [0.4746 0.385]; des.nrt = ones(1,2) * 0.6823;           % Winglet hot wires
% des.mrtp = [0.4746 0.44 0.4 0.36]; des.nrt = ones(1,4) * 0.6823;   % Bronwyn stator test

% Set geometry
% geom_name = 'SLA_S100_L100_N75';
geom_name = 'SLA_150B_SS'; 

% geom_name = 'SLA_S30_L30_N70'; 
% geom_name = 'SLA_S100_L100_N70'; 
% geom_name = 'SLA_S170_L170_N70'; 

% geom_name = 'SLA_S0_L0_N65';
% geom_name = 'SLA_S100_L100_N65'; 
% geom_name = 'SLA_S200_L200_N65';

%% Determine traverse type and load calibration map

if strcmp(probe_name(1:3),'FHP') == 1
    
    % Five hole probe traverse
    trav_type = 'FH';
    load([directory calib_name]);
    dsa.nchan = 16; tc08.nchan = 4; pxie.nchan = 1;
    pxie.chan = N.V_shaft; pxie.time = 0.2;
    
elseif strcmp(probe_name(1:2),'HW') == 1
    
    % Hot wire traverse
    trav_type = 'HW';
    load([directory calib_name]);
    dsa.nchan; tc08.nchan = 3; pxie.nchan = 2; 
    pxie.chan = [N.V_shaft N.V_hw]; pxie.time = probe.(probe_name).time;
    
elseif strcmp(probe_name(1:3),'ToP') == 1
    
    % Total temperature traverse
    trav_type = 'To';
    dsa.nchan = 16; tc08.nchan = 6;
    
end

%% Initialise all instruments

% Connect and power up stepper motors
s.park = exp_stepper_initialise(park,0);

% Initialise and zero DSA
s.dsa = exp_dsa_initialise(dsa);

% Initialise thermocouple logger
s.tc08 = exp_tc08_initialise(tc08);

% Initialise once per ref logging
s.pxie = exp_pxie_initialise(pxie,'DC');

% Point to low speed lab PTH sensor
s.ptu = ptu;

% Home all axes simultaneously
exp_stepper_home(s.park,park.home,1)

% Narrow limit range based on traverse location
park.limit = park.limit.(['T' num2str(st)]);

% Alter limit range based upon probe offset
park.limit(1,:) = park.limit(1,:) + park.d(1) * (10.82e-3 + probe.(probe_name).r_off); 
% park.limit(3,:) = park.limit(3,:) + park.d(3) * (0 + probe.(probe_name).Iota_off); 
park.limit = round(park.limit);

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
for a = 1:length(phi_loop)

if isfield(des,'phi') == 0
    des_current.mrtp = des.mrtp(a);
    des_current.nrt = des.nrt(a);
else
    des_current.phi = des.phi(a);
end
    
% Assemble r coordinates for traverse clustering on endwalls if necessary
if isempty(r_fac) == 0
    r = hyperbolic_bunch(nj,r_fac/nj,r_fac/nj) * (r_max - r_min) + r_min;
else
    r = linspace(r_min,r_max,nj);
end
r = repmat(reshape(r,[nj 1]),[1 nk]);

% Determine theta clustering as a function of r
if exist('t_fac','var') ~= 0 && length(t_fac) > 1
    t_facs = ones(nj,1);
    for j = 1:nj
        t_facs(j) = 2*abs((j-round(nj/2))/nj) * (t_fac(1) - t_fac(2)) + t_fac(2);
    end
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
    
    % Check if wakes have already been found
    A = dir([directory geom_name '*Ph0' strrep(num2str(phi_loop(a)*100),'.','')...
        '*T' num2str(st) '_*.mat']);
    if isempty(A) == 0
        i = 0; time_temp = 0;
        for n = 1:length(A)
            if A(n).datenum > time_temp
                time_temp = A(n).datenum;
                i = n;
            end
        end
        load([directory A(i).name]);
    else
        e.tr_wake = [];
    end
    
    % Determine where to get wake location information from
    if isfield(e,'tr_wake') == 1 && isempty(e.tr_wake) == 0 && search.override == 0
        % Read existing wake location
        tr_wake = e.tr_wake;
        clear('e');
    elseif exist('tr_wake','var') == 0 || search.override == 1
        % Search for and non-uniformly cluster on wakes
        tr_wake = exp_traverse_wakes(s,search,park,probe.(probe_name),N_rig,des_current,err,r_mid);
    end
    e.tr_wake = tr_wake; t_pitch = search.t_pitch;
    
    % Calculate quadraticly varying traverse grid from searched wakes
    p = polyfit(tr_wake(:,2),tr_wake(:,1),2); clear tr_wake;
    t = zeros(nj,nk);
    for j = 1:nj
        t_mid = polyval(p,r(j,1)); 
        t_max = t_mid + t_pitch/2;
        t_min = t_mid - t_pitch/2;
        t(j,:) = unique([hyperbolic_bunch(round((nk+1)/2),t_facs(j),1) * (t_mid - t_min) + t_min ...
            hyperbolic_bunch(round((nk+1)/2),1,t_facs(j)) * (t_max - t_mid) + t_mid]);
    end
    
elseif strcmp(t_clust,'linear') == 1
    
    % Uniform distribution
    t = repmat(reshape(linspace(t_min,t_max,nk),[1 nk]),[nj 1]);
    
end

% Repeat t and r coordinates for multiple passage traverses
t_1 = t; r_1 = r;
t_2 = zeros(nj,nk*length(bl)-(length(bl)-1));
r_2 = zeros(nj,nk*length(bl)-(length(bl)-1));
for n = length(bl):-1:1
    t_2(:,n*nk-nk+1:n*nk) = t_1 - (bl(n)-2) * t_pitch + n*1e-6;
    r_2(:,n*nk-nk+1:n*nk) = r_1;
end

t = []; r = [];
for j = 1:size(r_2,1)
    [t_temp,i] = unique(t_2(j,:));
    t = [t ; t_temp]; r = [r ; r_2(j,i)];
end
   
% Define angle matrix
if ni == 1
    Iota = Iota .* ones(nj,nk*length(bl));
    dim_coord = 3;
else
    Iota = repmat(reshape(linspace(Iota_min,Iota_max,ni),[1 1 ni]),[nj nk 1]);
    r = repmat(r,[1 1 ni]); t = repmat(t,[1 1 ni]);
    dim_coord = 4;
end

% Smooth tangenetial grid based upon limits
t_lim = (park.limit(2,1)+200 - park.c(2)) / park.d(2);
for j = 1:nj
    if sum(t(j,:) < t_lim) > 0
        dt = t_lim - t(j,1);
        t(j,[1 end]) = t(j,[1 end]) + dt;
        t(j,[2 end-1]) = t(j,[2 end-1]) + 0.5*dt;
        t(j,[3 end-2]) = t(j,[3 end-2]) + 0.2*dt;
    end
%     elseif sum(t(j,:) < t_lim) > 1
%         disp('Two Points out of Range');
%     end
end

% Calculate traverse grid in steps
q = round(cat(dim_coord,(r + probe.(probe_name).r_off) * park.d(1) + park.c(1),...
    t*park.d(2) + park.c(2),...
    (Iota + probe.(probe_name).Iota_off) * park.d(3) + park.c(3)));

% Create index vectors for either circumferential or radial stepping
if strcmp(dir_type,'radial') == 1
    J = repmat(reshape(1:nj,[nj 1]),[1 nk*length(bl)]); 
    K = repmat(reshape(1:nk*length(bl),[1 nk*length(bl)]),[nj 1]); 
    J = reshape(J,[],1); K = reshape(K,[],1);
elseif strcmp(dir_type,'tang') == 1
    J = repmat(reshape(nj:-1:1,[nj 1]),[1 nk*length(bl)]); 
    K = repmat(reshape(1:nk*length(bl),[1 nk*length(bl)]),[nj 1]); 
    J = reshape(permute(J,[2 1 3]),[],1); K = reshape(permute(K,[2 1 3]),[],1);
end

% Display minimum and maximum coordinates in each dimension
if ni == 1
    disp(reshape(max(max(q,[],1),[],2),[1 3]))
    disp(reshape(min(min(q,[],1),[],2),[1 3]))
else
    disp(reshape(max(max(max(q,[],1),[],2),[],3),[1 3]))
    disp(reshape(min(min(min(q,[],1),[],2),[],3),[1 3]))
end

%% Initialise data structure

% Initialise raw pressure matrices
if strcmp(trav_type,'FH') == 1
    e.P = nan(nj,nk*length(bl),dsa.nchan);
    e.time = nan(nj,nk*length(bl),6);
elseif strcmp(trav_type,'HW') == 1
    e.V_av = nan(nj,nk*length(bl),ni);
    e.time = nan(nj,nk*length(bl),ni,6);
end

% Initialise raw rig operating point information
rignames = [{'P' 'T' 'Pa' 'Ha' 'R' 'cp' 'ga' 'w' 'time'} fieldnames(des)'];
for v = 1:length(rignames)
    e.rig.(rignames{v}) = [];
end

% Save grid steps and real coordinates
e.q = q; e.r = r; e.t = t; e.Iota = Iota;
e.J = J; e.K = K;

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

% Determine new filename - no overwriting
n = 1;
file_prefix = [geom_name '_' trav_suffix '_Ph0' strrep(num2str(phi_loop(a)*100),'.','') '_' date...
    '_T' num2str(st) '_' trav_type];
while exist([directory file_prefix '_' num2str(n) '.mat'],'file') ~= 0
    n = n+1;
end
filename = [directory file_prefix '_' num2str(n)];

% Save initial traverse data
save([filename '.mat'],'e')

% Open new matlab instance with monitoring function
save('current_traverse.mat')
% if strcmp(trav_type,'HW') == 0
    system(['matlab -nosplash -nodesktop -r exp_traverse_monitor(''' filename ''') &'])
% end

% Save hot wire raw voltage files
if strcmp(trav_type,'HW') == 1
    nt = pxie.rate*pxie.time+1;
    V_hw = nan(nj, nk*length(bl), ni, nt, 'single');
    V_shaft = true(nj, nk*length(bl), ni, nt);
end

%% Start rig and set flow coefficient

% Record zero flow voltage
if strcmp(trav_type,'HW') == 1
    uiwait(msgbox('Rig off zero flow voltage'));
    [V_temp,t_temp] = exp_pxie_read(s.pxie);
    disp(mean(V_temp));
    e.V0_start = single(V_temp(:,2));
end

uiwait(msgbox('Start Rig & Set Traverse Support'));

% Set rig operating point for the first time
[e_rig,on_point] = exp_rig_set(s,N_rig,des_current,err,r_mid);

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
for n = 1:nj*nk*length(bl)
    
    % Get current coordinates
    j = J(n); k = K(n);
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
                exp_stepper_coord(s.park,[0 q_temp(2)-200 q_temp(3)],park.limit);
            end
            
            exp_stepper_coord(s.park,[q_temp(1) q_temp(2)-200 q_temp(3)],park.limit);
        end
        
        % Position HW at mid radius first
        if strcmp(trav_type,'HW') == 1
            exp_stepper_move(s.park,1,0)
            exp_stepper_position(s.park,1,0)
        end
        
        % Move to new tangential coordinates
        exp_stepper_move(s.park,2,q_temp(2))
        
        % Query tangential position
        exp_stepper_position(s.park,2,q_temp(2))
        
        % Reposition HW at intended radius
        if strcmp(trav_type,'HW') == 1
            exp_stepper_move(s.park,1,q_temp(1))
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
        exp_stepper_move(s.park,1,q_temp(1))
        
        % Query radial position
        exp_stepper_position(s.park,1,q_temp(1))
    end
    
    % Rotate the probe if desired angle has changed
    if ni == 1 && n > 1 && q(j,k,3) ~= q(J(n-1),K(n-1),3)
        exp_stepper_coord(s.park,q(j,k,:),park.limit);
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
            exp_stepper_move(s.park,3,q(j,k,i,3))
            exp_stepper_position(s.park,3,q(j,k,i,3))
        end
        
        % Check if rig is on point while probe is settling
        tic; 
        [e.rig,on_point,e_temp] = exp_rig_monitor(s,e.rig,N_rig,des_current,err,r_mid); 

        % Reset rig operating point if not on
        if on_point == 0 && i == 1
            exp_rig_set(s,N_rig,des_current,err,r_mid);
            [e.rig,on_point,e_temp] = exp_rig_monitor(s,e.rig,N_rig,des_current,err,r_mid); 
        end

        % Pause for remaining probe settling
        settled_time = toc;
        pause(max(probe.(probe_name).time - settled_time,0));

        % Log current instrument
        if strcmp(trav_type,'FH') == 1
            % Log DSA pressures
            e.P(j,k,:) = exp_dsa_read(s.dsa);
        elseif strcmp(trav_type,'HW') == 1
            % Record hotwire voltages from rig monitoring data logging
            V_hw(j,k,i,:) = e_temp.V(:,2);
            e.V_av(j,k,i) = mean(e_temp.V(:,2));
            V_shaft(j,k,i,e_temp.V(:,1) < V_thresh) = 0;
            e.dtime = e_temp.dtime;
            disp(['Vav = ' num2str(e.V_av(j,k,i))]);
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
    
    % Save row of hotwire voltages
    if strcmp(trav_type,'HW') == 1
        if k == nk*length(bl) && strcmp(dir_type,'tang') == 1
            e_temp.V_hw = V_hw(j,1:nk*length(bl),1:ni,1:nt); 
            e_temp.V_shaft = V_shaft(j,1:nk*length(bl),1:ni,1:nt);
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

% Save all data
if strcmp(trav_type,'HW') == 1
    e.V_shaft = V_shaft; clear V_shaft;
    e.V_hw = V_hw; clear V_hw;
    
    uiwait(msgbox('Rig off zero flow voltage'));
    [V_temp,t_temp] = exp_pxie_read(s.pxie);
    e.V0_end = single(V_temp(:,2));
    disp(mean(V_temp));
end
save([filename '.mat'],'e')
clear e;

% Finish looping over flow coefficients
end

%% Stop rig and close all instruments down

uiwait(msgbox('Stop Rig'));

% Disconnect dsa
fclose(s.dsa); delete(s.dsa); s = rmfield(s,'dsa');

% Close thermocouple sensor
calllib('usbtc08','usb_tc08_close_unit',s.tc08.h);
s = rmfield(s,'tc08');

% Disconnect motors
for n = 1:3
    exp_stepper_command(s.park,n,'OFF');
    exp_stepper_command(s.park,n,'E0');
end
fclose(s.park); delete(s.park); s = rmfield(s,'park');

% Close pxie connection
delete(s.pxie);

% Clear up connection structure
clear('s'); instrreset; delete(instrfindall);


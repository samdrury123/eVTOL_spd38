% Calibrate 5 hole probe in bingo machine or nodding donkey

% Stepper motors
% Channel 1 is blank
% Channel 2 is blank
% Channel 3 is yaw angle
% Channel 4 is pitch angle

% DSA
% Channel 1 is nozzle stagnation
% Channel 2 is jet static
% Channel 3 is probe centre
% Channel 4 is probe up         ] 
% Channel 5 is probe down       ] (When looking upstream with casing above)
% Channel 6 is probe right      ]
% Channel 7 is probe left       ]

% TC08
% Channel 1 is atmospheric temperature
% Channel 2 is probe temperature

% Clear workspace
close all; clear; clc; instrreset; delete(instrfindall);

%Load eVTOL data
eVTOL; dc = 0; settle_fac = 2; V_des = 30;

% Load C106 data
% C106; dc = -18000; dsa = dsa{2}; settle_fac = 7; V_des = 45;

% Constants for bingo machine & 5 hole probe
park.port = 'COM1'; park.rate = 9600; park.term = 'CR'; 
Iota_d = -200; Iota_c = 0 + dc; Tau_d = 59.2; Tau_c = 0; machine_name = 'Bingo';
limit = [0 72000 ; -1480 1480 ; -36000 36000 ; -1480 1480]; park.nchan = 4;
park.current = [0.6 1.5 0.6 1.5]; park.resistance = [4.2 2 4.2 2]; park.inductance = [2.5 1 2.5 1];
park.V = [10 2 10 2]; park.A = [100 20 100 20];

% Constants for nodding donkey & 5 hole probe
% Iota_d = 1000; Iota_c = 90000; Tau_d = -4068; Tau_c = 0;  machine_name = 'Donkey';
% limit = [0 360000 ; -100000 100000];
    
% Select probe
probe_name = 'FHP_M';

% Custom five hole
% probe_name = 'FHP_FERG';
% probe.FHP_FERG.Iota_off = 0; probe.FHP_FERG.time = 0.5; probe.FHP_FERG.arrange = 'plus';
% N.P_cn = 1; N.P_up = 2; N.P_dn = 4; N.P_lf = 3; N.P_rt = 5; probe.FHP_FERG.N = N;

% Fine calibration grid
Iota_min = -36; Iota_max = 36; Iota_n = 37;
Tau_min = -20; Tau_max = 20; Tau_n = 21;
calib_name = [num2str(V_des) '_fine'];
rev_iota = 0; 

% Test coarse calibration grid
% Iota_min = -10; Iota_max = 10; Iota_n = 11;
% Tau_min = -10; Tau_max = 10; Tau_n = 11;
% calib_name = [num2str(V_des) '_testcoarse'];
% rev_iota = 0; 

% Coarse calibration grid to test probe
% Iota_min = -32; Iota_max = 36;
% Tau_min = -20; Tau_max = 20;
% Iota_n = 5; Tau_n = 5;
%  C_des = [0 0 0]; calib_name = [num2str(V_des) '_coarse']; 
% rev_iota = 0;

% Yaw to check settling time
% Iota_min = -36; Iota_max = 36; Iota_n = 19;
% Tau_min = 0; Tau_max = 0; Tau_n = 1; 
% C_des = [0 0 0]; calib_name = [num2str(V_des) '_yaw']; rev_iota = 1;

% Pressure channels
N = probe.(probe_name).N; varnames = fieldnames(N); 
for v = 1:length(varnames); N.(varnames{v}) = N.(varnames{v}) + 2; end;
N.Po = 1; N.P = 2; 


%% Determine calibration type
if strcmp(probe_name(1:3),'FHP') == 1
    calib_fhp = 1;
    calib_top = 0;
    p_chan = 7; T_chan = 0;
elseif strcmp(probe_name(1:3),'ToP') == 1
    calib_fhp = 0;
    calib_top = 1;
    p_chan = 2; T_chan = 2; tc08.nchan = 2;
    settle_fac = 1;
end


%% Set atmospheric data

[Pa,Ta,Ha] = getCurrentPTH(ptu.loc);
ro = Pa ./ (air.R * Ta); 
% 

%% Take readings from calibration tunnel

% Define angle matrices
Iota = repmat(reshape(linspace(Iota_min,Iota_max,Iota_n),[Iota_n 1]),[1 Tau_n]);
Tau = repmat(reshape(linspace(Tau_min,Tau_max,Tau_n),[1 Tau_n]),[Iota_n 1]);

% Drive backwards to check repeatability
if rev_iota == 1
    Iota = [Iota ; flip(Iota,1)];
    Tau = [Tau ; flip(Tau,1)];
    Iota_n = size(Iota,1); Tau_n = size(Tau,2);
end

% Calculate calibration grid in steps
q = round(cat(3,+Iota*Iota_d + Iota_c,Tau * Tau_d + Tau_c));

% Initialise pressure and temperature arrays
P_raw = zeros([Iota_n Tau_n length(V_des) p_chan]);
if calib_top == 1
    T_raw = zeros([Iota_n Tau_n length(V_des) T_chan]);
end

% Initialise stepper motors and move to central position
s.park = exp_stepper_initialise(park,0);
exp_stepper_command(s.park,2,'V2');
exp_stepper_command(s.park,2,'A5');

exp_stepper_coord(s.park,[0 0 squeeze(q(round(Iota_n/2),round(Tau_n/2),:))'],limit)
disp('%% Stepper connected, zeroed and moved to central position %%')

% Initialise DSA with tunnel off
s.dsa = exp_dsa_initialise(dsa);
disp('%% DSA zeroed and initialised %%')

% Initialise thermocouple logger
if calib_top == 1
    s.tc08 = exp_tc08_initialise(tc08);
    disp('%% TC08 initialised %%')
end

% Open data logging figure window
h.track = figure('position',[1 29 1920 985]); hold on; grid on; box on;
xlabel('Yaw Coefficient'); ylabel('Pitch Coefficient'); C = lines(length(V_des));

% Start tunnel
uiwait(msgbox('Turn Tunnel On'));

% Loop over all desired tunnel speeds recording pressures and moving probe
for n = 1:length(V_des)

    % Set tunnel speed
%     fprintf(t, 'SET AVG 1\n');

    % Bring up jet velocity trace
    h.speed = figure(); hold on; grid on; xlabel('Time / s'); ylabel('Jet Velocity / ms^{-1}');
    time = clock;
    time_plot = 0;
    
    % Initialise looping variables for axes and termination
    V = 0; V_min = 1e3; V_max = -1e3; V_flag = 0; V_log = [];
    disp(['%% Set tunnel speed to ' num2str(V_des(n)) ' meters per second'])
    
    g = 0; b = uicontrol('style','push','string','V Set','callback','g=g+1');

    while V_flag == 0 && g == 0
        % Read DSA and time
        P = exp_dsa_read(s.dsa);
        time_new = clock;
        time_plot_new = (time_new(end) - time(end)) + (time_new(end-1) - time(end-1)) * 60 + ...
            (time_new(end-2) - time(end-2)) * 3600;

        % Calculate jet velocity
        Po = P(N.Po); P_s = P(N.P);
        V_new = ( 2 * (Po - P_s) / ro ) .^ 0.5;

        % Plot jet velocity
        figure(h.speed)
        plot([time_plot ; time_plot_new],[V ; V_new],'k-');
        time_plot = time_plot_new; V = V_new; V_log = [V_log ; V_new];
        V_min = min(V_min,min(V_new)); V_max = max(V_max,max(V_new)); 
        axis([time_plot_new - 300 time_plot_new V_min - 1 V_max + 1]);

        % Compare jet velocity to target
        if length(V_log) > 15
            if sum(abs(V_log(end-6:end) - V_des(n)) / V_des(n) < 0.005) > 4
                V_flag = 1;
            end
        end
        
        % Check if velocity set button has been pressed
        drawnow;
    end
    close();

    % Loop over all calibration points and record pressures on all channels
    figure(h.track);
    for j = 1:Tau_n
        for i = 1:Iota_n

            % Move stepper motors to coordinates
            exp_stepper_coord(s.park,[0 0 q(i,j,1) q(i,j,2)],limit)

            % Pause double at starts of cycles
            if i == 1
                pause(probe.(probe_name).time * settle_fac);
            end
            
            % Pause for stepper ringing and probe settling constant
            pause(probe.(probe_name).time * settle_fac);

            % Log DSA
            P = exp_dsa_read(s.dsa);
            P_raw(i,j,n,:) = P(1:p_chan);
            
            % Log thermocouples
            if calib_top == 1
                T = exp_tc08_read(s.tc08,tc08,T_chan);
                T_raw(i,j,n,:) = T(1:p_chan);
            end
            
            % Display point
            fprintf(['%% Logged point:\t\tYaw = ' num2str(Iota(i,j)) ...
                ' \t\tPitch = ' num2str(Tau(i,j)) ' %%\n'])
            
            % Plot point
            c_temp = exp_fhp_calibrate(Iota(i,j),Tau(i,j),P_raw(i,j,n,:),...
                probe.(probe_name),Ta,'hodson',0);
            plot(c_temp.C_Alpha,c_temp.C_Beta,'o','color',C(n,:),'markerfacecolor',C(n,:));
        end
    end
    
    % Plot spiders web for last tunnel speed
%     if calib_fhp == 1
%         c_temp = exp_fhp_calibrate(Iota,Tau,P_raw(:,:,n,:),probe.(probe_name),Ta,'hodson',1);
%     end
    
    % Plot temperature map for last tunnel speed
    if calib_top == 1
        figure(); 
        subplot(1,2,1); surf(Iota,Tau,T_raw(:,:,1,1)); colorbar;
        subplot(1,2,2); surf(Iota,Tau,T_raw(:,:,1,2)); colorbar;
        figure();
        surf(Iota,Tau,diff(T_raw,1,4)); colorbar;
    end
end

% Return motors to zeros 
exp_stepper_coord(s.park,[0 0 0 0],limit);

% Disconnect motors
for n = 1:park.nchan
    exp_stepper_command(s.park,n,'OFF');
    exp_stepper_command(s.park,n,'E0');
end
fclose(s.park); s = rmfield(s,'park');
disp('%% Stepper zeroed and disconnected %%')

% Power down tunnel
uiwait(msgbox('Turn Tunnel Off'))

% Disconnect dsa
fclose(s.dsa); s = rmfield(s,'dsa');
disp('%% DSA logged for zero and disconnected %%')

% Disconnect tc08
if calib_top == 1
    calllib('usbtc08','usb_tc08_close_unit',s.tc08);
    s = rmfield(s,'tc08');
    disp('%% TC08 disconnected %%')
end

% Delete connections
delete('s'); clear('s'); instrreset;


%% Calculate angle coefficients and plot map
if calib_fhp == 1
    % Calculate probe coefficients and Reynolds number
    c = exp_fhp_calibrate(Iota,Tau,P_raw,probe.(probe_name),Ta,'hodson',1);
end


%% Calculate temperature data
if calib_top == 1
    c.T_raw = T_raw;
end


%% Record calibrations 

% Search for unique filename
n = 1;
file_prefix = [probe_name '_' calib_name '_' date '_' machine_name];
while exist([file_prefix '_' num2str(n) '.mat'],'file') ~= 0
    n = n+1;
end

% Save figure window
if ishandle(h) == 1
    saveas(h,[file_prefix '_' num2str(n) '.fig'],'fig');
end

% Save calibration data
c.Iota = Iota; c.Tau = Tau;
c.P_raw = P_raw; c.Pa = Pa; c.Ta = Ta; c.Ha = Ha;
save([file_prefix '_' num2str(n) '.mat'],'c')


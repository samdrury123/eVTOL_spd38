% Run a set of characteristics on the Gibbons rig

% Clear workspace
close all; clear all; clc; instrreset;

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

N.V_shaft = 1;

% Length of data arrays to allocate for
run_time = 5; % Minutes

% Sector configuration
% N.P_4_sec = 1;
% N.P_5_sec = 2;
% chic_type = 'sector';
% log_daq = 0;

% Five hole probe incidence
% N.FHP = 1:5; probe_name = 'FHP_G'; N.P_5_sec = 15;
% calib_name = 'FHP_G_61_fine_05-Jun-2015_Bingo_3.mat'; load([directory calib_name]);
% chic_type = 'incidence'; Iota = 55; nk = 7; e.probe = probe.(probe_name);
% t_min = -2.5; t_max = 2.3;        % 75 Blades
% % t_min = -2; t_max = 3.14;         % 70 Blades
% % t_min = -2.5; t_max = 3.04;       % 65 Blades
% log_daq = 0;

% Normal characteristic
chic_type = 'annulus';
log_daq = 0;

% Torque characteristic
% chic_type = 'torque';
% log_daq = 1;
% N.V_torque = 5;
% pxie.rate = 10000;
% pxie.chan = [N.V_shaft N.V_torque];
% nv = run_time * 60 * pxie.rate;

% Kulite characteristic
% chic_type = 'kulite';
% log_daq = 1;
% N.V_torque = 2;
% N.V_kulite = 3:8;
% pxie.rate = 30000;
% pxie.chan = [N.V_shaft N.V_torque N.V_kulite];
% nv = run_time * 60 * pxie.rate;

% Kulite characteristic
% chic_type = 'kulite_row';
% log_daq = 1;
% N.V_kulite = 1:8;
% pxie.rate = 60000;
% pxie.chan = [N.V_kulite N.V_shaft];
% nv = run_time * 60 * pxie.rate;

% Overall channel numbers
P_chan = 16; T_chan = 4; 

% Current geometry
% geom_name = 'Machined_R12_C4';
% geom_name = 'Machined_R11_C4';

geom_name = 'SLA_S100_L100_N75'; 

% geom_name = 'SLA_S30_L30_N70'; 
% geom_name = 'SLA_S100_L100_N70'; 
% geom_name = 'SLA_S170_L170_N70'; 

% geom_name = 'SLA_S0_L0_N65'; 
% geom_name = 'SLA_S100_L100_N65'; 
% geom_name = 'SLA_S200_L200_N65'; 

% Preallocate data arrays
P_rate = 6; 
np = run_time * 60 * P_rate; 

% Select non-dimensional blade speed
nrt = 0.6823;

% Plotting options
plot_chic = 1; plot_phi = 0; t_delay = 200;
plot_statics = 1;


%% Initialise instruments

% Open thermocouple logger
tc08.nchan = T_chan; 
s.tc08 = exp_tc08_initialise(tc08);

% Open dsa connection and reduce average count
s.dsa = exp_dsa_initialise(dsa);
fprintf(s.dsa, 'SET AVG 20\n');

% Use stepper motors to position five hole probe
if strcmp(chic_type,'incidence') == 1
    % Start motors and home
    s.park = exp_stepper_initialise(park);
    exp_stepper_home(s.park,park.home,1)
    park.limit = park.limit.T4;
    
    % Calculate tangential positions
    t_temp = linspace(t_min,t_max,nk+1)'; t_temp = t_temp(1:end-1);
    r_temp = r_mid * ones(nk,1); Iota = Iota * ones(nk,1);
    
    % Calculate steps to centre of passage
    q = round([ (r_temp + probe.(probe_name).r_off) * park.d(1) + park.c(1) ...
        t_temp*park.d(2) + park.c(2) Iota*park.d(3) + park.c(3) ]);
    
    % Move motors to first position
    exp_stepper_coord(s.park,q(1,:),park.limit);
end

% Point to low speed lab PTH sensor
s.ptu = ptu;

% Preallocate plotting and data arrays
phi = zeros(np,1); psi = zeros(np,1); psi_rotor = zeros(np,1); psi_stator = zeros(np,1);
phi(1) = nan; psi(1) = nan; psi_rotor(1) = nan; psi_stator(1) = nan;
P = zeros(np,P_chan); T = zeros(np,T_chan); time_t = zeros(np,6); time_p = zeros(np,6);
P(1,:) = nan; T(1,:) = nan; time_t(1,:) = nan; time_p(1,:) = nan;

% Initialise coordinate arrays
if strcmp(chic_type,'incidence') == 1
    t = zeros(np,1); t(1,:) = nan;
end

% Preallocate global NI-DAQ variables
if log_daq == 1
    global V time_v 
    V = [];
    time_v = [];
end

%% Get atmospheric pressure and humidity
[Pa,~,Ha] = getCurrentPTH(ptu.loc);

%% Open plotting window for overall characteristic
if plot_chic == 1
    
% Open figure and label axes
h.window = figure(); set(h.window,'Position',[1 29 1920 985]); 
h.axes = axes('Position',[0.13 0.11 0.775 0.815]); hold on; grid on; box on;
xlabel('Flow Coefficient','FontSize',16); ylabel('Pressure Rise Coefficient','FontSize',16); 
axis([0.38 0.7 0 0.45]);

% Plot first lines
l.psi = plot(phi,psi,'k-','LineWidth',2); 
if plot_statics == 1
    l.psi_rotor = plot(phi,psi_rotor,'r-'); l.psi_stator = plot(phi,psi_stator,'b-');
    legend('Overall T-S','Rotor S-S','Stator S-S','Location','NorthEast','FontSize',16);
end

% Plot markers
l.psi_mark = plot(phi(1),psi(1),'ko','MarkerFaceColor','k','MarkerSize',10);
if plot_statics == 1
    l.psi_rotor_mark = plot(phi(1),psi_rotor(1),'ro','MarkerFaceColor','r','MarkerSize',7);
    l.psi_stator_mark = plot(phi(1),psi_stator(1),'bo','MarkerFaceColor','b','MarkerSize',7);
end

% % Set automatic properties to manual
% pn = {'ALimMode','CameraPositionMode','CameraTargetMode',...
%     'CameraUpVectorMode','CLimMode','TickDirMode','XLimMode',...
%     'YLimMode','ZLimMode','XTickMode','YTickMode',...
%     'ZTickMode','XTickLabelMode','YTickLabelMode','ZTickLabelMode'};
% pv = cell(size(pn)); for n = 1:length(pn); pv(n) = {'manual'}; end;
% set(h.axes,pn,pv);

end


%% Open plotting window for raw pressures and flow coefficient
if plot_phi == 1
    
% Open flow coefficient setting window
h.window = figure(); set(h.window,'Position',[1 29 1920 985]);
h.flow = axes('Position',[0.13 0.4 0.775 0.55]); hold on; grid on; box on; grid minor;
xlabel('Time / s'); ylabel('Flow Coef.');  axis([-t_delay 0 phi_min-0.01 phi_max+0.01]);
h.P = axes('Position',[0.13 0.1 0.775 0.25]); hold on; grid on; box on;
xlabel('Time / s'); ylabel('Raw Pressure');  axis([-t_delay 0 -5000 10000]);

% Initialise plotting variables
time_diff = -200*ones(np,1);

% Pressure and temperature legends
P_legend = cell(P_chan,1); T_legend = cell(T_chan,1);
for n = 1:P_chan; P_legend{n} = num2str(n); end;

% Plot initialised values
axes(h.flow); l.flow = line(time_diff,phi);
axes(h.P); l.P = line(time_diff,P);
legend(P_legend,'Location','West','FontSize',6);
    
end

% Log pxie zeros
if log_daq == 1
    pxie.time = 1;
    s.pxie = exp_pxie_initialise(pxie,'DC');
    V_zero = exp_pxie_read(s.pxie);
end


%% Start rig and begin logging

uiwait(msgbox('Start Rig'));

% Set non-dimensional speed
pxie_temp = pxie; pxie_temp.chan = N.V_shaft; pxie_temp.rate = 200e3; pxie_temp.time = 0.2;
s.pxie = exp_pxie_initialise(pxie_temp,'DC');
des.nrt = nrt; des.mrtp = 1;
err.phi = 0.0025; err.nrt = 0.0025; err.mrtp = 0.0025; 
N_rig.P = N.P_1_cas; N_rig.T = N.To_2_mid(2); N_rig.Po = N.Po_1; N_rig.Pa = N.Pa; N_rig.V = 1;
e_rig = exp_rig_set(s,N_rig,des,err,r_mid);

% Calculate blade speed
e.w = e_rig.w; e.nrt = e_rig.nrt;
U_mid = e_rig.w * r_mid; e.U_mid = U_mid;

% Open NI-DAQ connection for chic measurements
if log_daq == 1
    s.pxie = exp_pxie_initialise(pxie,'DC');
    s.pxie.IsContinuous = true;
    s.pxie.NotifyWhenDataAvailableExceeds = round(pxie.rate/5);

    % Start a listener with global variables
    lh = addlistener(s.pxie,'DataAvailable',@exp_pxie_background);
end

% Set stop button
g = 0; uicontrol('style','push','string','Finished','callback','g=g+1');
pause(0.5);

% Put in motor positioning button
if strcmp(chic_type,'incidence') == 1
    k = 1; k_old = 1;
    uicontrol('style','push','Position',[100 20 60 20],...
        'string','Next Pos','callback','k=k+1');
    pause(0.5);
end

% Calculate DSA and thermocouple period
t_sample = 1/P_rate;

% Begin logging from pxie
if log_daq == 1
    startBackground(s.pxie);
end

% Log until stop request
while g == 0 
    
    % Begin timing
    tic;
    
    % Move motors to next position
    if strcmp(chic_type,'incidence') == 1 && k ~= k_old    
        if k > nk
            break
        end
        exp_stepper_coord(s.park,q(k,:),park.limit);
        k_old = k;
    end
    
    % Record coordinates
    if strcmp(chic_type,'incidence') == 1
        t = [t_temp(k) ; t(1:end-1)];
    end
    
    % Log time and pressure
    time_p_read = clock;
    P_read = exp_dsa_read(s.dsa);
    
    % Log time and temperature
    time_t_read = clock;
    T_read = exp_tc08_read(s.tc08,1);
    
    % Replace temperature with old value if new reading not available
    if sum(isnan(T_read)) == 0
        T_old = T_read;
        time_old = time_t_read;
    elseif sum(isnan(T_read)) > 0
        T_read = T_old;
        time_t_read = time_old;
    end   
    
    % Record raw data
    P = [P_read ; P(1:end-1,:)];
    T = [T_read ; T(1:end-1,:)];
    time_t = [time_t_read ; time_t(1:end-1,:)];
    time_p = [time_p_read ; time_p(1:end-1,:)];
    if plot_phi == 1
        for n = 1:length(time_diff)
            time_diff(n,:) = etime(time_t(n,:),time_t_read);
        end
    end
    
    % Calculate quick flow and pressure rise coefficients
    Po_1 = P_read(N.Po_1); P_2_cas = P_read(N.P_2_cas); 
    P_3_cas = P_read(N.P_3_cas); P_6_cas = P_read(N.P_6_cas);
    T_2 = mean(T_read(N.To_2_mid));
    if strcmp(chic_type,'sector') == 1
        P_4_cas = P_read(N.P_4_sec); P_5_cas = P_read(N.P_5_sec);
    else
        P_4_cas = P_read(N.P_4_cas); P_5_cas = P_read(N.P_5_cas);
    end
    
    ro_2 = (P_2_cas + Pa) / (air.R * T_2);

    phi_new = real((2 * (Po_1 - P_2_cas) / ro_2).^0.5 / U_mid);
    psi_new = (P_5_cas - Po_1) / (ro_2 * U_mid^2);
    psi_rotor_new = (P_4_cas - P_3_cas) / (ro_2 * U_mid^2);
    psi_stator_new = (P_5_cas - P_4_cas) / (ro_2 * U_mid^2);
    
    phi = [phi_new ; phi(1:end-1,:)];
    psi = [psi_new ; psi(1:end-1,:)];
    psi_rotor = [psi_rotor_new ; psi_rotor(1:end-1,:)];
    psi_stator = [psi_stator_new ; psi_stator(1:end-1,:)];
    
    % Update characteristic plot
    if plot_chic == 1
        set(l.psi,'XData',phi,'YData',psi);
        set(l.psi_mark,'XData',phi(1),'YData',psi(1));
        if plot_statics == 1
            set(l.psi_rotor,'XData',phi,'YData',psi_rotor);
            set(l.psi_stator,'XData',phi,'YData',psi_stator);
            set(l.psi_stator_mark,'XData',phi(1),'YData',psi_stator(1));
            set(l.psi_rotor_mark,'XData',phi(1),'YData',psi_rotor(1));
        end
    end
    
    % Update raw pressures and flow coefficient plot
    if plot_phi == 1
        set(l.flow,'XData',time_diff,'YData',phi);
        for n = 1:length(l.P)
            set(l.P(n),'XData',time_diff,'YData',P(:,n));
        end
    end
    
    % Update finished button
    drawnow;
    
    % Build in delay for constant sample rate
    t_current = toc;
    if t_current < t_sample
        pause(t_sample - t_current);
    end
end

% Stop pxie logging
if log_daq == 1
    stop(s.pxie);
end

%% Disconnect and stop rig

% Move probe back to home
if strcmp(chic_type,'incidence') == 1
    exp_stepper_coord(s.park,[0 0 0],park.limit);
    e.Iota = Iota; e.probe_name = probe_name; e.c = c; e.r = r_temp;
    
    % Close steppers
    for n = 1:3
        exp_stepper_command(s.park,n,'OFF');
        exp_stepper_command(s.park,n,'E0');
    end
    fclose(s.park); delete(s.park); s = rmfield(s,'park');
end

% Close instruments
uiwait(msgbox('Stop Rig'));

% Retake zero, set average back and disconnect dsa
fprintf(s.dsa, ['SET AVG ' num2str(dsa.avg) '\n']);
fclose(s.dsa); s = rmfield(s,'dsa');

% Close TC08 thermocouples
calllib('usbtc08','usb_tc08_close_unit',s.tc08.h);
s = rmfield(s,'tc08');

delete('s'); clear('s'); instrreset;

%% Rearrange data before saving

% Remove nans
q = find(isnan(P(:,1)),1,'first'):size(P,1);
P(q,:) = []; T(q,:) = []; time_t(q,:) = []; time_p(q,:) = [];
if log_daq == 1
    q = find(isnan(V(:,1)),1,'first'):size(V,1);
    V(q,:) = []; time_v(q,:) = [];
end

% Flip into increasing time order
P = flipud(P);
T = flipud(T);
time_t = flipud(time_t);
time_p = flipud(time_p);
if log_daq == 1
    time_v = flipud(time_v);
    V = flipud(V);
end

% Operate on theta coordinates
if strcmp(chic_type,'incidence') == 1
    t(q,:) = [];
    t = flipud(t);
    e.t = t;
end


%% Save data

% Save raw pressure and temperature data
e.P = P; e.T = T; e.time_t = time_t; e.time_p = time_p;

% Save raw voltage and time data
if log_daq == 1
    e.V = V;
    e.time_v = time_v;
end

% Save atmospheric information
e.Pa = Pa; e.Ha = Ha;

% Save configuration information
e.N = N; e.geom_name = geom_name;

% Save kulite zeros
if log_daq == 1
    e.V_zero = V_zero;
end

% Save to mat file - no overwritting
n = 1;
file_prefix = [geom_name  '_' chic_type '_' date];
while exist([directory file_prefix '_' num2str(n) '_chic.mat'],'file') ~= 0
    n = n+1;
end
save([directory file_prefix '_' num2str(n) '_chic.mat'],'e')
%   Exp_Characteristic                              
%                               
%   A script to run a characteristic on any rig
%
%   Written by James Taylor                
%   October 2016
%
%   Choose type of call with chic_type - 'annulus', 'sector', 'incidence', 'torque' or 'kulite'

% Default to not demo
if exist('demo_mode','var') == 0
    demo_mode = 0;
end


%% Initialise instruments

% Open thermocouple logger
s.tc08 = exp_tc08_initialise(tc08);

% Open dsa connection and reduce average count
s.dsa = exp_dsa_initialise(dsa);
fprintf(s.dsa, 'SET AVG 20\n');

% Open connection to read voltages
pxie.time = 0.2; s.pxie = exp_pxie_initialise(pxie);

% Open connection to read RTDs
rtd.time = 0.2; s.rtd = exp_pxie_initialise(rtd);

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
if isempty(timerfind) == 1; PTUstart(ptu.loc); end
s.ptu = ptu;

% Preallocate plotting and data arrays
phi = zeros(np,1); psi = zeros(np,1); phi(1) = nan; psi(1) = nan;
P = zeros(np,dsa.nchan); Q = zeros(np,length(tc08.chan)); time_q = zeros(np,6); time_p = zeros(np,6);
P(1,:) = nan; Q(1,:) = nan; time_q(1,:) = nan; time_p(1,:) = nan;

% Initialise coordinate arrays
if strcmp(chic_type,'incidence') == 1
    t = zeros(np,1); t(1,:) = nan;
end

% Preallocate global NI-DAQ variables
global V time_v T time_t
V = []; time_v = []; T = []; time_t = []; 


%% Get atmospheric pressure and humidity

% Get atmospherics locally or from Vaisala box
if strcmp(N_rig.PTH,'Vaisala') == 1
%     if isempty(timerfind) == 1; PTUstart(ptu.loc); end;    
    [Pa,Ta,Ha] = getCurrentPTH(ptu.loc);
end


%% Open plotting window for overall characteristic

% Open figure and label axes
h.window = figure();
if demo_mode == 0
    h.axes = axes('Position',[0.13 0.11 0.775 0.815]); hold on; grid on; box on;
    set(h.window,'Position',[1 29 1920 985]); axis([0.36 0.67 0 0.45]);
    xlabel('Flow Coefficient','FontSize',16); ylabel('Pressure Rise Coefficient','FontSize',16); 
else
    h.axes = axes('Position',[0.13 0.11 0.775 0.815]); hold on; grid on; box on;
    
    set(h.window,'Position',[-3839 -1079 3840 2084]); axis([0.13 0.6  0.068  0.32]);
    xlabel('Flow Rate','FontSize',48,'color','k'); ylabel('Pressure Rise','FontSize',48,'color','k'); 
end
set(gca,'xticklabel','','yticklabel','');

% Plot first line
col = lines(1);
l.psi = plot(phi,psi,'-','LineWidth',4,'color',col); 

% Plot markers
l.psi_mark = plot(phi(1),psi(1),'ko','MarkerFaceColor','k','MarkerSize',20);


%% Start rig and begin logging

% Log zero voltages for one second
s.pxie.DurationInSeconds = 1; V_zero = exp_pxie_read(s.pxie); s.pxie.DurationInSeconds = 0.2; 

% Turn the rig on
uiwait(msgbox('Start Rig'));

% Set non-dimensional speed by measurment or default for demonstrations
if demo_mode == 0
    
    % Add the listeners to monitor the rig
    s.pxie.NotifyWhenDataAvailableExceeds = s.pxie.NumberOfScans;
    s.rtd.NotifyWhenDataAvailableExceeds = s.rtd.NumberOfScans;
    lv = addlistener(s.pxie,'DataAvailable',@exp_pxie_background);
    lt = addlistener(s.rtd,'DataAvailable',@exp_pxie_background);
    
    % Set the rig operating point
    des.nrt = nrt; des.mrtp = 1;
    err.phi = 0.0025; err.nrt = 0.0025; err.mrtp = 0.0025; 
    e_rig = exp_rig_set(s,N_rig,des,err,r_mid);
    
    % Delete the listeners and re-initialise the arrays
    delete(lv); delete(lt); 
    V = []; time_v = []; T = []; time_t = []; 
    
    % Calculate the blade speed and record
    e.w = e_rig.w; e.nrt = e_rig.nrt;
    U_mid = e_rig.w * r_mid; e.U_mid = U_mid;
else
    
    % Select default values of operating point
    e.w = 366; e.nrt = 0.6823; U_mid = e.w * r_mid; e.U_mid = U_mid;
end

% Set continuous voltage and temperature measurements
s.pxie.IsContinuous = true; s.pxie.NotifyWhenDataAvailableExceeds = round(pxie.rate/5);
s.rtd.IsContinuous = true; s.rtd.NotifyWhenDataAvailableExceeds = round(rtd.rate/5);

% Start a listener with global variables
addlistener(s.pxie,'DataAvailable',@exp_pxie_background);
addlistener(s.rtd,'DataAvailable',@exp_pxie_background);

% Set stop button
g = 0; uicontrol('style','push','string','Finished','callback','g=g+1');
pause(0.5);

% Put in motor positioning button
if strcmp(chic_type,'incidence') == 1
    k = 1; k_old = 1;
    uicontrol('style','push','Position',[100 20 60 20],'string','Next Pos','callback','k=k+1');
    pause(0.5);
end

% Calculate DSA and thermocouple period
t_sample = 1/P_rate;

% Begin logging from pxie
startBackground(s.pxie); startBackground(s.rtd);
pause(0.2);
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
    time_q_read = clock;
    Q_read = exp_tc08_read(s.tc08,1);
    
    % Replace temperature with old value if new reading not available
    if sum(isnan(Q_read)) == 0
        Q_old = Q_read;
        time_old = time_q_read;
    elseif sum(isnan(Q_read)) > 0
        Q_read = Q_old;
        time_q_read = time_old;
    end   
    
    % Record raw data
    P = [P_read ; P(1:end-1,:)];
    Q = [Q_read ; Q(1:end-1,:)];
    time_q = [time_q_read ; time_q(1:end-1,:)];
    time_p = [time_p_read ; time_p(1:end-1,:)];
    
    % Extract measurements for performance calculations
    Po_in = P(1,N_plot.Po_in); P_in = P(1,N_plot.P_in); P_out = P(1,N_plot.P_out); P_ref = P(1,N_rig.Pa); 
    if strcmp(N_plot.T_type,'rtd') == 1
        T_in = mean(T(1,N_plot.T_in)); 
    else
        T_in = mean(Q(1,N_plot.T_in)); 
    end
    
    % Calculate rough inlet density
    ro = (P_in + Pa - P_ref) / (air.R * T_in);

    % Calculate rough pressure rise and flow coefficients
    phi_new = real((2 * (Po_in - P_in) / ro).^0.5 / U_mid);
    psi_new = (P_out - Po_in) / (ro * U_mid^2);
    
    % Update flow and pressure rise coefficient vectors
    phi = [phi_new ; phi(1:end-1,:)];
    psi = [psi_new ; psi(1:end-1,:)];
    
    % Update characteristic plot
    set(l.psi,'XData',phi,'YData',psi);
    set(l.psi_mark,'XData',phi(1),'YData',psi(1));
    
    % Update finished button
    drawnow;
    
    % Build in delay for constant sample rate
    t_current = toc;
    if t_current < t_sample
        pause(t_sample - t_current);
    end
end

% Stop pxie logging
stop(s.pxie); stop(s.rtd);


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

% Stop PTU logger
if strcmp(N_rig.PTH,'Vaisala') == 1; delete(timerfind); end;

delete('s'); clear('s'); instrreset;


%% Rearrange data before saving

% Remove nans
q = find(isnan(P(:,1)),1,'first'):size(P,1);
P(q,:) = []; Q(q,:) = []; time_q(q,:) = []; time_p(q,:) = [];
q = find(isnan(V(:,1)),1,'first'):size(V,1);
V(q,:) = []; time_v(q,:) = [];
q = find(isnan(T(:,1)),1,'first'):size(T,1);
T(q,:) = []; time_t(q,:) = [];

% Flip into increasing time order
P = flipud(P); Q = flipud(Q); V = flipud(V); T = flipud(T);
time_p = flipud(time_p); time_q = flipud(time_q); time_v = flipud(time_v); time_t = flipud(time_t);

% Operate on theta coordinates
if strcmp(chic_type,'incidence') == 1
    t(q,:) = [];
    t = flipud(t);
    e.t = t;
end


%% Save data

% Save raw pressure and temperature data
e.P = P; e.Q = Q; e.T = T; e.time_t = time_t; e.time_p = time_p; e.time_q = time_q;

% Save raw voltage and time data
e.V = V; e.time_v = time_v;

% Save atmospheric information
e.Pa = Pa; e.Ha = Ha; e.Ta = Ta;

% Save configuration information
e.N = N; e.geom_name = geom_name;

% Save kulite zeros
e.V_zero = V_zero;

% Save to mat file - no overwritting
n = 1;
file_prefix = [geom_name  '_' chic_type '_' date];
while exist([dr.exp file_prefix '_' num2str(n) '_chic.mat'],'file') ~= 0
    n = n+1;
end
save([dr.exp file_prefix '_' num2str(n) '_chic.mat'],'e')


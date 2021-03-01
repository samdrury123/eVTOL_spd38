% Calibrate yaw of hotwire in bingo machine and speed in Gibbons rig

clear; close all; clc; instrreset; delete(instrfindall);

iGibbons

% Choose calibration type
% calib_type = 'yaw_gibbons';
% calib_type = 'yaw_bingo';
calib_type = 'velocity';

% Probe name
probe_name = 'HW_E';

% In Gibbons tests, push until dot on holder is ~5 mm outside of aluminium collet


%% Calibrate hot wire for yaw angle offset in Gibbons bellmouth
if strcmp(calib_type,'yaw_gibbons') == 1
    
    % Yaw settings - negative direction is anti-clockwise from zero
    Iota_min = -40; Iota_max = 40; ni = 9;
    
    % Hotwire settings
    pxie.nchan = 1; pxie.chan = probe.(probe_name).N.V_hw_av; pxie.time = 0.5;
    
    % Start PXIE
    s.pxie = exp_pxie_initialise(pxie);
    
    % Define angle and initialise voltage vector
    c.Iota = linspace(Iota_min,Iota_max,ni)';
    c.V = zeros(ni,pxie.time * pxie.rate + 1); c.V_mean = zeros(ni,1);
        
    % Open figure window to show readings
    figure(); hold on; grid on; box on; xlabel('Yaw Angle / ^\circ'); ylabel('Mean Voltage / V');
    cols = lines(3);
    
    % Turn on the rig
    uiwait(msgbox('Operate wire and turn rig on'));
    
    % Loop over all calibration points and record voltages
    for n = 1:length(c.Iota)
        
        % Move stepper motors to next coordinates
        uiwait(msgbox(['Move wire to Iota = ' num2str(c.Iota(n))]));
        
        % Read hotwire voltage
        [V,t] = exp_pxie_read(s.pxie);
        c.V(n,:) = V;
        
        % Plot average voltage
        c.V_mean(n) = trimmean(V,10);
        plot(c.Iota(n),c.V_mean(n),'.','color',cols(1,:),'markersize',12);
    end
    
    % Calculate quadratic fit on symmetric points only
    V_min_1 = min(c.V_mean(1:round(ni/2))); V_min_2 = min(c.V_mean(round(ni/2):end));
    q = c.V_mean >= max(V_min_1,V_min_2); p = polyfit(c.Iota(q),c.V_mean(q),2);
    
    % Plot polynomial fit
    Iota_fit = linspace(Iota_min,Iota_max,100);
    plot(Iota_fit,polyval(p,Iota_fit),'-','color',cols(2,:));
    
    % Calculate wire angle offset
    c.Iota_off = -p(2) / (2 * p(1));
    plot(c.Iota_off,polyval(p,c.Iota_off),'.','color',cols(3,:));
    
    % Disconnect instruments
    delete(s.pxie); clear('s'); instrreset; delete(instrfindall);   
    
end


%% Calibrate hot wire for yaw angle offset in bingo machine
if strcmp(calib_type,'yaw_bingo') == 1
    
    % Yaw motor settings
    park.current = 0.6; park.resistance = 4.2; park.inductance = 2.5;
    limit = [-72000 72000]; park.nchan = 1;
    Iota_d = 200; Iota_c = -36000;

    % Hotwire settings
    Iota_min = -0; Iota_max = 20; ni = 41;
    pxie.nchan = 1; pxie.chan = 3; pxie.time = 0.5;
    
    % Get atmospheric data
    [c.Pa,c.Ta,c.Ha] = getCurrentPTH(ptu.loc);
    
    % Start stepper drivers
    s.park = exp_stepper_initialise(park,0);
    
    % Start DSA
    s.dsa = exp_dsa_initialise(dsa);
    
    % Start PXIE
    s.pxie = exp_pxie_initialise(pxie,'DC');
    
    % Define angle matrices
    c.Iota = linspace(Iota_min,Iota_max,ni)';
    c.Iota = [c.Iota ; flipud(c.Iota)];
    c.q = round((c.Iota + probe.(probe_name).Iota_off) * Iota_d + Iota_c);
    
    % Initialise voltage and pressure matrices
    V = single(zeros(length(c.Iota),pxie.time*pxie.rate+1));
    P = zeros(length(c.Iota),2);
    
    % Record no flow voltage
    uiwait(msgbox('Operate Wire'));
    V_temp = exp_pxie_read(s.pxie);
    e.V0_start = single(V_temp);
    
    % Turn on bingo machine
    uiwait(msgbox('Turn Tunnel On'));
        
    % Loop over all calibration points and record voltages and pressures
    for n = 1:length(c.Iota)
        
        % Move stepper motors to next coordinates
        exp_stepper_coord(s.park,c.q(n),limit);
        pause(0.1);
        
        % Read hotwire voltage
        [V,t] = exp_pxie_read(s.pxie);
        c.V(n,:) = V;
        
        % Read pressures
        P = exp_dsa_read(s.dsa);
        c.P(n,:) = P(1:2);
    end
    
    % Yaw back to 0
    exp_stepper_coord(s.park,0,limit);
    
    % Turn off bingo machine
    uiwait(msgbox('Turn Tunnel Off'));
    
    % Record no flow voltage
    V_temp = exp_pxie_read(s.pxie);
    e.V0_end = single(V_temp);
    uiwait(msgbox('Standby Wire'));
    
    % Disconnect instruments
    fclose(s.dsa); s = rmfield(s,'dsa');
    
    exp_stepper_command(s.park,1,'OFF');
    exp_stepper_command(s.park,1,'E0');
    fclose(s.park); s = rmfield(s,'park');
    
    delete(s.pxie);
    
    clear('s'); instrreset; delete(instrfindall);   
end


%% Calibrate hot wire for velocity variations
if strcmp(calib_type,'velocity') == 1

    % Initialise and zero DSA
    s.dsa = exp_dsa_initialise(dsa);

    % Initialise once per rev and hot wire logging
    N.V_hw = probe.(probe_name).N.V_hw_av;
    pxie.nchan = 2; pxie.chan = [N.V_shaft N.V_hw]; pxie.time = 1;
    s.pxie = exp_pxie_initialise(pxie);
    
    % Initialise thermocouples
    s.rtd = exp_pxie_initialise(rtd);

    % Get atmospheric data
    [c.Pa,c.Ta,c.Ha] = getCurrentPTH(ptu.loc);

    % Define angle and speed range
    rpm = linspace(0,3600,10)';
    nn = length(rpm);

    % Pre-allocate arrays
    c.time = zeros(nn,6);
    c.P = zeros(nn,16);
    c.V = zeros(nn,pxie.time * pxie.rate + 1,2);

    % Figure window to show progress
    figure(); hold on; grid on; box on;
    xlabel('Pressure Difference / Pa'); ylabel('Voltage / V');
    
    % Turn wire on
    uiwait(msgbox('Operate Wire'));

    % Collect data
    for n = 1:length(rpm)
        
        % Set the rig speed
        uiwait(msgbox(['Set RPM = ' num2str(rpm(n))]));

        % Log data
        t_temp = clock;
        P_temp = exp_dsa_read(s.dsa);
        [V_temp,t_v_temp] = exp_pxie_read(s.pxie);
        T_temp = exp_pxie_read(s.rtd);

        % Record in struct
        c.time(n,:) = t_temp;
        c.P(n,:) = P_temp;
        c.V(n,:,:) = V_temp;
        c.dtime = t_v_temp;
        c.T(n,:) = mean(T_temp,1);

        % Plot the current point
        plot(P_temp(N.Po_bell) - P_temp(N.P_bell),mean(V_temp(:,2)),'b.','markersize',12)
    end

    % Disconnect dsa
    fclose(s.dsa); delete(s.dsa); s = rmfield(s,'dsa');
    
    % Turn wire off
    uiwait(msgbox('Standby Wire'));

    % Close PXIE
    delete(s.pxie);

    % Record standard and over temperatures
    c.T_std = c.Ta; c.T_over = probe.(probe_name).ohr / probe.(probe_name).alpha;
    
    % Calculate calbiration coefficients
    N.To = 1; N.P = N.P_bell; N.Po = N.Po_bell;
    c = exp_hw_calibrate(c,N,1,0.96);

end


% Save calibration data in a unique filename
n = 1;
file_prefix = [probe_name '_' calib_type '_' date];
while exist([dr.exp file_prefix '_' num2str(n) '.mat'],'file') ~= 0
    n = n+1;
end
save([dr.exp file_prefix '_' num2str(n) '.mat'],'c')


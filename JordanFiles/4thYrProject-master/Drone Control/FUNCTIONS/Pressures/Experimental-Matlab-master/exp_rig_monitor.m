function [e_rig,on_point,e_read] = exp_rig_monitor(s,e_rig,N,des,err,r_mid)
% EXP_RIG_MONITOR  Check if a rig is on the right operating point
%
%   on_point = EXP_RIG_MONITOR(s,N,des,err,r_mid)
%   
%   s - data structure with open instrument connections
%   e_rig - past history of operating points
%   N - data structure of required channels
%   des - data structure of desired operating point
%   err - data structure of allowable error fractions
%   r_mid - midheight radius for flow coefficient
%   e_read - read instantaneous readings of raw data

% Initialise on_point flag
on_point = 1;

% Decide whether to check moving average of rig operating point
if isempty(e_rig) == 0
    check_avg = 1; n_avg = 9;
else
    check_avg = 0;
end

% Initialise voltage and temperature variables
global V time_v T time_t
V = []; time_v = []; T = []; time_t = []; 


%% Get all raw rig readings

% Get rig logging time
time = clock;

% Get atmospheric readings from Vaisala PTU sensor
if strcmp(N.PTH,'Vaisala') == 1
    [Pa,Ta,Ha] = getCurrentPTH(s.ptu.loc);
elseif strcmp(N.PTH,'Local') == 1
    [Pa,Ta,Ha] = exp_read_PTH(s,N);
end

% Start scan on DSA
exp_dsa_read(s.dsa,'scan');

% Start logging NI-DAQ in background
if isfield(s,'rtd') == 1; stop(s.rtd); startBackground(s.rtd); end;
stop(s.pxie); startBackground(s.pxie);

% Get pressure readings from DSA
P = exp_dsa_read(s.dsa,'read');

% Wait for temperature readings to arrive
while isempty(T) == 1; pause(0.01); end;

% Get temperature readings from TC08
if isfield(s,'rtd') == 0
    T = exp_tc08_read(s.tc08);
else
    T = [mean(T,1) exp_tc08_read(s.tc08)];
end

% Wait for voltage readings to arrive
while isempty(V) == 1; pause(0.01); end;
V = flip(V,1); dtime = linspace(0,s.pxie.DurationInSeconds,size(V,1))';

% Get shaft frequency
w = exp_shaft_freq(s.pxie,N.V,0,V,dtime);

% Record raw data
e_read.P = P; e_read.T = T; e_read.Pa = Pa; e_read.Ta = Ta; e_read.Ha = Ha; 
e_read.w = w; e_read.N = N; e_read.time = time; e_read.dtime = dtime; e_read.V = V;


%% Compare current operating point to desired values

% Calculate rig operating point
e_read = exp_rig_calculate(e_read,r_mid);

% Update rig data or record single instance
if isempty(e_rig) == 0
    rignames = [{'P' 'T' 'Pa' 'Ta' 'Ha' 'w' 'time'} fieldnames(des)'];
    for v = 1:length(rignames)
        e_rig.(rignames{v}) = [e_rig.(rignames{v}) ; e_read.(rignames{v})];
    end
else
    e_rig = e_read;
end
    
% Check moving average or not
if check_avg == 0
    if isfield(des,'phi') == 0
        nrt_check = e_rig.nrt(end);
        mrtp_check = e_rig.mrtp(end);
    else
        phi_check = e_rig.phi(end);
    end
else
    
    % Average compressibly or incompressibly
    if isfield(des,'phi') == 0
        
        % Get all reasonable last points to average over
        q = length(e_rig.nrt)- min(n_avg,length(e_rig.nrt))+1 : length(e_rig.nrt);
        q(e_rig.nrt(q) < 0.95 * mean(e_rig.nrt(q))) = [];
        q(e_rig.nrt(q) > 1.05 * mean(e_rig.nrt(q))) = [];

        % Average values
        nrt_check = trimmean(e_rig.nrt(q),15);
        mrtp_check = trimmean(e_rig.mrtp(q),15);
    else
        
        % Average flow coefficient
        phi_check = trimmean(e_rig.phi(end-min(n_avg,length(e_rig.phi))+1:end),15);
    end
end

% Determine whether to check compressible flow function or flow coefficient
if isfield(des,'phi') == 0
    
    % Calculate error in desired values
    mrtp_error = abs((mrtp_check - des.mrtp)/des.mrtp);
    nrt_error = abs((nrt_check - des.nrt)/des.nrt);
    mrtp_last_error = abs((e_rig.mrtp(end) - des.mrtp)/des.mrtp);
    nrt_last_error = abs((e_rig.nrt(end) - des.nrt)/des.nrt);
    
    % Check inlet flow function and non-dimensional blade speed
    if (mrtp_error > err.mrtp && mrtp_last_error > err.mrtp) || ...
            (nrt_error > err.nrt && nrt_last_error > err.nrt)
        on_point = 0;
    end

else
    
    % Error in flow coefficient
    phi_error = abs((phi_check - des.phi)/des.phi);
    phi_last_error = abs((e_rig.phi(end) - des.phi)/des.phi);
    
    % Check flow coefficient
    if phi_error > err.phi && phi_last_error > err.phi
        on_point = 0;
    end
end




end


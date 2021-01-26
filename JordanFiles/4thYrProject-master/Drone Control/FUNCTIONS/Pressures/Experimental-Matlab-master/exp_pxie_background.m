function [] = exp_pxie_background(src,event,log_single)
% Collect data from a background event running on a NI-DAQ
    
% Default to log all data
if exist('log_single','var') == 0
    log_single = 0;
end

% Global time, temperature and voltage arrays
global V time_v T time_t
persistent time_temp V_rate T_rate

% Record streaming voltages
if strcmp(src.Channels(1).MeasurementType,'Voltage') == 1 && log_single == 0

    % Concatenate voltage array
    V_rate = src.NotifyWhenDataAvailableExceeds;
    V = [flipud(event.Data) ; V];

    % Calculate exact time vectors for each sample
    time_temp = repmat(datevec(event.TriggerTime),[V_rate 1]);
    time_temp(:,6) = time_temp(:,6) + flipud(event.TimeStamps);
    time_v = [datevec(datenum(time_temp)) ; time_v];
    
elseif strcmp(src.Channels(1).MeasurementType,'RTD') == 1 && log_single == 0

    % Concatenate temperature array
    T_rate = src.NotifyWhenDataAvailableExceeds;
    T = [flipud(event.Data) ; T];

    % Calculate exact time vectors for each sample
    time_temp = repmat(datevec(event.TriggerTime),[T_rate 1]);
    time_temp(:,6) = time_temp(:,6) + flipud(event.TimeStamps);
    time_t = [datevec(datenum(time_temp)) ; time_t];
    
elseif strcmp(src.Channels(1).MeasurementType,'Voltage') == 1 && log_single == 1

    % Get voltage array
    V = flipud(event.Data);
    
elseif strcmp(src.Channels(1).MeasurementType,'RTD') == 1 && log_single == 1

    % Get temperature array
    T = flipud(event.Data);
    
end

end
function [T,T_all,t_all] = exp_tc08_read(s,skip_fail)
% EXP_TC08_READ  Read temperatures continuously from TC08 logger
%
%   [T,T_all,t_all] = EXP_TC08_READ(s,tc08,nchan,skip_fail)
%
%   s - input data structure detailing USB connection
%   skip_fail - flag an error if temperature fails to read
%   T - output temperature vector
%   T_all - output all temperatures read in this function call
%   t_all - output all time values at each reading in this function call


if exist('skip_fail','var') == 0
    skip_fail = 0;
end

% Initialise pointers to contain data
buffer_size = 600;
data_pointer = libpointer('singlePtr', NaN(buffer_size, 1));
time_pointer = libpointer('longPtr', NaN(buffer_size, 1));
overflow_pointer = libpointer('int16Ptr', 0);

% Preallocate data arrays
T = nan(1,s.nchan);
T_all = cell(1,s.nchan);
t_all = cell(1,s.nchan);

% Loop over all channels and get temperature readings (data is FIFO) so if
% it hasn't been called for a while old data will be overwritten as the
% buffer fills up
for n = int16(1:s.nchan)
    % Preallocate nans for time and temperature
    T_all{n} = nan(1,buffer_size);
    t_all{n} = nan(1,buffer_size);
    
    % Retrieves data into pointers
    read_flag = calllib('usbtc08', 'usb_tc08_get_temp', s.h, data_pointer,...
        time_pointer, buffer_size, overflow_pointer, n, s.units, int16(0));
    
    % Check for a successful read
    if read_flag == 0 && skip_fail == 0
        error('Temperature data not read')
    elseif read_flag == 0 && skip_fail == 1
        return
    end
    
    % Store data from pointers
    T_all{n} = get(data_pointer, 'Value');
    t_all{n} = get(time_pointer, 'Value');
    
%      % Store last temperature
       T_good = T_all{n}(isnan(T_all{n}) == 0);
       T(n) = T_good(end);
end

end

%% Single call to TC08

% % Define the pointers for the data
% data_pointer = libpointer('singlePtr', nan(tc08.nchan+1, 1));
% overflow_pointer = libpointer('int16Ptr', 0);
% 
% % Fetch the data into DLL
% read_flag = calllib('usbtc08', 'usb_tc08_get_single', s, data_pointer, ...
%     overflow_pointer, int16(tc08.units));
% 
% % Check for a successful read
% if read_flag == 0
%     error('Temperature data not read')
% else
%     
% % Fetch all data into matlab including the cold junction
% T = get(data_pointer, 'Value');
% 
% % Discard the cold junction
% T = reshape(T(2:end),1,[]);
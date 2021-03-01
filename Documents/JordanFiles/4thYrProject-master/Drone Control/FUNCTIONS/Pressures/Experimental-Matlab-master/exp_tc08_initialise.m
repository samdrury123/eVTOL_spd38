function s = exp_tc08_initialise(tc08)
% Open a connection to TC-08 thermocouple logger

% Default all channels
if isfield(tc08,'chan') == 0
    tc08.chan = 1:tc08.nchan;
end

% Load library
if ~libisloaded('usbtc08')
    loadlibrary(tc08.library{1}, tc08.library{2})
end

% Open connection with temperature logger
s.h = int16(calllib('usbtc08', 'usb_tc08_open_unit'));

% Set up 8 channels for all the same type of thermocouples
for n = tc08.chan
    calllib('usbtc08', 'usb_tc08_set_channel', s.h, int16(n), int8(tc08.type));
end

% Set period and begin streaming
calllib('usbtc08', 'usb_tc08_run', s.h, int32((n+1)*tc08.period));

% Set units to read
s.units = tc08.units;

% Set number of channels
s.nchan = n;

end
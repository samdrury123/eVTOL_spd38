function s = exp_h309_initialise(h309)
% EXP_H309_INITIALISE  Open a connection to omega thermocouple logger
%
%   Example h309 structure:
%   h309.port = 'COM5'; h309.rate = 9600; h309.term = 'CR/LF';


% Open a serial connection to the data logger
s = exp_serial_open(h309);

% Set defaults

end
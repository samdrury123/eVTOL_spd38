function e = exp_stepper_command(s,n,command,wait_time)
% EXP_STEPPER_COMMAND  Send a command to a parker stepper motor drive unit
%
%   e = EXP_STEPPER_COMMAND(s,n,comm)
%
%   s - input data structure detailing serial connection
%   n - channel number to send command to
%   command - string of command
%   wait_time - optional extended wait time for increased reliability
%   e - output message from drive unit

% Default to short wait time
if exist('wait_time','var') == 0
    wait_time = 0.012;
end

% Clear buffer of data to be read or sent
flushinput(s); flushoutput(s);

% Send command
fprintf(s,[num2str(n) command '\n']);

% Wait for whole packet to arrive
nbytes_old = 0; 
while s.BytesAvailable ~= nbytes_old || s.BytesAvailable == 0
    nbytes_old = s.BytesAvailable;
    pause(wait_time);
end

% Read echo and return function
e = [];
while s.BytesAvailable > 1
    e = [e fgets(s)];
end

end
function [] = exp_serial_wait(t,time_max)
% Wait until data is ready to transmit, if times out quit operation

% Default waiting for 10 seconds
if exist('time_max','var') == 0
    time_max = 10; 
end

% Set increments to wait in
dtime = time_max / 1000;

% While bytes are unavailable pause matlab
time = 0;
while t.BytesAvailable == 0
    pause(time)
    time = time+dtime; 
    if time > time_max; 
        disp('Error - Device not reading'); 
        fclose(t); delete(t); clear t; 
        gtg; 
    end 
end
pause(time)

end
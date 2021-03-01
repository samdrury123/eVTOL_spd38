function T = read_h309(s)
% Log temperatures from omega thermocouple sensor

% Send 'A' command to log everything
flushinput(s); flushoutput(s);
fwrite(s,65);

% Pause for communication
n = 1;
while s.BytesAvailable ~=45
    pause(0.1);
    n = n+1;
    if mod(n,5) == 0
        flushinput(s); flushoutput(s);
        fwrite(s,65);
    end
end

% Read binary data and interpret temperatures
temp = fread(s,7);
T_temp = fread(s,8,'uchar');
temp = fread(s,30);

T = zeros(1,4);
for n = 1:4
    T(n) = (256*T_temp(2*n-1) + T_temp(2*n))/10 + 273.15;
end

end
function [P,T] = exp_dsa_read(s,scan)
% EXP_DSA_READ  Read a single scan from a DSA
%
%   [P,T] = EXP_DSA_READ(s,wait_read)
% 
%   s - input data structure detailing TCP connection or cell array containing multiple connections
%   scan - input string specifying operation - 'scan', 'wait' or 'read'
%   P - output pressure vector in Pascals
%   T - output temperature vector in Kelvin

% Default to scan, wait and read data
if exist('scan','var') == 0
    scan = 'wait';
end

% Check whether multiple or single DSAs are to be read
if iscell(s) == 0
    s = {s};
end

% Initialise temperature and pressure vectors
P = []; T = [];

% Loop over all DSAs to be read and start scan if data isn't already available
for n = 1:length(s)
    if any(strcmp(scan,{'scan','wait'})) == 1
    
        % Clear buffer of data to be read or sent
        flushinput(s{n}); flushoutput(s{n});

        % Send scan command
        fprintf(s{n}, 'SCAN\n'); 
        
    end
end

% Read scan on each DSA as it comes back
for n = 1:length(s)
    if any(strcmp(scan,{'read','wait'})) == 1
    
        % Pause until whole packet is available
        while s{n}.BytesAvailable ~= 104
            pause(0.01);    
        end

        % Check packet type
        a = fread(s{n},2,'uint8');
        packet_type = hex2dec([dec2hex(a(2)) dec2hex(a(1))]);

        % Error if the packet type is incorrect
        if packet_type ~= 5
            error('Packet Type not Scan EU')
        end

        % Discard two bytes
        fread(s{n},1,'int16');

        % Frame number
        fread(s{n},4,'uint8');

        % Read pressures
        a = fread(s{n},64,'uint8');
        P = [P reshape(hexsingle2num([dec2hex(a(4:4:64)) dec2hex(a(3:4:63)) ...
            dec2hex(a(2:4:62)) dec2hex(a(1:4:61))]),1,[])];

        % Read temperatures
        a = fread(s{n},32,'uint8');
        T = [T reshape(hex2dec([dec2hex(a(2:2:32)) dec2hex(a(1:2:31))]),1,[])];
        
    end
end

end
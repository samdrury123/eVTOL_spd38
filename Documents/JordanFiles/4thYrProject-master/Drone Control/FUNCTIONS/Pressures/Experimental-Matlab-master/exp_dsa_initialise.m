function s = exp_dsa_initialise(dsa)
% Initialise a DSA with desired configuration variables and calibrate zero

% Check whether multiple or single DSAs are to be initialised
if iscell(dsa) == 0
    dsa = {dsa};
end

% Initialise connection cell array
s = cell(length(dsa),1);

% Loop over all DSAs
for n = 1:length(dsa)

    % Create TCP/IP object and connect to the server
    s{n} = tcpip(dsa{n}.ip, dsa{n}.port); 
    s{n}.InputBufferSize = 30000; 
    fopen(s{n});

    % Set configuration variables
    % fprintf(s{n}, 'SET BIN 0\n');                              % ASCII output
    fprintf(s{n}, 'SET BIN 1\n');                              % Binary output
    fprintf(s{n}, ['SET AVG ' num2str(dsa{n}.avg) '\n']);         % Number of samples per channel
    fprintf(s{n}, ['SET PERIOD ' num2str(dsa{n}.period) '\n']);   % Period between samples
    fprintf(s{n}, 'SET FPS 1\n');                              % 1 frame only
    fprintf(s{n}, 'SET EU 1\n');                               % Engineering units
    % fprintf(s{n}, 'SET EU 0\n');                               % Raw P and T counts
    fprintf(s{n}, 'SET UNITSCAN PA\n');                        % Return pascals

    % Save configuration
    fprintf(s{n}, 'SAVE\n');

    % Long calibrate 0
    % fprintf(s, 'CALZ 240 65000 15\n');
    % pause(15);

    % Short calibrate 0
    fprintf(s{n}, 'CALZ\n');
    pause(10);
end

% If only one DSA is used save in a struct
if length(s) == 1
    s = s{1};
end

end
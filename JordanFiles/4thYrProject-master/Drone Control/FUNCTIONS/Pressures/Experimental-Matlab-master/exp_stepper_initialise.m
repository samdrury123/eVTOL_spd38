function s = exp_stepper_initialise(park,set_address)
% EXP_STEPPER_INITIALISE  Open a serial connection to Parker drive unit
%
%   s = EXP_STEPPER_INITIALISE(park)
%
%   park - input data structure of connection definition
%   set_address - 0 or 1 to set drive addressing for the first time
%   s - output data structure detailing serial connection
%
%   Example park data structure for connection:
%       park.port = 'COM1'; park.rate = 9600; 
%       park.term = 'CR'; park.nchan = 4;
%       park.current = [1.8 0.6 0.6]; 
%       park.resistance = [0.72 6.2 4.2]; 
%       park.inductance = [1 8.8 2.5];

% Default to not set address
if exist('set_address','var') == 0
    set_address = 0;
end

% Choose drive type
drive_type = 'VIX';
% drive_type = 'DG';

% Open connection
s = exp_serial_open(park);

% Set addressing
if set_address == 1
    
    % Set automatic addressing of daisy-chained boxes
    fprintf(s,'#1\n'); pause(5); fprintf(s,'0SV\n'); pause(5); fprintf(s,'0Z\n');
    
    % Explicitly set addressing of all axes
    for n = 1:park.nchan
        fprintf(s,[num2str(n) '#' num2str(n) '\n']);
    end
    pause(5); fprintf(s,'0SV\n'); pause(5); fprintf(s,'0Z\n'); pause(5);

    % Save and reset
    pause(5);
    exp_stepper_command(s,0,'SV'); pause(5);
    exp_stepper_command(s,0,'Z'); pause(5);
end

% Enable and check motor currents
current_set = 1;
for n = 1:park.nchan
    
    % Enable motors
    exp_stepper_command(s,n,'E1');
    
    % Query motor settings
    motor = exp_stepper_command(s,n,'MOTOR',0.03);

    % Check currents
    mc = regexp(motor,'MC(\w+.\w+)','tokens'); mc = str2double(mc{1});
    if mc ~= park.current(n)
        current_set = 0;
    end
end

% Apply motor settings
if current_set == 0
    
    % Loop over each axis and set motors
    for n = 1:park.nchan
        
        % Enable motors
        exp_stepper_command(s,n,'E1');

        % Configure motor types - 800 steps/rev & series winding
        exp_stepper_command(s,n,['MOTOR(513,' num2str(park.current(n)) ',800,720,4,'...
            num2str(park.resistance(n)) ',' num2str(park.inductance(n)) ')']);
        pause(5);
    end

    % Save and reset
    pause(5);
    exp_stepper_command(s,0,'SV'); pause(5);
    exp_stepper_command(s,0,'Z'); pause(5);
end

% Initialise motors
for n = 1:park.nchan
    
    % Power motors up
    exp_stepper_command(s,n,'ON');

    % Set zeros and enable absolute position referencing
    exp_stepper_command(s,n,'W(PA,0)');
    exp_stepper_command(s,n,'MA');

    % Set velocities and accelerations
    exp_stepper_command(s,n,['V' num2str(park.V(n))]);
    exp_stepper_command(s,n,['A' num2str(park.A(n))]);

    % Disable limits
    exp_stepper_command(s,n,'LIMITS(3,0,0)');
end

% Enable brake
if park.nchan == 3
    exp_stepper_command(s,3,'O(X1X)');
end

% Setup digiplan drive unit
if strcmp(drive_type,'DG') == 1
    for n = 1:park.nchan
        exp_stepper_command(s,n,'E');
        exp_stepper_command(s,n,'LD3');
        exp_stepper_command(s,n,'MPA');
        exp_stepper_command(s,n,'MN');
        exp_stepper_command(s,n,'PZ');
        exp_stepper_command(s,n,'V5'); % 5 rev s^(-1)
        exp_stepper_command(s,n,'A10'); % 10 rev s^(-2)
    end
end


end
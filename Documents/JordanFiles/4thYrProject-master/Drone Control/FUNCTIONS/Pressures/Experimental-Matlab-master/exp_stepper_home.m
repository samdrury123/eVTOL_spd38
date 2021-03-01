function [] = exp_stepper_home(s,home,drive_all)
% EXP_STEPPER_HOME  Move stepper motors to home positions
%
%   [] = EXP_STEPPER_HOME(s,home,drive_all,V_h)
%
%   s - input data structure detailing serial connection
%   home - homing data structure to define directions and speeds
%   drive_all - 0 or 1 to drive motors simultaneously
%
%   home has fields:
%       e - edge to stop on '-' or '+'
%       t - switching type: normally open - 0 or normally closed - 1 
%       V - homing velocity and direction in rev s-1

% Set default to drive one motor at a time
if exist('drive_all','var') == 0
    drive_all = 0;
end

% Set final homing velocity
for n = 1:length(home.t)
    exp_stepper_command(s,n,['W(HF,' num2str(home.Vh(n)) ')']);
end

% Choose whether to drive all motors at once or one at a time
if drive_all == 1
    
    % Send home command on all channels
    for n = 1:length(home.t)
        start_home(s,n,home.e(n),home.t(n),home.V(n));
    end
    
    % Check in home position
    for n = 1:length(home.t)
        is_home(s,n);
    end
    
else
    
    % Loop over all channels, start homing and wait until finished
    for n = 1:length(home.t)
        start_home(s,n,home.e(n),home.t(n),home.V(n));
        is_home(s,n);
    end
end

end

function [] = start_home(s,n,e,t,V)
% Start home maneuvre

% Arm command
exp_stepper_command(s,n,'ARM1');

% Configure home options
exp_stepper_command(s,n,['HOME1(' e ',' num2str(t) ',' num2str(V) ,',10,2)']);

% Execute move
exp_stepper_command(s,n,'GH');

end

function [] = is_home(s,n)
% Check if motor is home yet

% Query 7th status bit continuously until no longer homing
homing = 1;
while homing == 1
    pause(0.1)
    status = regexp(exp_stepper_command(s,n,'R(ST)'),'\n','split');
    homing = str2double(status{end}(7));
end

end
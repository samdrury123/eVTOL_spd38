function [] = exp_stepper_coord(s,q,limit,drive_all)
% EXP_STEPPER_COORD  Move all stepper motors to desired coordinate vector
%
%   [] = EXP_STEPPER_COORD(s,q,limit,drive_all)
%
%   s - input data structure detailing serial connection
%   q - coordinate vector in raw step numbers
%   limit - compulsory limit array to avoid overstepping
%   drive_all - 0 or 1 to drive motors simultaneously
%
%   Example limit structure, one line for each channel:
%       limit = [-18799 23921 ; -35000 150000 ; -180000 180000];

% Set default to drive one motor at a time
if exist('drive_all','var') == 0
    drive_all = 0;
end

% Check limits
for n = 1:length(q)
    if q(n) < limit(n,1) || q(n) > limit(n,2)
        error(['Channel ' num2str(n) ' limit exceeded']);
    end
end

% Choose whether to drive all motors at once or one at a time
if drive_all == 1
    
    % Loop over all channels and start movements
    for n = 1:length(q)
        exp_stepper_move(s,n,q(n),limit)
    end
    
    % Loop over all channels and query positions
    for n = 1:length(q)
        exp_stepper_position(s,n,q(n))
    end

else 
    
    % Loop over all channels and drive individually
    for n = 1:length(q)
        
        % Begin movement
        exp_stepper_move(s,n,q(n),limit)
        
        % Check if in position
        exp_stepper_position(s,n,q(n))
    end
    
end

end
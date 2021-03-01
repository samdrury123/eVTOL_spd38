function [] = exp_stepper_position(s,n,q)
% EXP_STEPPER_POSITION  Check coordinates of stepper motors
%
%   [] = EXP_STEPPER_POSITION(s,n,q)
%
%   s - input data structure detailing serial connection
%   n - channel number of motor to check position
%   q - coordinate vector in raw step numbers

% Choose position reporting string
% pos_type = 'PA'; % Parker Paragon / L50i
pos_type = 'PT'; % Parker VIX500
% pos_type = 'DG'; % Digiplan

% Initialise in position flag
in_pos = 0;

% Query position until flag is high
while in_pos == 0
        
    % Send position request
    pos = exp_stepper_command(s,n,['R(' pos_type ')']);

    % Extract position number from string, pause if failed
    if isempty(pos) == 0
        if strcmp(pos_type,'DG') == 0
            pos = regexp(pos,'*','split');
            if length(pos) >= 2
                pos = str2double(pos{2});
            end
        else
            if length(pos) > 12
                pos = str2double(pos(end-11:end-1));
            else
                pos = nan;
            end
        end
    end
    
    % Check position is a single value
    if numel(pos) ~= 1
        pos = nan;
    end

    % Compare to step coordinate
    if iscell(pos) == 0 && pos == round(q)
        in_pos = 1;
    end
    
end

end
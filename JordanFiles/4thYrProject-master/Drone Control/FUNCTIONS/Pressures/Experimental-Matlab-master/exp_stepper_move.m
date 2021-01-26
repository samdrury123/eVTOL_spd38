function [] = exp_stepper_move(s,n,q,limit)
% EXP_STEPPER_MOVE  Command a single motor to move to a new position
%
%   [] = EXP_STEPPER_MOVE(s,n,q)
%
%   s - input data structure detailing serial connection
%   n - channel number of motor to move
%   q - coordinate in raw step number

% Check limit
if q < limit(n,1) || q > limit(n,2)
    error(['Channel ' num2str(n) ' limit exceeded']);
end

% Program new absolute step position
fprintf(s,[num2str(n) 'D' num2str(round(q)) '\n']);

% Send move command
fprintf(s,[num2str(n) 'G\n']);

end
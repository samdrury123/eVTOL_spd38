function [q,t] = exp_pxie_read(s)
% EXP_PXIE_READ  Read already set up channels from PXIE ADC
%
%   [q,t] = exp_pxie_read(s)
%
%   s - connection data structure
%   q - raw voltage matrix
%   t - time vector in seconds

% Excecute samples
[q,t] = s.startForeground();

end
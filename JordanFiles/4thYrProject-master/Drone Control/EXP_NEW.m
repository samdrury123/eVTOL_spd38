function [S] = EXP_NEW(type, intake, exit, sigma)
% EXP_NEW generate new data set with correct structure fields
%  

S.FOM = [];
S.P = [];
S.T = [];
S.mass = [];
S.rpmmean = [];
S.nomV = [];
S.pwm = [];
S.sig = sigma;
S.intake = intake;
S.exit = exit;
S.type = type;

end


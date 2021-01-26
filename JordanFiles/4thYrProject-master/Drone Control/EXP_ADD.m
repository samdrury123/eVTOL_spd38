function S = EXP_ADD(S, data)
% Add data to end of dataset
%

S.FOM = [S.FOM data.FOM];
S.P = [S.P data.P];
S.T = [S.T data.T];
S.mass = [S.mass data.mass];
S.rpmmean = [S.rpmmean data.rpmmean];
S.nomV = [S.nomV data.nomV];
S.pwm = [S.pwm data.pwm];

end

% Use DSA-343 (eVTOL DSA)
dsa.ip = '169.254.96.89';

% Set default DSA logging variables
dsa.port = 22; dsa.nchan = 16;
dsa.period = 150; dsa.avg = 150;

fan_dsa = exp_dsa_initialise(dsa);

%%
P_chan = dsa.nchan;

for i = 1:10
    e.P(i,:) = exp_dsa_read(fan_dsa);
end
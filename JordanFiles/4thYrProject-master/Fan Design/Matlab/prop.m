function [rpm] = proprpm();
% Return required rpm for prop thrust = weight
% 

load('vars', 'W', 'rpm');
%% Determine weight of drone

weight = W.body + W.motor;

testrpm = linspace(2000,16000,1000);

%% Load APC data
load('APC_FIT.mat', 'APC_FIT');

testthrust = APC_FIT.rpm_thr(testrpm);

for i = 1:length(testrpm)
    if testthrust(i) > weight
        idx = i;
        rpm = testrpm(idx);
    end
end

end


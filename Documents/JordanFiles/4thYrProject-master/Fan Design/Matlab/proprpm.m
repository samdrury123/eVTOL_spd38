function [rpm, weight, power] = proprpm()
% Return required rpm for prop thrust = weight and weight of prop design
% 
load('efficiencies.mat');
load('vars', 'W', 'rpm');
%% Determine weight of drone

weight = W.body + W.prop;

testrpm = linspace(2000,16000,1000);

%% Load APC data
load('APC_FIT.mat', 'APC_FIT');

testthrust = APC_FIT.rpm_thr(testrpm);

for i = 1:1000
    if testthrust(i) > weight
        idx = i;
        rpm = testrpm(idx);
        break
    end
end

power = APC_FIT.rpm_pwr(rpm)/(eta.control * eta.motor);

end


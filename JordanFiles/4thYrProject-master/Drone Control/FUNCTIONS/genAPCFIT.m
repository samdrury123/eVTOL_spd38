% Import APC data for 10x3.8SF prop
% Fit curves for RPM PWR THRUST

% HP to W
pwrconv = 745.7;

% LBF to N
thrconv = 4.44822;

% Import
[RPMAPC, V5APC, JAPC, PeAPC, CtAPC, CpAPC, PWRAPC, TorqueAPC, THRAPC] = importAPC("/Users/jordaneriksen/Documents/Uni/Part 2B/Project/Drone Control/APCPROPDATA.txt", [4, 18]);

PWRAPC = pwrconv*PWRAPC;
THRAPC = thrconv*THRAPC;

APC_FIT.rpm_pwr = APC_RPM_PWR(RPMAPC, PWRAPC);
APC_FIT.rpm_thr = APC_RPM_THR(RPMAPC, THRAPC);

save('APC_FIT', 'APC_FIT');

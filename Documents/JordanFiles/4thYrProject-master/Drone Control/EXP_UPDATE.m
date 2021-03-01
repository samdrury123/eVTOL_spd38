%% Update Data 

load('EXP_DATA.mat','S08S','S10S','S12S','S08L','S10L','S12L','S08','S10','S12','R10','R17','R20');
load('EXP_DATA_PROP.mat', 'APC', 'BASE');

load('EXP_META.mat')
data = load('EXP_OP.mat');
data.nomV = nomV;
data.pwm = pwm;

switch SIGMANAME
    case '8'
        S08 = EXP_ADD(S08, data);
    case '08S'
        S08S = EXP_ADD(S08S, data);
    case '08L'
        S08L = EXP_ADD(S08L, data);
    case '10'
        S10 = EXP_ADD(S10, data);
    case '10S'
        S10S = EXP_ADD(S10S, data);
    case '10L'
        S10L = EXP_ADD(S10L, data);
    case '12'
        S12 = EXP_ADD(S12, data);
    case '12S'
        S12S = EXP_ADD(S12S, data);
    case '12L'
        S12L = EXP_ADD(S12L, data);
    case 'R10'
        R10.FOM = [R10.FOM data.FOM];
        R10.P = [R10.P data.P];
        R10.T = [R10.T data.T];
        R10.mass = [R10.mass data.mass];
        R10.rpmmean = [R10.rpmmean data.rpmmean];
        R10.pwm = [R10.pwm data.pwm];
        R10.nomV = [R10.nomV data.nomV];
        R10.sig = [R10.sig sigma];
        R10.intake = [R10.intake, intake];
        R10.exit = [R10.exit, exit];
    case 'R17'
        R17.FOM = [R17.FOM data.FOM];
        R17.P = [R17.P data.P];
        R17.T = [R17.T data.T];
        R17.mass = [R17.mass data.mass];
        R17.rpmmean = [R17.rpmmean data.rpmmean];
        R17.pwm = [R17.pwm data.pwm];
        R17.nomV = [R17.nomV data.nomV];
        R17.sig = [R17.sig sigma];
        R17.intake = [R17.intake, intake];
        R17.exit = [R17.exit, exit];
    case 'R20'
        R20.FOM = [R20.FOM data.FOM];
        R20.P = [R20.P data.P];
        R20.T = [R20.T data.T];
        R20.mass = [R20.mass data.mass];
        R20.rpmmean = [R20.rpmmean data.rpmmean];
        R20.pwm = [R20.pwm data.pwm];
        R20.nomV = [R20.nomV data.nomV];
        R20.sig = [R20.sig sigma];
        R20.intake = [R20.intake, intake];
        R20.exit = [R20.exit, exit];
    case 'APC'
        PROP.FOM = [PROP.FOM data.FOM];
        PROP.P = [PROP.P data.P];
        PROP.T = [PROP.T data.T];
        PROP.mass = [PROP.mass data.mass];
        PROP.pwm = [PROP.pwm data.pwm];
        PROP.nomV = [PROP.nomV data.nomV];
    case 'BASE'
        BASE.FOM = [BASE.FOM data.FOM];
        BASE.P = [BASE.P data.P];
        BASE.T = [BASE.T data.T];
        BASE.mass = [BASE.mass data.mass];
        BASE.pwm = [BASE.pwm data.pwm];
        BASE.nomV = [BASE.nomV data.nomV];
end

save('EXP_DATA.mat','S08S','S10S','S12S','S08L','S10L','S12L','S08','S10','S12','R10','R17','R20');
save('EXP_DATA_PROP.mat', 'APC', 'BASE');

clear all;

% EXP_PLOT;
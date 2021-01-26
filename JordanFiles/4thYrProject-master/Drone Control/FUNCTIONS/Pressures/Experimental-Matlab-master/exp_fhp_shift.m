function P = exp_fhp_shift(P,r,rt,probe,plot_stuff)
% EXP_FHP_SHIFT  Correct pressures by interpolating onto the centre hole location

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Plot pressure values before for comparison
if plot_stuff == 1
    figure(); cols = lines(5);
    subplot(1,2,1); hold on; grid on; box on; xlabel('Pitch'); ylabel('Pressure');
    j_mid = round(size(P,1) / 2); for n = 1:5; plot(rt(j_mid,:),P(j_mid,:,n),'color',cols(n,:)); end;
    subplot(1,2,2); hold on; grid on; box on; xlabel('Pressure'); ylabel('Radius');
    k_mid = round(size(P,2) / 2); for n = 1:5; plot(P(:,k_mid,n),r(:,k_mid),'color',cols(n,:)); end;
end

% Hole offset value
if strcmp(probe.arrange,'cross') == 1
    d = probe.D / (3 * sqrt(2));
else
    d = probe.D / 3;
end

% Channel names for correction
N = probe.N;
if strcmp(probe.arrange,'cross') == 1
    N_up = {'P_ur' 'P_ul'}; N_dn = {'P_dr' 'P_dl'};
    N_rt = {'P_ur' 'P_dr'}; N_lf = {'P_ul' 'P_dl'};
else
    N_up = {'P_up'}; N_dn = {'P_dn'}; N_rt = {'P_rt'}; N_lf = {'P_l'};
end

% Interpolate in radial direction
for k = 1:size(P,2)
    for v = 1:length(N_up)
        P(:,k,N.(N_up{v})) = interp1(r(:,k),P(:,k,N.(N_up{v})),r(:,k) - d,'pchip','extrap');
    end
    for v = 1:length(N_dn)
        P(:,k,N.(N_dn{v})) = interp1(r(:,k),P(:,k,N.(N_dn{v})),r(:,k) + d,'pchip','extrap');
    end    
end

% Interpolate in pitchwise direction
for j = 1:size(P,1)
    for v = 1:length(N_lf)
        P(j,:,N.(N_lf{v})) = interp1(rt(j,:),P(j,:,N.(N_lf{v})),rt(j,:) + d,'pchip','extrap');
    end
    for v = 1:length(N_rt)
        P(j,:,N.(N_rt{v})) = interp1(rt(j,:),P(j,:,N.(N_rt{v})),rt(j,:) - d,'pchip','extrap');
    end    
end

% Plot pressure values after correction
if plot_stuff == 1
    subplot(1,2,1); for n = 1:5; plot(rt(j_mid,:),P(j_mid,:,n),'--','color',cols(n,:)); end;
    subplot(1,2,2); for n = 1:5; plot(P(:,k_mid,n),r(:,k_mid),'--','color',cols(n,:)); end;
end


end
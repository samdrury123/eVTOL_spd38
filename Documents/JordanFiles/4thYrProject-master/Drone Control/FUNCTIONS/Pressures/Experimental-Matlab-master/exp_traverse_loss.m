function [l,G,h,E] = exp_traverse_loss(directory,trav_names,calib_name,N,r,x,des,h,C,plot_stuff)
% Plot interesting results from a set of Gibbons rig traverses

% Default to plot figures
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

% Default colour to black
if exist('C','var') == 0 || isempty(C) == 1
    C = [0 0 0];
end

% Open figure windows
if exist('h','var') == 0 || isempty(h) == 1
    h.chic = figure(); hold on; grid on; box on;
    xlabel('Flow Coefficient'); ylabel('Pressure Rise Coefficient');
    
    h.flow = figure(); hold on; grid on; box on;
    xlabel('Point'); ylabel('Flow / Pressure Rise Coefficients');
        
    h.mid = figure(); hold on; grid on; box on;
    xlabel('Tangential Position'); ylabel('Stagnation Pressure');
    
    h.foils = figure(); 
    subplot(1,3,1); hold on; grid on; box on;
    xlabel('Flow Coefficient'); ylabel('Radius');
    subplot(1,3,2); hold on; grid on; box on;
    xlabel('Whirl Angle'); ylabel('Radius');
    subplot(1,3,3); hold on; grid on; box on;
    xlabel('Pitch Angle'); ylabel('Radius');
    
    h.loss = figure(); hold on; grid on; box on;
    xlabel('Stagnation Pressure Loss'); ylabel('Radius');
    axis([0 0.2 0 1]);
    
    h.Po = figure(); hold on; grid on; box on;
    xlabel('Stagnation Pressure'); ylabel('Radius');
    
    h.error = figure(); hold on; grid on; box on;
    xlabel('FHP Error / ^{\circ}'); ylabel('Radius');
end

% Load in all raw data
E = cell(length(trav_names),1); G = cell(length(trav_names),1);
for n = 1:length(trav_names)
    load([directory trav_names{n}]);
    E{n} = e;
end
    
% Load in calibration maps
if isempty(calib_name) == 0
    c_temp = cell(length(calib_name));
    for n = 1:length(calib_name);
        load([directory calib_name{n}]);
        c_temp{n} = c;
    end
    c = c_temp;
else
    c = [];
end

%% Plot characteristic and traverse operating points

% Process all traverse data and plot operating points
M = {'o' 'x' '^' 'd' 'p' 's'};
if length(E) == 10
    n_trav = [1 4];
    n_chic = 7:8;
elseif length(E) == 3
    n_trav = 1:2; n_chic = 3;
end
for n = n_trav
    G{n} = exp_traverse_calculate(E{n},c,x(n),r,N,des);
    if plot_stuff == 1
        figure(h.chic);
        plot(mean(G{n}.mrtp),mean(G{n}.Rp),'-','Color',C,'Marker',M{n},'MarkerFaceColor',C);
    end
end

% Calculate static-static characteristics and plot on top of traverses
N.Po_in = N.rig.Po; N.P_in = N.rig.P; N.To_in = N.rig.T; N.P_out = N.rig.P_out;
for n = n_chic
    G{n} = exp_chic_statics(directory,trav_names{n},N,r.mid,[],[],0);
    if plot_stuff == 1
        figure(h.chic);
        plot(G{n}.mrtp,G{n}.Rp,'-','Color',C,'LineWidth',2);
    end
end

%% Plot five hole probe error calculations
if plot_stuff == 1
    figure(h.error)
    for n = n_trav
        r_nondim = ts_mass_average(G{n},'r_nondim',2);
        A = ts_mass_average(G{n},'A',2);
        plot(A,r_nondim,'-','Color',C);
    end
end

%% Plot midspan total pressures to check non-dimensionalisation
if plot_stuff == 1
    figure(h.mid);
    for n = n_trav
        j = round(size(G{n}.t,1)/2);
        plot(G{n}.t(j,:),G{n}.Po(j,:),'-','Color',C,'Marker',M{n});
        plot(G{n}.t(j,:),G{n}.P_corr(j,:,N.P_out),'--','Color',C,'Marker',M{n});
        plot(G{n}.t(j,:),G{n}.P_corr(j,:,N.Po_in),'--','Color',C,'Marker',M{n});
    end
end

% Plot flow and pressure rise coefficient traces for rig drift
if plot_stuff == 1
    figure(h.flow)
    for n = n_trav
        plot(G{n}.phi,'-','Color',C);
        plot(G{n}.psi,'-','Color',C);
    end
end

%% Plot circumferentially averaged foils for inlet and exit

% Loop over all traverses
for n = n_trav
    
    % Calculate tangentially averaged properties
    r_nondim = ts_mass_average(G{n},'r_nondim',2);
    phi = ts_mass_average(G{n},'Vx',2) / G{n}.U_mid;
    Alpha = ts_mass_average(G{n},'Alpha',2);
    Beta = ts_mass_average(G{n},'Beta',2);
    
    if plot_stuff == 1
        figure(h.foils);
        subplot(1,3,1); plot(phi,r_nondim,'-','Color',C);
        subplot(1,3,2); plot(Alpha,r_nondim,'-','Color',C);
        subplot(1,3,3); plot(Beta,r_nondim,'-','Color',C);
    end
end


%% Plot radial loss coefficients

if plot_stuff == 1
    % Select inlet and outlet traverses
    n_out = n_trav(1); n_in = n_trav(2);
    
    % Loop over both passages
%     for n = 1:2
%         g = G{n_out};
%         varnames = fieldnames(g);
%         for v = 1:length(varnames)
%             if isa(g.(varnames{v}),'single') == 1
%                 if n == 1
%                     g.(varnames{v}) = g.(varnames{v})(:,1:27);
%                 else
%                     g.(varnames{v}) = g.(varnames{v})(:,28:end);
%                 end
%             end
%         end
    g = G{n_out};
    
    % Average pressures and radii
    [Po_in, mass_in] = ts_mass_average(G{n_in},'Po',2);
    r_in = ts_mass_average(G{n_in},'r_nondim',2);
    [Po_out, mass_out] = ts_mass_average(g,'Po',2);
    r_out = ts_mass_average(g,'r_nondim',2);
    dP = ts_mass_average(G{n_in},'Po',3) - ts_mass_average(G{n_in},'P',3);
    
    % Calculate mass fractions
    mass_in = [0 cumsum(mass_in)]; mass_out = [0 cumsum(mass_out)];
    mass_in = mass_in / max(mass_in); mass_out = mass_out / max(mass_out);
    
    % Radial loss coefficient based on inlet average
%     Yp = (ts_mass_average(G{n_in},'Po',3) - Po_out) / dP;
    
    % Radial loss coefficient based on cumilative mass fraction
    Yp = (interp1(mass_in,Po_in,mass_out,'pchip') - Po_out) / dP;
    
    % Radial loss coefficient based on finite mass fractions
%     m_frac = 0.1; Po_in_frac = zeros(length(Po_in),1);
%     for j = 1:length(Po_out)
%         m_temp = linspace(max(mass_in(j)-m_frac,0),min(mass_in(j)+m_frac,1),10);
%         Po_in_frac(j) = mean(interp1(mass_in,Po_in,m_temp,'pchip'));
%     end
%     Yp = (interp1(mass_in,Po_in_frac,mass_out,'pchip') - Po_out) / dP;
%     
    
    % Total loss
    l.Yp_total = (ts_mass_average(G{n_in},'Po',3) - ts_mass_average(g,'Po',3)) / dP;

    % Correct for rotor wake mixing loss
%     Yp = Yp - 250/dP;

    l.Yp = Yp; l.r = r_out;
    
    % Calculate exit based free stream total pressure
    Po_max = sort(g.Po,2); Po_max = trimmean(Po_max(:,end-10:end),20,'round',2);
    
    % Plot loss profile
    figure(h.loss);
    plot(Yp,r_out,'-','Color',C);
    v = axis; axis([v(1:2) 0 1]);

    % Plot total pressure profiles
    figure(h.Po);
    plot(Po_in,r_in,'-','Color',C);
    plot(Po_out,r_out,'-','Color',C);
    plot(Po_max,r_out,'-','Color',C);
    
    % Record freestream total pressures
    G{n_out}.Po_max = repmat(Po_max,[1 size(G{n_out}.Po,2)]);
%     end
end



%% Contour plot of loss

if plot_stuff == 1
    % Select inlet and outlet traverses
    n_out = n_trav(1); n_in = n_trav(2);
    
    h.wake = figure(); hold on; axis equal;
        
    % Mass average coefficients at inlet
    Po_in = ts_mass_average(G{n_in},'Po',3);
    P_in = ts_mass_average(G{n_in},'P',3);

    % Calculate loss coefficient at exit
    Yp = (Po_in - G{n_out}.Po) ./ (Po_in - P_in);
    
    % Plot contours
    Yp_min = -0.05; Yp_max = 0.4;
    contourf(G{n_out}.z,G{n_out}.y,Yp,linspace(Yp_min,Yp_max,12));
    colormap(flipud(bone)); colorbar; caxis([Yp_min Yp_max]);

    % Plot hub and casing lines
    [z_hub,y_hub] = pol2cart(linspace(G{n_out}.t(1,1),G{n_out}.t(1,end),100)+pi/2,...
        ones(1,100)*r.hub);
    [z_cas,y_cas] = pol2cart(linspace(G{n_out}.t(end,1),G{n_out}.t(end,end),100)+pi/2,...
        ones(1,100)*r.cas);

    plot(-z_hub,y_hub,'k--'); plot(-z_cas,y_cas,'k--');
end

%% Calculate endwall and profile loss coefficients

% p.Yp = ts_break_loss_mass_2(G{1},G{2},0.25);
% p = [];

%% Calculate mixed out loss coefficient
% figure(h.loss);
% 
% % Mix out flow condition
% c_mix = ts_mix_out_average(G{2});
% Po_mix = c_mix.Po_2;
% 
% % Turn mixed out stagnation pressure into a coefficient
% P_ref = mean(G{2}.P_corr(:,:,N_ref.P),2);
% T_ref = mean(G{2}.T_corr(:,:,N_ref.T),2);
% ro_ref = P_ref ./ (air.R * T_ref);
% 
% Cpo_mix = (Po_mix - P_ref) ./ (ro_ref * U_mid^2);
% 
% % Based on mass
% [Po_in, mass_in] = ts_mass_average(G{1},'Cpo',2);
% 
% % Calculate mass flow fractions at inlet and exit
% [~, mass_out] = ts_mass_average(G{2},'Cpo',2);
% mass_in = [0 cumsum(mass_in)]; mass_out = [0 cumsum(mass_out)];
% mass_in = mass_in / max(mass_in); mass_out = mass_out / max(mass_out);
% 
% % Calculate entropy loss coefficient
% % p.Yp_mix = (max(max(G{2}.Cpo)) - Cpo_mix.') ./...
% %     (ts_mass_average(G{1},'Cpo',3) - ts_mass_average(G{1},'Cp',3));
% l.Yp_mix = (interp1(mass_in,Po_in,mass_out,'pchip') - Cpo_mix.') ./...
%     (ts_mass_average(G{1},'Cpo',3) - ts_mass_average(G{1},'Cp',3));
% 
% plot(l.Yp_mix,r_out,'--','Color',C);
% 
% % Calculate total mixed out loss coefficient
% [~, mass_out] = ts_mass_average(G{2},'Cpo',2);
% l.Yp_mix_total = sum(0.5*(l.Yp_mix(1:end-1) + l.Yp_mix(2:end)) .* mass_out) / sum(mass_out);


%% Calculate incidence
l.Alpha = interp1(ts_mass_average(G{n_in},'r_nondim',2),ts_mass_average(G{n_in},'Alpha',2),0.5,'pchip');

end


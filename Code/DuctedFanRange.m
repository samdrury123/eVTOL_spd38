%   DuctedFanRange
%                               
%   Investigate operating range for different diameter fans
%
%   Written by James Taylor & Sam Grimshaw               
%   February 2021

% Clear workspace and load eVTOL environment
clear; close all; clc; eVTOL;

% Script operation
plot_diameter = 1;


%% Investigate fan diameter on operating range for different designs
if plot_diameter == 1
   
    % Operating points
    Vf = [55.4 0]; F = [1200 6000]; alt = [5000 0]; 
    
    % Design point at cruise
    phi_des = 0.7; sig = 1; r = linspace(0.5,0.9,21); lam = 0.4; eta = 0.85;
    
    % Air conditions at altitude
    for n = 1:length(Vf)
        [Ta(n),Pa(n),ro(n),~] = Altitude(alt(n));
    end
    
    % Figure window for chics
    h.chic = figure(); hold on; grid on; box on; ns = [1 11 21]; cols = lines(length(ns));
    xlabel('Flow Coefficient'); ylabel('Stage Loading Coefficient'); axis([0 1 0 1]);
    
    % Loop over all radii and calculate required fan operating range 
    p = zeros(size(r)); psi_des = zeros(size(r)); psi_sta = zeros(size(r)); 
    for n = 1:length(r)
        
        % Fan disc area
        A = pi * r(n)^2 * (1 - lam^2);
        
        % Calculate axial velocity from thrust quadratic formula 
        Vx = (ro .* A .* Vf + (ro.^2 .* A.^2 .* Vf.^2 + 4 * ro .* A .* F / sig).^0.5) ./ ...
            (2 * ro .* A / sig);
        phi(1) = phi_des; Vj = Vx / sig;
        
        % Blade speed at cruise point - note double as using relative blade speed
        U(1) = Vx(1) / phi(1);
        Cptt(1) = (Vj(1)^2 - Vf(1)^2) / U(1)^2;
        
        % Compressor inlet and outlet conditions
        Po_1 = Pa + 0.5 * ro .* Vf.^2; Po_2 = Pa + 0.5 * ro .* Vj.^2;
        To_1 = Ta; To_2 = To_1 .* (Po_2 ./ Po_1).^((ga-1)/(ga * eta));
        
        % Stage loading coefficient at cruise
        psi(1) = cp * (To_2(1) - To_1(1)) / U(1)^2;
        
        % Calculate blade speed at take-off
%         U_list = linspace(0.5*U(1),2*U(1),100);
%         U_list = linspace(0.1,1.2,100) * Vx(1) / phi(1);
        U_list = [linspace(0.5 * Vx(1) / phi(1),2*U(1),100) 1e5];
        psi_list = cp * (To_2(2) - To_1(2)) ./ U_list.^2;
        phi_list = Vx(2) ./ U_list;
        psi_chic = (psi(1) - 0.5) * phi_list / phi(1) + 0.5;
        U(2) = interp1(psi_chic - psi_list,U_list,0,'spline');
        
        % Calculate flow and stage loading coefficients at take-off
        phi(2) = Vx(2) / U(2);
        psi(2) = cp * (To_2(2) - To_1(2)) / U(2)^2;
        
        % Record flow range requirement and stage loading coefficients
        p(n) = (phi(1) - phi(2)) / phi(1);
        
        % Record stage loading coefficients
        psi_des(n) = psi(1); psi_sta(n) = psi(2);
        
        % Plot the chic and operating points for selected cases
        if any(n == ns) == 1
            figure(h.chic); l(n == ns) = plot(phi_list,psi_chic,'-','color',cols(n == ns,:));
            plot(phi(1),psi(1),'s','color',cols(n == ns,:))
            plot(phi(2),psi(2),'x','color',cols(n == ns,:))
            leg{n == ns} = ['D = ' num2str(r(n) * 2)];
        end

    end
    
    % Add legend to chics
    figure(h.chic); legend(l,leg); 
    
    % Plot flow range
    figure(); hold on; grid on; box on; 
    xlabel('Fan Diameter / m'); ylabel('Required Flow Range Reduction');
    plot(r*2,p,'-');
%     export_pres(gcf,[dr.fig 'ideal_chic_flow.png'])
    
    % Plot stage loading coefficients
    figure(); hold on; grid on; box on; xlabel('Fan Diameter / m'); ylabel('Stage Loading Coefficient');
    plot(r*2,psi_des,'-'); plot(r*2,psi_sta,'-'); legend('Cruise','Static','location','southwest');
    v = axis; axis([v(1:2) 0 v(4)]);
%     export_pres(gcf,[dr.fig 'ideal_chic_loading_1.png'])
    
    % Plot stage loading increase
    figure(); hold on; grid on; box on; 
    xlabel('Fan Diameter / m'); ylabel('Required Stage Loading Increase');
    plot(r*2,(psi_sta - psi_des) ./ psi_des,'-');
%     export_pres(gcf,[dr.fig 'ideal_chic_loading_2.png'])
    
end


%%% Processing Fan Characteristic Data %%%
% 01/12/2020 - Sam Drury %
% This takes the RPM, Phi, Psi_ts, and Psi_tt data from the EDF, and
% outputs the fan characteristic, showing stall

clear; clc; close;
%% INPUTS

%Declare variables
rho = 1.225;
mu = 1.81e-5;

% Prop 1 (grey) = P1
% Prop 2 (green) = P2
% EDF

P1.ri = 35/2/1000;  P1.ro = 10*25.4/2/1000;   P1.c = 29/1000;
P2.ri = 35/2/1000;  P2.ro = 10*25.4/2/1000;   P2.c = 27.18/1000;
EDF.ri = 20/1000;    EDF.ro = 60/1000;      EDF.c = 20/1000;
EDF.sigma = 1.1;

%% DATA ENTRY

% Start from slowest entry - ensure consistent dimensions
% P1.FOM = [.462, .432, .422];
% P2.FOM = [.430, .410, .402];
% EDF.FOM   = [.795, .725, .691, .671];
% 
% P1.RPM = [5221, 5595, 5811];
% P2.RPM = [5253, 5550, 5697];
% EDF.RPM   = [5813, 6083, 6375, 6830];
% 
% P1.T   = [503, 590, 629]*9.81/1000;
% P2.T   = [491, 579, 611]*9.81/1000;
% EDF.T  = [301, 386, 454, 481]*9.81/1000;

P1.FOM = [.432, .422];
P2.FOM = [.410, .402];
EDF.FOM   = [.725, .691, .671];

P1.RPM = [5595, 5811];
P2.RPM = [5550, 5697];
EDF.RPM   = [6083, 6375, 6830];

P1.T   = [590, 629]*9.81/1000;
P2.T   = [579, 611]*9.81/1000;
EDF.T  = [386, 454, 481]*9.81/1000;

%% DATA PROCESSING

P1.Ax = pi*(P1.ro^2 - P1.ri^2);
P2.Ax = pi*(P2.ro^2 - P2.ri^2);
EDF.Ax = pi*(EDF.ro^2 - EDF.ri^2);

P1.Vx = sqrt(P1.T./rho./P1.Ax);
P2.Vx = sqrt(P2.T./rho./P2.Ax);
EDF.Vx = sqrt(EDF.T.*EDF.sigma./rho./EDF.Ax);

P1.rm = rms([P1.ri, P1.ro]);
P2.rm = rms([P2.ri, P2.ro]);
EDF.rm = rms([EDF.ri, EDF.ro]);

P1.U = P1.rm*P1.RPM*2*pi/60;
P2.U = P2.rm*P2.RPM*2*pi/60;
EDF.U = EDF.rm*EDF.RPM*2*pi/60;

P1.V1_rel = sqrt(P1.U.^2 + P1.Vx.^2);
P2.V1_rel = sqrt(P2.U.^2 + P2.Vx.^2);
EDF.V1_rel = sqrt(EDF.U.^2 + EDF.Vx.^2);

P1.Re = rho*P1.V1_rel*P1.c/mu;
P2.Re = rho*P2.V1_rel*P2.c/mu;
EDF.Re = rho*EDF.V1_rel*EDF.c/mu;

mins = [min(P1.Re),min(P2.Re),min(EDF.Re)];
maxs = [max(P1.Re),max(P2.Re),max(EDF.Re)];
Re_span = [min(mins), max(maxs)];

P1.FOM_m = rms(P1.FOM)*ones(2);
P2.FOM_m = rms(P2.FOM)*ones(2);
EDF.FOM_m = rms(EDF.FOM)*ones(2);

err = [std(P2.FOM), std(P1.FOM), std(EDF.FOM)];

%% PLOTTING

% Plotting all data
% figure(1); hold on; grid on; %box on;
% plot(P1.Re, P1.FOM, 'ko'); plot(P2.Re, P2.FOM, 'bo'); plot(EDF.Re, EDF.FOM, 'ro');
% plot(Re_span, P1.FOM_m, '-k'); plot(Re_span, P2.FOM_m, '-b'); plot(Re_span, EDF.FOM_m, '-r');
% xlabel('Re'); ylabel('M_F');
% legend('Prop 1', 'Prop 2', 'EDF');%, 'Prop 1 mean', 'Prop 2 mean', 'EDF mean');

figure(2);
FOMs = [P2.FOM_m(1), P1.FOM_m(1), EDF.FOM_m(1)]
labels = {'Amazon Prop', 'APC VTOL Prop', 'Ducted Fan'};
x = categorical(labels);
x = reordercats(x,labels);
b=bar(x,FOMs,0.5,'FaceColor','flat', 'EdgeColor', 'none');
ylabel('M_F');

% Apply error bars
hold on;
er = errorbar(x,FOMs,err, 'LineWidth', 1);    
er.Color = [0 0 0]; er.LineStyle = 'none';  

% Set colours of bars
b.CData(1,:) = [0 0.4470 0.7410]; b.CData(2,:) = [0.8500 0.3250 0.0980]; 
b.CData(3,:) = [0.9290 0.6940 0.1250];

set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.6 0.7]);
hold off;

% Save fig
%saveas(gcf, 'FOM_comp.eps', 'epsc'); 
%exportgraphics(gcf, 'FOM_comp.png', 'Resolution', 600);

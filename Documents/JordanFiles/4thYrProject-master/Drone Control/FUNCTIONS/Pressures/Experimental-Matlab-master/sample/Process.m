% Compare retests of design point losses

clear all
close all
clc

global Pref Tref
Pref =  101325.8; 
Tref = 288.15; 

Gibbons_Exp

%% List of complete area traverses and chics

% Design        100% nrt = 1.079        100% mrtp = 0.4746          100% phi = 0.5891
% High mrtp     100% nrt = 1.079        102% mrtp = 0.4841
% Low nrt        98% nrt = 1.057        mrtp = 0.4635               100% phi = 0.5891

% Includes:
    % x1 Fine 2 passage T5 traverse at ~100% nrt & 100% mrtp            #1
    % x1 Coarse 1 passage T5 traverse at ~98% nrt & 100% phi            #2
    % x1 Coarse 1 passage T5 traverse at ~100% nrt & 102% mrtp          #3
    % x1 Fine 2 passage T4 traverse at ~100% nrt & 100% mrtp            #4
    % x1 Coarse 1 passage T4 traverse at ~98% nrt & 100% phi            #5
    % x1 Coarse 1 passage T4 traverse at ~100% nrt & 102% mrtp          #6
    % x1 Sector chics at ~100% nrt                                      #7
    % x1 Sector chics at ~98% nrt                                       #8
    % x1 T4 incidence chics at ~100% nrt                                #9
    % x1 T4 incidence chics at ~98% nrt                                 #10
    
% Design point data
des.mrtp = 0.4746; des.nrt = 1.079; des.phi = 0.5891; des.psi = 0.243;
des.P = 101210; des.T = 287.15; des.U_mid = r_mid * 2 * pi * 3500 / 60;

% Coordinates
r.hub = r_hub; r.mid = r_mid; r.cas = r_cas;
  
% Test array size
p_3d = [0 6 14 20 26 34 40]; p_sc = [100 108 116];
nn = length(p_sc); nm = length(p_3d);
f.tests = cell(nn,nm,3);

% Colour data
Re = [linspace(1,0,round(nm/2)) zeros(1, nm-round(nm/2))]';
Gr = [zeros(1, nm-round(nm/2)) linspace(0,1,round(nm/2))]';
Bl = linspace(0,1,round(nm/2)); Bl = [Bl(1:end) Bl(end-1:-1:1)]';
Col = repmat(reshape([Re Gr Bl],[1 nm 3]),[nn 1 1]) .* repmat([0.4 0.7 1]',[1 nm 3]);
Col(1,4,:) = [0 0 0];
Col(3,4,:) = [0.7 0.7 0];

% Useful channels
N.Po_in = 6; N.To_in = 2; N.P_in = 8; N.P_out = 14; N.P_out_fhp = 15; N.P_in_fhp = 10; 
N.P_Sin = 1; N.P_Rin = 9; N.Pa = 16; 
N.FHP = 1:5; N.T_fhp = 5;
N.P_ref = 10; N.T_ref = 4; 
N.rig.P = 8; N.rig.Po = 6; N.rig.T = 2; N.rig.P_out = 14; N.rig.Pa = 16;

%% Current design loss testing

% N75 Tests
% f.tests(1,4,:) = {'SLA_S100_L100_N75_54x27_Ph059_02-Jun-2015_T5_FH_2'...
%     'SLA_S100_L100_N75_15x27_Ph059_02-Jun-2015_T4_FH_5'...
%     'SLA_S100_L100_N75_incidence_02-Jun-2015_1_chic'};                      % Probe calibration error
f.tests(1,4,:) = {'SLA_S100_L100_N75_54x27_Ph059_07-Jun-2015_T5_FH_1'...
    'SLA_S100_L100_N75_15x27_Ph059_07-Jun-2015_T4_FH_1'...
    'SLA_S100_L100_N75_incidence_07-Jun-2015_1_chic'};

% N70 Tests
f.tests(2,2,:) = {'SLA_S30_L30_N70_54x27_Ph059_05-Jun-2015_T5_FH_2'...
    'SLA_S30_L30_N70_15x27_Ph059_05-Jun-2015_T4_FH_1'...
    'SLA_S30_L30_N70_incidence_05-Jun-2015_3_chic'};

% f.tests(2,4,:) = {'SLA_S100_L100_N70_54x27_Ph059_02-Jun-2015_T5_FH_2'...
%     'SLA_S100_L100_N70_15x27_Ph059_02-Jun-2015_T4_FH_1'...
%     'SLA_S100_L100_N70_incidence_03-Jun-2015_2_chic'};                      % Lots of leakage
f.tests(2,4,:) = {'SLA_S100_L100_N70_54x27_Ph059_04-Jun-2015_T5_FH_1'...
    'SLA_S100_L100_N70_15x27_Ph059_04-Jun-2015_T4_FH_1'...
    'SLA_S100_L100_N70_incidence_04-Jun-2015_1_chic'};

% f.tests(2,6,:) = {'SLA_S170_L170_N70_54x27_Ph059_04-Jun-2015_T5_FH_1'...
%     'SLA_S170_L170_N70_15x27_Ph059_04-Jun-2015_T4_FH_1'...
%     'SLA_S170_L170_N70_incidence_04-Jun-2015_1_chic'};                      % Odd steps in operating point
f.tests(2,6,:) = {'SLA_S170_L170_N70_54x27_Ph059_04-Jun-2015_T5_FH_2'...
    'SLA_S170_L170_N70_15x27_Ph059_04-Jun-2015_T4_FH_2'...
    'SLA_S170_L170_N70_incidence_04-Jun-2015_1_chic'};
% f.tests(2,6,:) = {'SLA_S170_L170_N70_54x27_Ph059_07-Jun-2015_T5_FH_1'...    % Retest with leakage
%     'SLA_S170_L170_N70_15x27_Ph059_04-Jun-2015_T4_FH_2'...
%     'SLA_S170_L170_N70_incidence_04-Jun-2015_1_chic'};

% N65 Tests
f.tests(3,4,:) = {'SLA_S100_L100_N65_54x27_Ph059_05-Jun-2015_T5_FH_1'...
    'SLA_S100_L100_N65_15x27_Ph059_05-Jun-2015_T4_FH_1'...
    'SLA_S100_L100_N65_incidence_05-Jun-2015_1_chic'};

f.calib{1,4} = 'FHP_G_61_fine_05-Jun-2015_Bingo_3.mat';
f.calib{2,2} = 'FHP_G_61_fine_05-Jun-2015_Bingo_3.mat';
f.calib{2,4} = 'FHP_G_61_fine_02-Jun-2015_Bingo_2.mat';
f.calib{2,6} = 'FHP_G_61_fine_05-Jun-2015_Bingo_2.mat';
f.calib{3,4} = 'FHP_G_61_fine_05-Jun-2015_Bingo_3.mat';

%% Current off design traverses

% SLA_S30_L30_N70_7x15_Ph046_05-Jun-2015_T4_FH_1
% SLA_S30_L30_N70_7x15_Ph049_05-Jun-2015_T4_FH_1
% SLA_S30_L30_N70_54x27_Ph046_05-Jun-2015_T5_FH_1
% SLA_S30_L30_N70_54x27_Ph049_05-Jun-2015_T5_FH_1

% SLA_S170_L170_N70_7x15_Ph046_04-Jun-2015_T4_FH_1
% SLA_S170_L170_N70_7x15_Ph048_04-Jun-2015_T4_FH_1
% SLA_S170_L170_N70_54x27_Ph046_04-Jun-2015_T5_FH_2
% SLA_S170_L170_N70_54x27_Ph048_04-Jun-2015_T5_FH_1

%% Test results

% Allocate storage cell arrays
C = cell(nn,nm); G = cell(nn,nm); L = cell(nn,nm);

% Coordinates
x = [x_5 x_4 nan(1,4)];

% Open chic figure window
h = figure(); 
subplot(2,1,1); hold on; grid on; box on;
subplot(2,1,2); hold on; grid on; box on;

% Initialise traverse figure windows
H = [];

% Loop over all geometries
for n = 1:nn
    for m = 1:nm
        
        % Check if chic data exists
        if isempty(f.tests{n,m,3}) == 0
            
            % Calculate incidence characteristics
            C{n,m} = exp_chic_incidence(directory,f.tests{n,m,3},[],N,r_mid,[],[0 0 0],0);
            
            % Plot pressure rise against incidence
            figure(h);
            subplot(2,1,1); plot(C{n,m}.phi_av,C{n,m}.psi_av,'-','Color',squeeze(Col(n,m,:)))
            subplot(2,1,2); plot(C{n,m}.Alpha_av,C{n,m}.psi_av,'-','Color',squeeze(Col(n,m,:)))
        end
        
        % Check if loss data exists
        if isempty(f.tests{n,m,1}) == 0
            
            % Calculate loss traverses
            [L{n,m},G{n,m},H] = exp_traverse_loss(directory,f.tests(n,m,1:3),f.calib(n,m),...
                N,r,x,des,H,squeeze(Col(n,m,:)),1);
        end
    end
end
gtg
%% Five hole probe rubbing testing
f.names = {'SLA_S100_L100_N75_7x15_Ph059_27-May-2015_T4_FH_1'...
    'SLA_S100_L100_N75_7x15_Ph059_27-May-2015_T4_FH_2'...
    'SLA_S100_L100_N75_15x15_Ph059_27-May-2015_T5_FH_1'...
    'SLA_S100_L100_N75_15x15_Ph059_27-May-2015_T5_FH_2'...
    'SLA_S100_L100_N75_7x15_Ph059_27-May-2015_T4_FH_3'...
    'SLA_S100_L100_N75_15x15_Ph059_27-May-2015_T5_FH_3'...
    'SLA_S100_L100_N75_7x15_Ph059_27-May-2015_T4_FH_4'...
    'SLA_S100_L100_N75_15x15_Ph059_27-May-2015_T5_FH_4'...
    'SLA_S100_L100_N75_7x15_Ph059_27-May-2015_T4_FH_5'...
    'SLA_S100_L100_N75_7x15_Ph059_01-Jun-2015_T4_FH_1'...
    'SLA_S100_L100_N75_7x15_Ph059_01-Jun-2015_T4_FH_2'...
    'SLA_S100_L100_N75_15x15_Ph059_01-Jun-2015_T5_FH_1'...
    'SLA_S100_L100_N75_7x15_Ph059_01-Jun-2015_T4_FH_3'...
    'SLA_S100_L100_N75_15x27_Ph059_02-Jun-2015_T4_FH_1'...
    'SLA_S100_L100_N75_54x27_Ph059_02-Jun-2015_T5_FH_1'...
    'SLA_S100_L100_N75_15x27_Ph059_02-Jun-2015_T4_FH_3'...
    'SLA_S100_L100_N75_15x27_Ph059_02-Jun-2015_T4_FH_4'...
    'SLA_S100_L100_N75_54x27_Ph059_02-Jun-2015_T5_FH_2'...
    'SLA_S100_L100_N75_15x27_Ph059_02-Jun-2015_T4_FH_5'};

N_plot = length(f.names)-4:length(f.names);
C = lines(length(N_plot));
h.probes = figure(); 
subplot(1,4,1); hold on; grid on; box on; xlabel('Po');
subplot(1,4,2); hold on; grid on; box on; xlabel('Yaw');
subplot(1,4,3); hold on; grid on; box on; xlabel('Pitch');
subplot(1,4,4); hold on; grid on; box on; xlabel('Error');
h.P_raw_radial = figure(); for m = 1:5; subplot(1,5,m); hold on; grid on; box on; end;
h.P_raw_tangential = figure(); for m = 1:5; subplot(5,1,m); hold on; grid on; box on; end;
for n = 1:length(N_plot)
    load([directory f.names{N_plot(n)}])
    g = exp_traverse_calculate(e,{e.c},x_4,r,N,des);
    figure(h.probes)
    subplot(1,4,1); plot(ts_mass_average(g,'Po',2),ts_mass_average(g,'r_nondim',2),'-','Color',C(n,:));
    subplot(1,4,2); plot(ts_mass_average(g,'Alpha',2),ts_mass_average(g,'r_nondim',2),...
        '-','Color',C(n,:));
    subplot(1,4,3); l(n) = plot(ts_mass_average(g,'Beta',2),ts_mass_average(g,'r_nondim',2),...
        '-','Color',C(n,:));
    subplot(1,4,4); plot(ts_mass_average(g,'A',2),ts_mass_average(g,'r_nondim',2),'-','Color',C(n,:));
    figure(h.P_raw_radial);
    for m = 1:5; subplot(1,5,m); plot(mean(g.P_corr(:,:,m),2),g.r(:,1),'-','Color',C(n,:)); end;
    figure(h.P_raw_tangential);
    for m = 1:5; subplot(5,1,m); plot(g.t(1,:),mean(g.P_corr(:,:,m),1),'-','Color',C(n,:)); end;
end
for n = 1:length(N_plot)
    L{n} = num2str(n);
end
legend(l,L)
gtg

%% Previous loss testing

% Useful channels
N.Po_in = 6; N.To_in = 2; N.P_in = 8; N.P_out = 2; 
N.P_Sin = 1; N.P_Rin = 9; N.Pa = 16; 
N.FHP = 1:5; N.T_fhp = 5; N.P_out_fhp = 14;
N.P_ref = 12; N.T_ref = 5; 
N.rig.P = 8; N.rig.Po = 6; N.rig.T = 2; N.rig.P_out = 14; N.rig.Pa = 16;
    
f.tests = cell(nn,nm,10);

% Calibration data
% f.calib{1} = 'FHP_G_27_55_fine_02-Apr-2015_Bingo_1.mat';
f.calib{1} = 'FHP_G_27_41_55_fine_15-Apr-2015_Bingo_1_flip.mat';

% S100_N75
f.tests(1,4,:) = {'SLA_S100_L100_N75_54x27_Ph059_10-Apr-2015_T5_FH_1'
    'SLA_S100_L100_N75_15x15_Ph059_10-Apr-2015_T5_FH_1'
    'SLA_S100_L100_N75_15x15_Ph060_10-Apr-2015_T5_FH_1'
    'SLA_S100_L100_N75_11x27_Ph059_10-Apr-2015_T4_FH_1'
    'SLA_S100_L100_N75_7x15_Ph059_10-Apr-2015_T4_FH_1'
    'SLA_S100_L100_N75_7x15_Ph060_10-Apr-2015_T4_FH_1'
    'SLA_S100_L100_N75_sector_10-Apr-2015_3_chic'
    'SLA_S100_L100_N75_sector_10-Apr-2015_4_chic'
    'SLA_S100_L100_N75_incidence_10-Apr-2015_1_chic'
    'SLA_S100_L100_N75_incidence_10-Apr-2015_2_chic'};

% S30_N70
f.tests(2,2,:) = {'SLA_S30_L30_N70_54x27_Ph059_09-Apr-2015_T5_FH_1'
    'SLA_S30_L30_N70_15x15_Ph059_09-Apr-2015_T5_FH_2'
    'SLA_S30_L30_N70_15x15_Ph060_09-Apr-2015_T5_FH_1'
    'SLA_S30_L30_N70_11x27_Ph059_09-Apr-2015_T4_FH_1'
    'SLA_S30_L30_N70_7x15_Ph059_09-Apr-2015_T4_FH_1'
    'SLA_S30_L30_N70_7x15_Ph060_09-Apr-2015_T4_FH_1'
    'SLA_S30_L30_N70_sector_09-Apr-2015_2_chic';...         % Wrong plumbing
    'SLA_S30_L30_N70_sector_09-Apr-2015_1_chic';...         % Wrong plumbing
    'SLA_S30_L30_N70_incidence_09-Apr-2015_1_chic'
    'SLA_S30_L30_N70_incidence_09-Apr-2015_2_chic'};

% S100_N70
f.tests(2,4,:) = {'SLA_S100_L100_N70_54x27_Ph059_10-Apr-2015_T5_FH_1'
    'SLA_S100_L100_N70_15x15_Ph059_10-Apr-2015_T5_FH_1'
    'SLA_S100_L100_N70_15x15_Ph060_10-Apr-2015_T5_FH_1'
    'SLA_S100_L100_N70_11x27_Ph059_10-Apr-2015_T4_FH_1'
    'SLA_S100_L100_N70_7x15_Ph059_10-Apr-2015_T4_FH_1'
    'SLA_S100_L100_N70_7x15_Ph060_10-Apr-2015_T4_FH_1'
    'SLA_S100_L100_N70_sector_10-Apr-2015_1_chic'
    'SLA_S100_L100_N70_sector_10-Apr-2015_2_chic'
    'SLA_S100_L100_N70_incidence_10-Apr-2015_2_chic'
    'SLA_S100_L100_N70_incidence_10-Apr-2015_1_chic'};

% S170_N70
f.tests(2,6,:) = {'SLA_S170_L170_N70_54x27_Ph059_09-Apr-2015_T5_FH_2'
    'SLA_S170_L170_N70_15x15_Ph059_09-Apr-2015_T5_FH_1'
    'SLA_S170_L170_N70_15x15_Ph060_09-Apr-2015_T5_FH_1'
    'SLA_S170_L170_N70_11x27_Ph059_09-Apr-2015_T4_FH_1'
    'SLA_S170_L170_N70_7x15_Ph059_09-Apr-2015_T4_FH_1'
    'SLA_S170_L170_N70_7x15_Ph060_09-Apr-2015_T4_FH_1'
    'SLA_S170_L170_N70_sector_09-Apr-2015_2_chic';...         % Wrong plumbing
    'SLA_S170_L170_N70_sector_09-Apr-2015_1_chic';...         % Wrong plumbing
    'SLA_S170_L170_N70_incidence_09-Apr-2015_1_chic'
    'SLA_S170_L170_N70_incidence_09-Apr-2015_2_chic'};


% S0_N65
f.tests(3,1,:) = {'SLA_S0_L0_N65_54x27_Ph059_09-Apr-2015_T5_FH_1'
    'SLA_S0_L0_N65_15x15_Ph059_09-Apr-2015_T5_FH_1'
    'SLA_S0_L0_N65_15x15_Ph060_09-Apr-2015_T5_FH_1'
    'SLA_S0_L0_N65_11x27_Ph059_09-Apr-2015_T4_FH_1'
    'SLA_S0_L0_N65_7x15_Ph059_09-Apr-2015_T4_FH_1'
    'SLA_S0_L0_N65_7x15_Ph060_09-Apr-2015_T4_FH_1'
    'SLA_S0_L0_N65_sector_09-Apr-2015_3_chic'
    'SLA_S0_L0_N65_sector_09-Apr-2015_2_chic'
    'SLA_S0_L0_N65_incidence_09-Apr-2015_1_chic'
    'SLA_S0_L0_N65_incidence_09-Apr-2015_2_chic'};

% S100_N65
f.tests(3,4,:) = {'SLA_S100_L100_N65_54x27_Ph059_08-Apr-2015_T5_FH_3'
    'SLA_S100_L100_N65_15x15_Ph059_08-Apr-2015_T5_FH_1'
    'SLA_S100_L100_N65_15x15_Ph060_08-Apr-2015_T5_FH_1'
    'SLA_S100_L100_N65_11x27_Ph059_08-Apr-2015_T4_FH_1'
    'SLA_S100_L100_N65_7x15_Ph059_08-Apr-2015_T4_FH_1'
    'SLA_S100_L100_N65_7x15_Ph060_08-Apr-2015_T4_FH_1'
    'SLA_S100_L100_N65_sector_08-Apr-2015_2_chic'
    'SLA_S100_L100_N65_sector_08-Apr-2015_1_chic'
    'SLA_S100_L100_N65_incidence_08-Apr-2015_1_chic'
    'SLA_S100_L100_N65_incidence_08-Apr-2015_2_chic'};

% S200_N65
f.tests(3,7,:) = {'SLA_S200_L200_N65_54x27_Ph059_08-Apr-2015_T5_FH_3'
    'SLA_S200_L200_N65_15x15_Ph059_08-Apr-2015_T5_FH_1'
    'SLA_S200_L200_N65_15x15_Ph060_08-Apr-2015_T5_FH_1'
    'SLA_S200_L200_N65_11x27_Ph059_08-Apr-2015_T4_FH_1'
    'SLA_S200_L200_N65_7x15_Ph059_08-Apr-2015_T4_FH_1'
    'SLA_S200_L200_N65_7x15_Ph060_08-Apr-2015_T4_FH_1'
    'SLA_S200_L200_N65_sector_08-Apr-2015_2_chic'
    'SLA_S200_L200_N65_sector_08-Apr-2015_1_chic'
    'SLA_S200_L200_N65_incidence_08-Apr-2015_1_chic'
    'SLA_S200_L200_N65_incidence_08-Apr-2015_2_chic'};

% Replace N65 tests with retests
f.tests{3,7,1} = 'SLA_S200_L200_N65_54x27_Ph059_21-Apr-2015_T5_FH_2';
f.tests{3,7,4} = 'SLA_S200_L200_N65_17x27_Ph059_21-Apr-2015_T4_FH_1';
f.tests{3,1,1} = 'SLA_S0_L0_N65_54x27_Ph059_21-Apr-2015_T5_FH_1';
f.tests{3,1,4} = 'SLA_S0_L0_N65_17x27_Ph059_21-Apr-2015_T4_FH_1';
f.tests{3,4,1} = 'SLA_S100_L100_N65_54x27_Ph059_21-Apr-2015_T5_FH_1';
f.tests{3,4,4} = 'SLA_S100_L100_N65_17x27_Ph059_21-Apr-2015_T4_FH_1';

% Replace inlet traverses with datum traverse
% for n = 1:size(f.tests,1)
%     for m = 1:size(f.tests,2)
%         if isempty(f.tests{n,m}) == 0
%            f.tests{n,m,4} = 'SLA_S100_L100_N75_11x27_Ph059_10-Apr-2015_T4_FH_1';
%         end
%     end
% end

%% Incidence characteristic plots

% % Allocate chic storage cell array
% C = cell(nn,nm);
% 
% % Open figure window
% figure(); 
% subplot(2,1,1); hold on; grid on; box on;
% subplot(2,1,2); hold on; grid on; box on;
% 
% % Loop over all geometries
% for n = 1:nn
%     for m = 1:nm
%         
%         % Check if test data exists
%         if isempty(f.tests{n,m,9}) == 0
%             
%             % Choose calibration map
%             
%             
%             % Calculate incidence characteristics
%             C{n,m} = exp_chic_incidence(directory,f.tests(n,m,[7 9]),f.calib{1},...
%                 N,r_mid,[],[0 0 0],0);
%             
%             % Plot pressure rise against incidence
%             subplot(2,1,1); plot(C{n,m}.phi,C{n,m}.psi,'-','Color',squeeze(Col(n,m,:)))
%             subplot(2,1,2); plot(C{n,m}.Alpha,C{n,m}.psi,'-','Color',squeeze(Col(n,m,:)))
%         end
%     end
% end

%% Loss traverse plots

% Allocate traverse storage cell array
G = cell(nn,nm); L = cell(nn,nm);

% Coordinates
x = [x_5*ones(1,3) x_4 * ones(1,3) nan(1,4)];

% Loop over all geometries
H = [];
for n = 2%1:nn
    for m = 1:nm
        
        % Check if test data exists
        if isempty(f.tests{n,m,1}) == 0
            
            % Calculate loss traverses
            [L{n,m},G{n,m},H] = exp_traverse_loss(directory,f.tests(n,m,1:8),f.calib,...
                N,r,x,des,H,squeeze(Col(n,m,:)),1);
        end
    end
end

%% Plot loss traverses for ASME presentation

h = figure('position',[228 126 1511 721]); 
subplot(1,2,1); hold on; axis equal; 
subplot(1,2,2); hold on; axis equal; 
n = 2; ms = [2 6]; p = [1 2];
for o = 1:length(ms)
    Po_in = ts_mass_average(G{n,ms(o)}{4},'Po',3);
    P_in = ts_mass_average(G{n,ms(o)}{4},'P',3);
    Po_out = G{n,ms(o)}{1}.Po; Po_max = ts_mass_average(G{n,ms(o)}{1},'Po_max',3);
    z = G{n,ms(o)}{1}.z; y = G{n,ms(o)}{1}.y; t = G{n,ms(o)}{1}.t;
    Yp = (Po_max - Po_out) ./ (Po_in - P_in);
    subplot(1,2,p(o)); contourf(z,y,Yp,8); colormap(flipud(bone));
    axis([-0.04 0.04 0.305 0.365]);
    [z_hub,y_hub] = pol2cart(linspace(t(1,1),t(1,end),100)+pi/2,ones(1,100)*r.hub);
    [z_cas,y_cas] = pol2cart(linspace(t(end,1),t(end,end),100)+pi/2,ones(1,100)*r.cas);
    plot(-z_hub,y_hub,'k--'); plot(-z_cas,y_cas,'k--');
    (ts_mass_average(G{n,ms(o)}{1},'Po_max',3) - ts_mass_average(G{n,ms(o)}{1},'Po',3)) ./ (Po_in - P_in)
end
colorbar('WestOutside');
gtg
export_fig(h,'Figures\Exp_Loss_Wakes.eps','-painters');
gtg
h = figure(); hold on; grid on; box on; xlabel('Loss Coefficient'); ylabel('Radius');
n = 2; ms = [2 4 6];
for m = ms
    plot(L{n,m}.Yp/0.0578,L{n,m}.r,'--','LineWidth',2,'Color',squeeze(Col(end,m,:)));
end

% Load in TURBOSTREAM solutions
p = [30 100 170]; C = [1 0 0 ; 0 0 1 ; 0 1 0];
for m = 1:length(p)
    c = gibbons_cut(cfd_directory,['RS_S_F2_S' num2str(p(m))...
        '_L' num2str(p(m)) '_N70_DRR_CRR_H0_Y1_p1s_avg.hdf5']);
    
    % Average pressures and radii
    [Po_in, mass_in] = ts_mass_average(c.S_in,'Po',2);
    [Po_out, mass_out] = ts_mass_average(c.S_out,'Po',2);
    r_out = ts_mass_average(c.S_out,'r_nondim',2);
    dP = ts_mass_average(c.S_in,'Po',3) - ts_mass_average(c.S_in,'P',3);
    
    % Calculate mass fractions
    mass_in = [0 cumsum(mass_in)]; mass_out = [0 cumsum(mass_out)];
    mass_in = mass_in / max(mass_in); mass_out = mass_out / max(mass_out);
    
    % Radial loss coefficient based on finite mass fractions
%     m_frac = 0.1; Po_in_frac = zeros(length(Po_in),1);
%     for j = 1:length(Po_out)
%         m_temp = linspace(max(mass_in(j)-m_frac,0),min(mass_in(j)+m_frac,1),10);
%         Po_in_frac(j) = mean(interp1(mass_in,Po_in,m_temp,'pchip'));
%     end
%     Yp = (interp1(mass_in,Po_in_frac,mass_out,'pchip') - Po_out) / dP;
    
    % Radial loss coefficient based on inlet average
%     Yp = (ts_mass_average(c.S_in,'Po',3) - Po_out) / dP;

    % Radial loss coefficient based on cumilative mass fraction
    Yp = (interp1(mass_in,Po_in,mass_out,'pchip') - Po_out) / dP;
    
    plot(Yp/0.0578,r_out,'-','LineWidth',2,'Color',C(m,:));
end
axis([0 2.5 0 0.7]);

saveas(h,[directory '../Pres_CFD_Exp_Loss.fig'])
gtg
%% Save experimental data for Chris

% % Save matlab files
% G = G(2,[2 4 6]);
% for n = 1:length(G)
%     G{n} = G{n}([1 4 7]);
% end
% save([directory 'N70_Exp_Data.mat'],'G');
% 
% % Save excel files
% p_3d = [30 100 170];
% for n = 1:length(G)
%     
%     % Save traverses
%     varnames = {'x' 'r' 'rt' 'y' 'z' 'ro' 'rovx' 'rovr' 'rorvt' 'roe' 'P' 'Po'};
%     for m = 1:2
%         for v = 1:length(varnames)
%             % Processed data
%             xlswrite([directory 'N70_S' num2str(p_3d(n)) '.xls'],varnames(v),...
%                 m,['A' num2str(v*28 - 27)]);
%             xlswrite([directory 'N70_S' num2str(p_3d(n)) '.xls'],G{n}{m}.(varnames{v}),...
%                 m,['A' num2str(v*28 - 26)]);
%         end
%     end
%     
%     % Save characteristics
%     m = 3; varnames = fieldnames(G{n}{3}); l = 1;
%     for v = 1:length(varnames)
%         xlswrite([directory 'N70_S' num2str(p_3d(n)) '.xls'],varnames(v),...
%             m,[xlscol(l) '1']);
%         xlswrite([directory 'N70_S' num2str(p_3d(n)) '.xls'],G{n}{m}.(varnames{v}),...
%             m,[xlscol(l) '2']);
%         l = l + size(G{n}{m}.(varnames{v}),2) + 1;
%     end    
% end

%% Probe error plots

% Open figure window
figure(); hold on; grid on; box on; xlabel('Time'); ylabel('Average Probe Error');

% Loop over all traverses and average error data
time = []; A = []; C = [];
for n = 1:nn
    for m = 1:nm
        
        % Check if test data exists
        if isempty(f.tests{n,m,1}) == 0
            
            for o = [1 4]%1:6
                % Calculate mass averaged error
                A = [A ; ts_mass_average(G{n,m}{o},'A',3)];

                % Get file date
                file_info = dir([directory f.tests{n,m,o} '.mat']);
                time = [time ; file_info.datenum];
                
                % Record colour
                C = [C ; squeeze(Col(n,m,:))'];
            end
        end
    end
end

% Sort in time order
[time,i] = sort(time);
A = A(i); C = C(i,:);

% Plot error increase
plot(time,A,'k.-');

for n = 1:length(time)
    plot(time(n),A(n),'.','Color',C(n,:));
end


% exp_traverse_loss(directory,f_retest.tests(3,7,1:8),f.calib,N,r,x,des,H,[0 0.5 0],1);
% exp_traverse_loss(directory,f_retest.tests(3,4,1:8),f.calib,N,r,x,des,H,[0 0 0.5],1);
% exp_traverse_loss(directory,f_retest.tests(3,1,1:8),f.calib,N,r,x,des,H,[0.5 0 0],1);

%% Test on drift of inlet pressure measurements

load([directory f.calib])

names = {'SLA_S100_L100_N75_11x27_Ph059_10-Apr-2015_T4_FH_1' ;...
    'SLA_S100_L100_N75_11x27_Ph059_22-Apr-2015_T4_FH_1' ;...
    'SLA_S100_L100_N75_11x27_Ph059_23-Apr-2015_T4_FH_1' ;...
    'SLA_S100_L100_N75_11x27_Ph059_23-Apr-2015_T4_FH_2';...
    'SLA_S100_L100_N75_11x27_Ph059_23-Apr-2015_T4_FH_3'};

colours = [0 0 0 ; 1 0 0 ; 0 0 1 ; 0 1 0 ; 0 0.5 1];
legends = {'Original' 'Hot' 'Cold' 'Retest Cold #1' 'Retest Cold #2'};

h.profile = figure(); hold on; box on; grid on; xlabel('Stagnation Pressure'); ylabel('Radius');
h.flow = figure(); hold on; box on; grid on; xlabel('Flow Coefficient'); ylabel('Radius');
h.cen = figure(); 
for n = 1:5
    subplot(1,5,n); hold on; box on; grid on; xlabel('P probe'); ylabel('Radius');
end
for n = 1:length(names)
    g = exp_traverse_calculate([directory names{n}],c,x_5,r,N,des);
    figure(h.profile)
    plot(ts_mass_average(g,'Po',2),ts_mass_average(g,'r_nondim',2),'-','Color',colours(n,:))
    figure(h.flow)
    plot(g.rovx,g.r,'-','Color',colours(n,:))
    figure(h.cen)
    for m = 1:5
        subplot(1,5,m)
        plot(g.P_corr(:,:,m),g.r,'-','Color',colours(n,:))
    end
%     C_Alpha = 0.5 * (P_ur + P_dr - P_ul - P_dl) ./ (P_cen - P_av);
end
figure(h.profile)
legend(legends,'location','SouthWest');

% Plot individual pressure readings and pitch and yaw coefficients
h.cen = figure(); 
for n = 1:5
    subplot(1,5,n); hold on; box on; grid on; xlabel('P probe'); ylabel('Radius');
end
% h.calib = figure(); hold on; box on; grid on;
% plot(c.C_Alpha(:,:,end),c.C_Beta(:,:,end),'k-')
% plot(c.C_Alpha(:,:,end)',c.C_Beta(:,:,end)','k-')
for n = 1:size(G,1)
    for m = 1:size(G,2)
        if isempty(G{n,m}) == 0
            h.cen = figure();
            for p = 1:5
                figure(h.cen); subplot(1,5,p); hold on; box on; grid on;
%                 plot(mean(G{n,m}{4}.P_corr(:,:,p),2),mean(G{n,m}{4}.r,2),'-',...
%                     'Color',squeeze(Col(n,m,:)))
                plot(G{n,m}{4}.P_corr(:,:,p),G{n,m}{4}.r,'-','Color',squeeze(Col(n,m,:)))
%                 plot(mean(G{n,m}{1}.P_corr(:,1,p),2),mean(G{n,m}{1}.r,2),'-',...
%                     'Color',squeeze(Col(n,m,:)))
                plot(G{n,m}{1}.P_corr(:,:,14),G{n,m}{1}.r,'--','Color',squeeze(Col(n,m,:)))
            end
%             figure(h.calib);
%             plot(G{n,m}{4}.C_Alpha,G{n,m}{4}.C_Beta,'.','Color',squeeze(Col(n,m,:)))
        end
    end
end

figure(); hold on; grid on; box on;
N = [75 70 65];
for n = 1:size(G,1)
    for m = 1:size(G,2)
        if isempty(G{n,m}) == 0
            [Po_in,m_in] = ts_mass_average(G{n,m}{4},'Po',2);
            [Beta_in,m_in] = ts_mass_average(G{n,m}{4},'Beta',2);
            [r_in,m_in] = ts_mass_average(G{n,m}{4},'r_nondim',2);
            plot(interp1(r_in,Beta_in,0.5),interp1(r_in,Po_in,0.5),'.','Color',squeeze(Col(n,m,:)))
        end
    end
end

%% Old data

%%%%%%%% Pre Probe Rub Tue 31st Night - Regrind FHP G %%%%%%%

% f.calib = {'FHP_G_20_40_65_fine_14-Jan-2015_Donkey'};

% f.tests(1,4,:) = {'SLA_S100_L100_N75_54x27_Ph059_26-Mar-2015_T5_FH_1'
%     'SLA_S100_L100_N75_15x15_Ph059_27-Mar-2015_T5_FH_1'
%     'SLA_S100_L100_N75_15x15_Ph060_27-Mar-2015_T5_FH_1'
%     'SLA_S100_L100_N75_11x27_Ph059_26-Mar-2015_T4_FH_1'
%     'SLA_S100_L100_N75_7x15_Ph059_26-Mar-2015_T4_FH_1'
%     'SLA_S100_L100_N75_7x15_Ph060_26-Mar-2015_T4_FH_1'
%     'SLA_S100_L100_N75_sector_27-Mar-2015_1_chic'
%     'SLA_S100_L100_N75_sector_27-Mar-2015_2_chic'
%     'SLA_S100_L100_N75_incidence_27-Mar-2015_2_chic'};
% 
% f.tests(3,3,:) = {''
%     ''
%     ''
%     'SLA_S100_L100_N65_11x27_Ph059_31-Mar-2015_T4_FH_2'
%     'SLA_S100_L100_N65_7x15_Ph059_31-Mar-2015_T4_FH_1'
%     'SLA_S100_L100_N65_7x15_Ph060_31-Mar-2015_T4_FH_1'
%     'SLA_S100_L100_N65_sector_31-Mar-2015_1_chic'
%     'SLA_S100_L100_N65_sector_31-Mar-2015_2_chic'
%     'SLA_S100_L100_N65_incidence_31-Mar-2015_1_chic'};

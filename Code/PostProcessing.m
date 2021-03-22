%% Test function to post-process TS files
clear; close all; clc;
startup;
eVTOL;
% g = ts_read_hdf5([dr.ts 'NRF_design_4_2_avg.hdf5']);
% conv = ts_plot_conv([dr.ts 'log_NRF_design_12_2.txt']);

cids.cell_0 = 'Rotor mixing plane';
cids.cell_1 = 'Rotor trailing gedge';
cids.cell_2 = 'Inlet';
cids.cell_3 = 'Rotor leading edge';
cids.cell_4 = 'Rotor shrouds and TE';
cids.cell_5 = 'Rotor+';
cids.cell_6 = 'Rotor-';
cids.cell_7 = 'Exit';
cids.cell_8 = 'Stator trailing edge';
cids.cell_9 = 'Stator mixing plane';
cids.cell_10 = 'Stator leading edge';
cids.cell_11 = 'Stator shrouds and TE';
cids.cell_12 = 'Stator+';
cids.cell_13 = 'Stator-';

%% Questions

% When calculating phi,psi - should I calculate at each point then mass
% average across entire area?

% Which efficiency should I use?

% Blades are quite different - is this an issue, especially when looking at
% chis?
%% Generating contour variables

for o=4:22
g = ts_read_hdf5([dr.ts 'NRF_design_' num2str(o) '_2_avg.hdf5']);
c_in = ts_structured_cut(g,2,'en','en',1,'en',1,'en');
cs_in = ts_secondary(c_in);

c_out = ts_structured_cut(g,7,1,1,1,'en',1,'en');
cs_out = ts_secondary(c_out);

ga = cs_in.bv.ga;

% Efficiency & Ch'ic
Po1 = ts_mass_average(cs_in,'Po',3);
[Po3,mdot,A] = ts_mass_average(cs_out,'Po',3);
ro3 = ts_mass_average(cs_out,'ro',3);

To1 = ts_mass_average(cs_in,'To',3);
To3 = ts_mass_average(cs_out,'To',3);
To3s = To1 * (Po3/Po1)^((ga-1)/ga);
dho = ts_mass_average(cs_out,'ho',3) - ts_mass_average(cs_in,'ho',3);

eta(o-3) = (To3s - To1) / (To3 - To1);
Um = cs_in.bv.omega * sqrt(0.5*(cs_in.r(1,1)^2 + cs_in.r(end,1)^2)); % j=50 is mean line 0.0447
phi(o-3) = mdot/(ro3*A*Um); %%
psi(o-3) = dho/Um^2;
rpm(o-3) = cs_in.bv.omega*60/2/pi;


disp([num2str(o) ' done'])
end


%% Deviation
des = 12;
g = ts_read_hdf5([dr.ts 'NRF_design_' num2str(des) '_2_avg.hdf5']);

% Initialise figure
figure; hold on; grid on; box on;
set(gcf,'color','w','position',[200 100 400 600]); xlim([-10 10]);
xlabel('Flow angle (degrees)'); ylabel('% of span');
title(['Exit flow angle \alpha_3 for design ' num2str(des)]);

% Stator TE
c = ts_structured_cut(g,7,1,1,1,'en',1,'en');
cs = ts_secondary(c);

span = (cs.r(:,1)-cs.r(1,1)) / (cs.r(end,1)-cs.r(1,1));

alpha   = ts_mass_average(cs,'Alpha',2);

plot(alpha,span); 

% End of CV
c = ts_structured_cut(g,7,'en','en',1,'en',1,'en');
cs = ts_secondary(c);

alpha   = ts_mass_average(cs,'Alpha',2);

plot(alpha,span); 
plot(zeros(numel(span)),span,'-r','LineWidth',1.5)

% contourf(cs.z,cs.y,cs.Alpha,25,'edgecolor','none'); colorbar;
% axis equal; %caxis([-10 10]);

%% Rotor stag pressure loss coefficient

des = 16;
g = ts_read_hdf5([dr.ts 'NRF_design_' num2str(des) '_2_avg.hdf5']);

% Initialise figure
figure; hold on; grid on; box on;
set(gcf,'color','w','position',[650 100 450 600]); xlim([0 inf]);
xlabel('(P_{01} - P_{02}) / (P_{01} - P_1)'); ylabel('% of span');
% xlabel('Stagnation pressure loss coefficient'); ylabel('% of span');
title(['Rotor stagnation pressure loss coefficient Y_P for design ' num2str(des)]);

% Stator inlet and exit
c_in = ts_structured_cut(g,2,1,1,1,'en',1,'en');
cs_in = ts_secondary(c_in);
c_out = ts_structured_cut(g,0,'en','en',1,'en',1,'en');
cs_out = ts_secondary(c_out);

Po1 = ts_mass_average(cs_in,'Po_rel',2);
Po2 = ts_mass_average(cs_out,'Po_rel',2);
P1 = ts_mass_average(cs_in,'P',2);
Yp = (Po1-Po2) ./ (Po1-P1);

% ro1 = ts_mass_average(cs_in,'ro',2);
% V1 = ts_mass_average(cs_in,'V_rel',2);
% Yp = (Po1-Po2) ./ (.5*ro1.*V1.^2);

span = (cs_in.r(:,1)-cs_in.r(1,1)) / (cs_in.r(end,1)-cs_in.r(1,1));

plot(Yp,span);

%% Stator stag pressure loss coefficient

des = 16;
g = ts_read_hdf5([dr.ts 'NRF_design_' num2str(des) '_2_avg.hdf5']);

% Initialise figure
figure; hold on; grid on; box on;
set(gcf,'color','w','position',[1050 100 450 600]); xlim([0 inf]);
xlabel('(P_{02} - P_{03}) / (P_{02} - P_2)'); ylabel('% of span');
% xlabel('Stagnation pressure loss coefficient'); ylabel('% of span');
title(['Stator stagnation pressure loss coefficient Y_P for design ' num2str(des)]);

% Stator inlet and exit
c_in = ts_structured_cut(g,9,1,1,1,'en',1,'en');
cs_in = ts_secondary(c_in);
c_out = ts_structured_cut(g,7,'en','en',1,'en',1,'en');
cs_out = ts_secondary(c_out);

Po2 = ts_mass_average(cs_in,'Po',2);
Po3 = ts_mass_average(cs_out,'Po',2);
P2 = ts_mass_average(cs_in,'P',2);
Yp = (Po2-Po3) ./ (Po2-P2);

% ro2 = ts_mass_average(cs_in,'ro',2);
% V2 = ts_mass_average(cs_in,'V',2);
% Yp = (Po2-Po3) ./ (.5*ro2.*V2.^2);

span = (cs_in.r(:,1)-cs_in.r(1,1)) / (cs_in.r(end,1)-cs_in.r(1,1));

plot(Yp,span); 

%% Phi

des = 16;
g = ts_read_hdf5([dr.ts 'NRF_design_' num2str(des) '_2_avg.hdf5']);

% Initialise figure
figure; hold on; grid on; box on;
set(gcf,'color','w','position',[1050 100 400 600]); xlim([0 inf]);
xlabel('\phi'); ylabel('% of span');
title(['Flow coefficient \phi at inlet for design ' num2str(des)]);

c_in = ts_structured_cut(g,2,10,10,1,'en',1,'en'); % needed for omega
cs_in = ts_secondary(c_in);
c = ts_structured_cut(g,9,5,5,1,'en',1,'en'); % 2=inlet 7=exit
cs = ts_secondary(c);

Um = cs_in.bv.omega * sqrt(0.5*(cs_in.r(1,1)^2 + cs_in.r(end,1)^2)); % j=50 is mean line 0.0447
Vx = ts_mass_average(cs,'Vx',2); % change cs to cs_in to change point
phi_s = Vx/Um; %%
span = (cs_in.r(:,1)-cs_in.r(1,1)) / (cs_in.r(end,1)-cs_in.r(1,1));

plot(phi_s,span);

%% Plotting efficiency chart
load([dr.exp 'N.mat'])
sigmas = N.sigmas;
sigmas(:,1) = [];
etas = var2matrix(eta);
phis = var2matrix(phi);
psis = var2matrix(psi);
rpms = var2matrix(rpm);
figure; hold on; box on; set(gcf, 'color', 'w', 'Position',[350 200 400 420]);
contourf(phis,psis,etas,25,'edgecolor','none'); colorbar;
title("Fan efficiency - CFD");
xlabel("\phi"); ylabel("\psi");

% figure; hold on; box on; set(gcf, 'color', 'w', 'Position',[350 200 400 420]);
% contourf(phis,psis,rpms,25,'edgecolor','none'); colorbar;
% title("RPM - CFD");
% xlabel("\phi"); ylabel("\psi");

% figure; hold on; grid on; box on; axis equal;
% U = cs.bv.omega .* cs.r;
% phi = cs.Vx./U;
% contourf(cs.z,cs.y,cs.Ts,25,'edgecolor','none'); colorbar;

%% Plotting a property
g = ts_read_hdf5([dr.ts 'NRF_design_18_2_avg.hdf5']);

% Cut on r-theta plane
i=50;
c12 = ts_structured_cut(g,12,i,i,1,'en',1,'en');
cs12 = ts_secondary(c12);
c13 = ts_structured_cut(g,13,i,i,1,'en',1,'en');
cs13 = ts_secondary(c13);

figure; hold on; grid on; box on;
title([num2str(i)]);
contourf(cs12.z,cs12.y,(cs12.Vx),25,'edgecolor','none'); colorbar;
contourf(cs13.z,cs13.y,(cs13.Vx),25,'edgecolor','none'); colorbar;
axis equal; %caxis([-10 10]);

% Cut on x-r plane
k=25;
c = ts_structured_cut(g,13,1,'en',1,'en',k,k);
cs = ts_secondary(c);

figure; hold on; grid on; box on;
title([num2str(j)]);
contourf(cs.x,cs.y,(cs.Vx),25,'edgecolor','none'); colorbar;
axis equal; %caxis([-10 10]);

% Cut on x-theta plane
j=10;
c12 = ts_structured_cut(g,12,1,'en',j,j,1,'en');
cs12 = ts_secondary(c12);
c13 = ts_structured_cut(g,13,1,'en',j,j,1,'en');
cs13 = ts_secondary(c13);

figure; hold on; grid on; box on;
title([num2str(j)]);
contourf(cs12.x,-cs12.rt,(cs12.Vx),25,'edgecolor','none'); colorbar; cax12=caxis;
contourf(cs13.x,-cs13.rt,(cs13.Vx),25,'edgecolor','none'); colorbar; cax13=caxis;
caxis([min(cax12(1), cax13(1)) max(cax12(2),cax13(2))]);
contour(cs12.x, -cs12.rt, (cs12.Vx), [0 0], ':r', 'LineWidth',1.5);
contour(cs13.x, -cs13.rt, (cs13.Vx), [0 0], ':r', 'LineWidth',1.5);

axis equal; %caxis([-10 10]);
    
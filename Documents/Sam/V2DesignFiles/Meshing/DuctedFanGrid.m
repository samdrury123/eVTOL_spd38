%   DuctedFanGrid
%                               
%   Generate blades, meshes and CFD solutions for fan designs
%
%   Written by James Taylor & Sam Grimshaw               
%   May 2020

% Clear workspace and load eVTOL environment
clear; close all; clc;
eVTOL

% Script operation
create_geom = 0;
redesign_geom = 0;
run_des = 0;
create_bc = 0;
run_chic = 0;
plot_des = 1;

% Current job variables
job.priority = 'MILLER-SL3-GPU';
job.cards = 1;
job.walltime = 4;

% Design to run
% des_name = 'CRF_3';
des_name = 'CRF_3b';
% des_name = 'NRF_3b';


%% Make blades to specified angle distributions from DuctedFanDesign
if create_geom == 1
    
    % Load blade design from velocity triangles
    load([dr.geom des_name '.mat']);
    
    % Donor blade parameters
    b.sweep = [0 20e-3 20e-3 0]; b.lean = [0 -20e-3 -20e-3 0]; b.dcam_le = -1.1 * ones(1,4); 
    b.dcam_te = -0.4 * ones(1,4); b.qcam = 2 * ones(1,4); b.thick_te = 0.25 * ones(1,4); 
    b.s_thick_max = 0.3 * ones(1,4); b.rad_thick_max = 0.2 * ones(1,4); b.rad_le = 7 * ones(1,4); 
    b.wedge_te = 60 * ones(1,4); b.thick_max = [12 8 8 12] * 1e-3;
    
    % Meridional parameters
    A = pi * (d.g.rcas^2 - d.g.rhub^2); r = 0.5 * (d.g.rcas + d.g.rhub);
    b.A_le = A; b.A_te = A; b.r_le = r; b.r_te = r;
    b.drdx_hub_le = 0; b.drdx_hub_te = 0; b.drdx_cas_le = 0; b.drdx_cas_te = 0;
    
    % Loop over both blades of fan design
    x_ref = [0.025 0.2]; d.g.chi(:,:,1) = - d.g.chi(:,:,1);
    for n = 1:2
        
        % Fit with blade design from DuctedFanDesign
        r_nondim = (d.g.r - d.g.r(1)) / (d.g.r(end) - d.g.r(1));
        b.tchord = spap2(1,4,r_nondim,d.g.c(:,n));
        b.chi_le = spap2(1,4,r_nondim,d.g.chi(:,1,n));
        b.chi_te = spap2(1,4,r_nondim,d.g.chi(:,2,n));
        b.x_ref = x_ref(n);
        b.N = d.g.N(n);
            
        % Fit splines to numeric blade parameters
        B{n} = bl_spline_fit(b);
        
    end
                
    % Record blade speed
    I.rpm = [d.d.rpm1 d.d.rpm2];
    
    % Inlet and outlet coordinates
    B{1}.x_in = -0.1; B{2}.x_out = 0.4;

    % Plot all blades
    figure(); hold on; grid on; box on; axis equal; cols = lines(length(B));
    for n = 1:length(B)
        b = bl_construct_blade(B{n},0);
        surf(b.xrrt(:,:,1),-b.xrrt(:,:,3),b.xrrt(:,:,2),...
            repmat(reshape(cols(n,:),[1 1 3]),[601 43 1]),'linestyle','none');
        plot3(b.xrrt(:,2:2:end,1),-b.xrrt(:,2:2:end,3),b.xrrt(:,2:2:end,2),'k-');
    end
    
    % Save all blades
    save([dr.geom des_name '_splines.mat'],'B','I');
    
end


%% Tweak geometries for better design
if redesign_geom == 1
    
    % Inlet metal angles for NRF design 3
    load([dr.geom 'NRF_3_splines.mat']); B_old = B; r = 0:0.2:1;
    B{1}.chi_le.coefs = B{1}.chi_le.coefs + 2;
    B{2}.chi_le.coefs = B{2}.chi_le.coefs + 0.5;
    chi_te = fnval(B{2}.chi_te,r) - linspace(0.3,2.5,length(r));
    B{2}.chi_te = spap2(B{2}.chi_te.knots,B{2}.chi_te.order,r,chi_te);
    save([dr.geom 'NRF_3b_splines.mat'],'B','I');
    
    % Inlet metal angles for CRF design 3
%     load([dr.geom 'CRF_3_splines.mat']); B_old = B; r = 0:0.2:1;
%     B{1}.chi_le.coefs = B{1}.chi_le.coefs + 1.5;
%     chi_le = fnval(B{2}.chi_le,r) + linspace(1.8,0.8,length(r));
%     B{2}.chi_le = spap2(B{2}.chi_le.knots,B{2}.chi_le.order,r,chi_le);
%     chi_te = fnval(B{2}.chi_te,r) + linspace(5,-1,length(r));
%     B{2}.chi_te = spap2(B{2}.chi_te.knots,B{2}.chi_te.order,r,chi_te);
%     save([dr.geom 'CRF_3b_splines.mat'],'B','I');
    
    % Plot designs for both rows
    cols = lines(2); leg = {'Old','New'};
    for n = 1:length(B)
        h = bl_spline_plot(B_old{n},[],cols(1,:));
        bl_spline_plot(B{n},h,cols(2,:));
        legend(leg);
    end
    
end


%% Create mesh and run TURBOSTREAM on the design
if run_des == 1
    
    % Load design
    load([dr.geom des_name '_splines.mat'])
    
    % Mesh parameters for normal rotating fan
    if strcmp(des_name(1),'N') == 1
        M.stin = [1 0]; M.stout = [-1 0]; M.n_aid = [25 17]; M.n_aou = [25 17];
    else
        M.stin = [1 -1]; M.stout = [0 0]; M.n_bl = [81 49]; M.n_bu = [49 81]; M.n_so = 9;
    end
    
    % Mesh and run the design point in TURBOSTREAM
    nom_prefix = [des_name '_']; job.run_sa = 0;
    job = bl_create_machine(B,dr,nom_prefix,job,I,M,F,D,1,1,1);


end


%% Generate inlet boundary condition and guess files
if create_bc == 1
    
    % Generate inlet boundary condition with 2% high boundary layers
    r = [0 0.02 0.98 1]; Po = [-1.1e3 0 0 -1.1e3] + Pref; ni = 401; r_new = linspace(0,1,ni)';
    A = [r_new Tref * ones(ni,1) interp1(r,Po,r_new,'pchip') zeros(ni,3)];
    dlmwrite([dr.ts 'inlet_2percent.dat'],A,'\t',1,0);
    
end


%% Run inlet profile and characteristic of design
if run_chic == 1
        
    % Define pressures to run
    pouts = (-600:100:600) + Pref;
    
    % Read in previous design
    g = ts_read_hdf5([dr.ts des_name '_avg.hdf5']);

    % Loop over all pressures and submit jobs
    Jobs = cell(0,1);
    for pout = pouts
        for bid = ts_get_bids(g)
            for pid = ts_get_pids(bid,g)
                if g{bid+1}.patch{pid+1}.attribute.kind == 1
                    g{bid+1}.patch{pid+1}.pout = pout;
                end
            end
        end
        job.name = [des_name '_P' num2str(round(pout))];
        Jobs{end+1,1} = ts_sbatch(job,g);
    end
    
    % Monitor all jobs
    ts_smonitor(Jobs);
    
    
end


%% Plot output from fan designs
if plot_des == 1
    
    % Read the case
    g = ts_read_hdf5([dr.ts des_name '_avg.hdf5']);
%     ts_export_paraview(g,[dr.pv des_name '.hdf5'],'HighSpeed')
    ts_export_paraview(g,[dr.pv des_name '.hdf5'],'Track',[],1);
    load([dr.geom des_name '_xr_hub.mat']);
    load([dr.geom des_name '_xr_cas.mat']);
    load([dr.geom des_name '_splines.mat']);
    
    % Cut the rows
    C = ts_cut_rows(g,xr_hub,xr_cas);
    
    % Plot whirl and blade angles at each location
    figure(); r = linspace(0,1,100);
    subplot_tight(1,3,1,[0.1 0.05]); hold on; grid on; box on; title('Inlet');
    xlabel('Whirl Angle / ^\circ'); ylabel('Radius');
    plot(ts_mass_average(C{1}.inlet,'Alpha',2),ts_mass_average(C{1}.inlet,'r_nondim',2));
    plot(-fnval(B{1}.chi_le,r),r);
    plot(ts_mass_average(C{1}.inlet,'Alpha_rel',2),ts_mass_average(C{1}.inlet,'r_nondim',2));
    legend('Absolute Whirl','B1 Metal Inlet','B1 Relative Whirl');
    subplot_tight(1,3,2,[0.1 0.05]); hold on; grid on; box on; title('Betweeen Rows');
    xlabel('Whirl Angle / ^\circ'); set(gca,'yticklabel',[])
    plot(ts_mass_average(C{1}.outlet,'Alpha_rel',2),ts_mass_average(C{1}.outlet,'r_nondim',2));
    plot(-fnval(B{1}.chi_te,r),r);
    plot(ts_mass_average(C{2}.inlet,'Alpha_rel',2),ts_mass_average(C{2}.inlet,'r_nondim',2));
    plot(fnval(B{2}.chi_le,r),r);
    plot(ts_mass_average(C{2}.inlet,'Alpha',2),ts_mass_average(C{2}.inlet,'r_nondim',2));
    legend('B1 Relative Whirl','B1 Metal Outlet','B2 Relative Whirl','B2 Metal Inlet','Absolute Whirl');
    subplot_tight(1,3,3,[0.1 0.05]); hold on; grid on; box on; title('Outlet');
    xlabel('Whirl Angle / ^\circ'); set(gca,'yticklabel',[])
    plot(ts_mass_average(C{2}.outlet,'Alpha_rel',2),ts_mass_average(C{2}.outlet,'r_nondim',2));
    plot(fnval(B{2}.chi_te,r),r);
    plot(ts_mass_average(C{2}.outlet,'Alpha',2),ts_mass_average(C{2}.outlet,'r_nondim',2));
    legend('B2 Relative Whirl','B2 Metal Outlet','Absolute Whirl');    
    
    % Calculate fan 1D performance
    inlet = C{1}.inlet; outlet = C{2}.outlet; rpm = inlet.bv.rpm - outlet.bv.rpm;
    c = ts_calc_chic(inlet,outlet,rpm);
        
    % Calculate thrust and figure of merit
    [Vx,A] = ts_area_average(outlet,'Vx',3); [~,mdot] = ts_mass_average(outlet,'Vx',3);
    mdot = mdot * double(outlet.bv.nblade); A = A * double(outlet.bv.nblade);
    ro = ts_area_average(inlet,'ro',3);
    c.Wx = mdot * inlet.av.cp * (ts_mass_average(outlet,'To',3) - ts_mass_average(inlet,'To',3));
    c.T = mdot * Vx; c.Mf = c.T^1.5 / (c.Wx * (2 * ro * A)^0.5);    
    
    
end





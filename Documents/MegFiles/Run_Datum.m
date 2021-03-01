%    Run_Datum
%    
%    Run TURBOSTREAM on AUTOGRID mesh of datum geometry
%
%    Written by Megan July 2019

% Clear the workspace
close all; clear; clc;

% Script operation
patch_inlet = 1;
run_from_previous = 1;
run_all = 0;
run_single = 0;
run_datum = 0;
run_chic = 0;
extract_bc = 0;

% Load the eVTOL variables
eVTOL

% Boundary conditions
I.pout = 95000; I.rpm = 5000; I.inlet = [dr.ts 'inlet.dat'];
I.poisson_nstep = 10000; I.nstep = 60000;
I.m = 1.4532; I.k = 100;

% File to run
% file_prefix = 'datum_20_';
file_prefix = 'ducted_fan_1_';
run_name = 'diffuser_61_short_splitter_slip_Wall';
    
% Running options
job.priority = 'MILLER-SL2-GPU';
job.cards = 1;
job.walltime = 4;

% Interactive running
% job.node = 'gpu-e-2';

% Update block IDs
% baseline/short/ar/diffuser1/diffuser2/diffuser_21/splittered_diffusers/splittered_inlet
bid_nacelle_out = 22 ; % 24/25/21/23/23/22/18
bid_diffuser_in = 19; % 21/17/17/22/18/22/19/17
bid_in = [24]; %26/26/26/26/26/24/26
bid_out = [26]; %25/24/23/25/25/26/23
bid_farfield = [26]; %25/24/23/25/25/26/23
bid_centre = [24 30 20 23 26]; %[26 27 28 19 30 25]/[26 27 28 30 22 24]/[26 24 28 30 22 23]/[26 28 21 24 25]/[26 28 21 20 25]/[26 28 21 24 25]/[24 30 20 23 26]/[26 24 31 28 21 23]/;
bid_nose = 3; %11/11/11/4/4/4/3/11/14
bid_tail = 5; %7/8/4/6/6/6/5/4/1


%% Initial guess update

if run_from_previous == 1
    I.guess = [dr.ts 'ducted_fan_1_diffuser_54_splitte.hdf5'];
    
else
    I.guess = [1.2 40 0 0 2e5 0.0002];
    
end 

%% Patch the fan to the inlet and exit pointwise mesh
if patch_inlet == 1
    
    % Read both meshes g is Autogrid and g_inlet is Pointwise
    g = ts_read_hdf5([dr.ag 'datum_22_input.hdf5']);
    g_inlet = ts_read_hdf5([dr.pw 'input_diffuser_61_short_splitter.hdf5']);
    
    % Load hub and casing lines
    xr_hub = dlmread([dr.ag 'xr_hub.dat'],'\t',3,0);
    xr_cas = dlmread([dr.ag 'xr_cas.dat'],'\t',3,0);
    xr_stub = dlmread([dr.ag 'g2.dat'],'\t',3,0);
    xr_hub = xr_hub(:,[2 1]); xr_cas = xr_cas(:,[2 1]); xr_stub = xr_stub(:,[2 1]); 
    xr_rotating = xr_stub(1:2,:);
    xr_stationary = xr_stub(end-1:end,:);
    
    % Set rpms on blade rows and then cavities
    g = ts_set_rpms(g,[1 0],xr_hub,xr_cas,xr_stationary,xr_rotating);

    % Set stationary faces at bottom of shroud wells
    for bid = ts_get_bids(g)
        x = sum(sum(sum(g{bid+1}.x))) / numel(g{bid+1}.x);
        r = sum(sum(sum(g{bid+1}.r))) / numel(g{bid+1}.r);
        if r < interp1(xr_hub(:,1),xr_hub(:,2),x)
            g{bid+1}.bv.rpmj2 = 0;
        end
    end
    
    % Make the nose bullet spin
    g_inlet{bid_nose+1}.bv.rpmi2 = 1; g_inlet{bid_tail+1}.bv.rpmj1 = 1; 
    
    % Add periodics to all pointwise blocks
    for bid = ts_get_bids(g_inlet)
        ni = g_inlet{bid+1}.attribute.ni; nj = g_inlet{bid+1}.attribute.nj; 
        nk = g_inlet{bid+1}.attribute.nk;
        g_inlet = ts_add_periodic_patch(bid,bid,[0 ni ; 0 nj ; 0 1],[0 ni ; 0 nj ; nk-1 nk],...
            [0 1 2],[0 1 2],g_inlet,5);
        g{bid+1}.av.turbvis_lim = 100000000;

    end
  
    % Stitch both meshes together
    g = ts_append(g,g_inlet);

    % Display rotating surfaces
    ts_check_rpms(g);
   
    % Block and patch ids to add the mixing planes
    bid_rot_in = 2; pid_rot_in = 4; 
    bid_sta_out = 9; pid_sta_out = 1; 

    % Add the patch at the pointwise inlet
    p.bid = bid_nacelle_out; p.pid = length(g{p.bid+1}.patch); p.kind = 2; 
    p.ist = g{p.bid+1}.attribute.ni-1; p.ien = g{p.bid+1}.attribute.ni; 
    p.jst = 0; p.jen = g{p.bid+1}.attribute.nj; 
    p.kst = 0; p.ken = g{p.bid+1}.attribute.nk; 
    p.idir = 0; p.jdir = 1; p.kdir = 2;
    p.nxbid = bid_rot_in; p.nxpid = pid_rot_in;
    p.nface = 0; p.nt = 1;
    g{p.bid+1}.patch{end+1}.attribute = p; g{p.bid+1}.attribute.np = g{p.bid+1}.attribute.np + 1;
    
    % Convert the inlet patch to a mixing plane
    g{bid_rot_in+1}.patch{pid_rot_in+1}.attribute.kind = 2;
    g{bid_rot_in+1}.patch{pid_rot_in+1}.attribute.nxbid = p.bid;
    g{bid_rot_in+1}.patch{pid_rot_in+1}.attribute.nxpid = p.pid;
    
    % Add the patch at the pointwise inlet
    p.bid = bid_diffuser_in; p.pid = length(g{p.bid+1}.patch); p.kind = 2; 
    p.ist = 0; p.ien = 1; 
    p.jst = 0; p.jen = g{p.bid+1}.attribute.nj; 
    p.kst = 0; p.ken = g{p.bid+1}.attribute.nk; 
    p.idir = 0; p.jdir = 1; p.kdir = 2;
    p.nxbid = bid_sta_out; p.nxpid = pid_sta_out;
    p.nface = 0; p.nt = 1;
    g{p.bid+1}.patch{end+1}.attribute = p; g{p.bid+1}.attribute.np = g{p.bid+1}.attribute.np + 1;
    
    % Convert the inlet patch to a mixing plane
    g{bid_sta_out+1}.patch{pid_sta_out+1}.attribute.kind = 2;
    g{bid_sta_out+1}.patch{pid_sta_out+1}.attribute.nxbid = p.bid;
    g{bid_sta_out+1}.patch{pid_sta_out+1}.attribute.nxpid = p.pid;   
    
    % Add the new inlet patches
    for bid = bid_in
        
        % Create the patch
        p.bid = bid; p.pid = length(g{p.bid+1}.patch); p.kind = 0; 
        p.ist = 0; p.ien = 1; 
        p.jst = 0; p.jen = g{p.bid+1}.attribute.nj; 
        p.kst = 0; p.ken = g{p.bid+1}.attribute.nk; 
        p.idir = 0; p.jdir = 1; p.kdir = 2;
        p.nxbid = -1; p.nxpid = -1; p.nface = 0; p.nt = 1;
        g{bid+1}.patch{end+1}.attribute = p; g{p.bid+1}.attribute.np = g{p.bid+1}.attribute.np + 1;
        
         % Add periodicity to make full 360 degree
        g{bid+1}.bv.nblade = 36; g{bid+1}.bv.fblade = 36;
     
    end
 
%     Make farfield patch an inlet patch
%     for bid = bid_farfield
%         
%         Create the patch
%         p.bid = bid; p.pid = length(g{p.bid+1}.patch); p.kind = 0; 
%         p.ist = 0; p.ien = g{p.bid+1}.attribute.ni; 
%         p.jst = g{p.bid+1}.attribute.nj-1; p.jen = g{p.bid+1}.attribute.nj; 
%         p.kst = 0; p.ken = g{p.bid+1}.attribute.nk; 
%         p.idir = 0; p.jdir = 1; p.kdir = 2;
%         p.nxbid = -1; p.nxpid = -1; p.nface = 0; p.nt = 1;
%         g{bid+1}.patch{end+1}.attribute = p; g{p.bid+1}.attribute.np = g{p.bid+1}.attribute.np + 1;
%         
%          Add periodicity to make full 360 degree
%         g{bid+1}.bv.nblade = 36; g{bid+1}.bv.fblade = 36;
%      
%     end
      
     % Add patch to farfield as slip wall
    for bid = bid_farfield
        % Create the patch
        p.bid = bid; p.pid = length(g{p.bid+1}.patch); p.kind = 7; 
        p.ist = 0; p.ien = g{p.bid+1}.attribute.ni; 
        p.jst = g{p.bid+1}.attribute.nj-1; p.jen = g{p.bid+1}.attribute.nj;  
        p.kst = 0; p.ken = g{p.bid+1}.attribute.nk; 
        p.idir = 0; p.jdir = 1; p.kdir = 2;
        p.nxbid = -1; p.nxpid = -1; p.nface = 0; p.nt = 1;
        g{bid+1}.patch{end+1}.attribute = p; g{p.bid+1}.attribute.np = g{p.bid+1}.attribute.np + 1;
        
    end
    
    
    % Add the new outlet patches
    for bid = bid_out
        
        % Create the patch
        p.bid = bid; p.pid = length(g{p.bid+1}.patch); p.kind = 1; 
        p.ist = g{p.bid+1}.attribute.ni-1; p.ien = g{p.bid+1}.attribute.ni; 
        p.jst = 0; p.jen = g{p.bid+1}.attribute.nj; 
        p.kst = 0; p.ken = g{p.bid+1}.attribute.nk; 
        p.idir = 0; p.jdir = 1; p.kdir = 2;
        p.nxbid = -1; p.nxpid = -1; p.nface = 0; p.nt = 1;
        g{bid+1}.patch{end+1}.attribute = p; g{p.bid+1}.attribute.np = g{p.bid+1}.attribute.np + 1;
        
        % Set the boundary conditions
        g{bid+1}.patch{end}.pout = Pref;
        g{bid+1}.patch{end}.throttle_type = 0;
        g{bid+1}.patch{end}.ipout = -3;
    end
     
%      % Change upper boundary to farfield
%      for bid = bid_farfield
%         
%         % Create the patch
%         p.bid = bid; p.pid = length(g{p.bid+1}.patch); p.kind = 13; 
%         p.ist = 0; p.ien = g{p.bid+1}.attribute.ni; 
%         p.jst = g{p.bid+1}.attribute.nj-1; p.jen = g{p.bid+1}.attribute.nj; 
%         p.kst = 0; p.ken = g{p.bid+1}.attribute.nk; 
%         p.idir = 0; p.jdir = 1; p.kdir = 2;
%         p.nxbid = -1; p.nxpid = -1; p.nface = 0; p.nt = 1;
%         g{bid+1}.patch{end+1}.attribute = p; g{p.bid+1}.attribute.np = g{p.bid+1}.attribute.np + 1;
%      
%         % Add periodicity to make full 360 degree
%         g{bid+1}.bv.nblade = 36; g{bid+1}.bv.fblade = 36;
%         
%         % Set the boundary conditions
%         g{bid+1}.patch{end}.pstag_free = Pref;
%         g{bid+1}.patch{end}.pstat_free = Pref-5;
%         g{bid+1}.patch{end}.sf_free = 0.5;
%         g{bid+1}.patch{end}.tstag_free = Tref;
%         
%         % Adjust mixing length limit
%         g{bid+1}.xllim = 10000; 
%         
%      end
%     
%    % Change inlet to farfield
%      for bid = bid_in
%         
%         % Create the patch
%         p.bid = bid; p.pid = length(g{p.bid+1}.patch); p.kind = 13; 
%         p.ist = 0; p.ien = 1; 
%         p.jst = 0; p.jen = g{p.bid+1}.attribute.nj; 
%         p.kst = 0; p.ken = g{p.bid+1}.attribute.nk; 
%         p.idir = 0; p.jdir = 1; p.kdir = 2;
%         p.nxbid = -1; p.nxpid = -1; p.nface = 0; p.nt = 1;
%         g{bid+1}.patch{end+1}.attribute = p; g{p.bid+1}.attribute.np = g{p.bid+1}.attribute.np + 1;
%      
%         % Add periodicity to make full 360 degree
%         g{bid+1}.bv.nblade = 36; g{bid+1}.bv.fblade = 36;
%         
%         % Set the boundary conditions
%         g{bid+1}.patch{end}.pstag_free = Pref;
%         g{bid+1}.patch{end}.pstat_free = Pref;
%         g{bid+1}.patch{end}.sf_free = 0.5;
%         g{bid+1}.patch{end}.tstag_free = Tref;
%         
%         % Adjust mixing length limit
%         g{bid+1}.xllim = 10000; 
%         
%     end
    
     
    % Add patch to centre line boundary 
    for bid = bid_centre
        
        % Create the patch
        p.bid = bid; p.pid = length(g{p.bid+1}.patch); p.kind = 7; 
        p.ist = 0; p.ien = g{p.bid+1}.attribute.ni; 
        p.jst = 0; p.jen = 1; 
        p.kst = 0; p.ken = g{p.bid+1}.attribute.nk; 
        p.idir = 0; p.jdir = 1; p.kdir = 2;
        p.nxbid = -1; p.nxpid = -1; p.nface = 0; p.nt = 1;
        g{bid+1}.patch{end+1}.attribute = p; g{p.bid+1}.attribute.np = g{p.bid+1}.attribute.np + 1;
        
    end
     
    % Save patch data to file
    ts_check_patches(g,2,[dr.ts file_prefix 'patch.dat'])
    
    % Apply inlet boundary conditions and rotational speed
    g = ts_apply_inlet(g,[Pref Tref 0 0], 0, -max(I.rpm));

    % Apply initial settings
    g = ts_create_initial(g,10000,10000,I.guess);
    g{1}.av.nomatch_int = 0;

    % Write yplus values
    g{1}.av.write_yplus = 1;
    
    % Change averaging limit
    g{1}.av.nstep_save_start = 200000;
    
    % Set initial velocity in cavities and slots to zero
    for bid = ts_get_bids(g);
        x = mean(g{bid+1}.x(:)); r = mean(g{bid+1}.r(:));
        if r < interp1(xr_hub(:,1),xr_hub(:,2),x)
            g{bid+1}.rovx = zeros(size(g{bid+1}.rovx));
            g{bid+1}.rovr = zeros(size(g{bid+1}.rovx));
            g{bid+1}.rorvt = g{bid+1}.rorvt + 5;
        end
    end
    
     % Change rfin to 0
     for bid = [24]
         g{bid+1}.patch{8}.rfin = 0.3;
     end
     
%      for bid = [26]
%          g{bid+1}.patch{7}.rfin = 0.3;
%      end
%      
     
% %     % Set pitch angle at inlet
    for bid = 24
        
        for t = 1:g{bid+1}.attribute.nk
            alpha = transpose(linspace(0, -45, g{bid+1}.attribute.nj));  
            a = 1 + (g{bid+1}.attribute.nj*(t-1)); b = g{bid+1}.attribute.nj + (g{bid+1}.attribute.nj*(t-1));
            g{bid+1}.patch{8}.pitch(a:b) = alpha();
        end
     
    end
%     
%     
%     % Set pitch angle at farfield when set as inlet
%     for bid = 26
%         
%         for t = 1:g{bid+1}.attribute.nk
%             alpha = transpose(linspace(-45, -45, g{bid+1}.attribute.nj));
%             a = 1 + (g{bid+1}.attribute.nj*(t-1)); b = g{bid+1}.attribute.nj + (g{bid+1}.attribute.nj*(t-1));
%             g{bid+1}.patch{7}.pitch(a:b) = alpha();
%         end
%      
%     end
%     
     processed_name = [file_prefix run_name '_preprocess'];
     ts_write_hdf5(g,[dr.ts processed_name '.hdf5']);

end 

%% Run all cases 
if run_all == 1
    
    % Run poisson loop and coarse mixing length
    job.name = [file_prefix run_name 'poisson'];
    job = ts_smonitor(ts_sbatch(job,g));

%     % Force run of ML
%     job.finished = 1; job.name = [file_prefix run_name 'ML']; job.logname = ['log_' file_prefix run_name 'ML.txt'];

    % Run fine mixing length solution
     if isnan(job.finished) == 0
         g = ts_read_hdf5([job.directory job.name '.hdf5']);
         g = ts_smooth_poisson(g);
%          g = ts_create_restart(g,15000,1);
         g = ts_create_restart(g,400000,1);
         job.name = [file_prefix run_name 'ML'];
         job = ts_smonitor(ts_sbatch(job,g)); 
     end

   
    % Force run of SA only
    %job.finished = 1; job.name = [file_prefix run_name 'ML']; job.logname = ['log_' file_prefix run_name 'ML.txt'];
     
%     % Run SA restart
%     if isnan(job.finished) == 0
%         g = ts_read_hdf5([job.directory job.name '.hdf5']);
%         g = ts_create_restart(g,I.nstep,2);
%         %g{1}.av.poisson_nstep = 10000;
%         job.name = [file_prefix run_name(1:end-1)];
%         job = ts_smonitor(ts_sbatch(job,g));
%     end
     
    % Plot convergence history
    ts_plot_conv(job);
  
end

%% Run a single case with new initial guess
if run_single == 1
    
    % Read the pre-processed mesh
    g = ts_read_hdf5([dr.ts processed_name '.hdf5']);
    
    % Settings to run SA/ML and poisson a little bit more
    g = ts_create_initial(g,10000,10000,I.guess);
    g{1}.av.nomatch_int = 0;
    g{1}.av.write_yplus = 1;
    g{1}.av.nstep_save_start = 20000;
    g = ts_create_restart(g,400000,1); %change end entry to 2 if SA
    g{1}.av.poisson_nstep = 10000;
    
    % Run the case
    job.name = [file_prefix run_name(1:end-1)];
    job = ts_smonitor(ts_sbatch(job,g));
    
    
end

%% Pre-process the TURBOSTREAM input files
if run_datum == 1

% Read in converted file
g = ts_read_hdf5([dr.ag file_prefix 'input.hdf5']);

% Load hub and casing lines
xr_hub = dlmread([dr.ag 'xr_hub.dat'],'\t',3,0);
xr_cas = dlmread([dr.ag 'xr_cas.dat'],'\t',3,0);
xr_stub = dlmread([dr.ag 'g2.dat'],'\t',3,0);
xr_hub = xr_hub(:,[2 1]); xr_cas = xr_cas(:,[2 1]); xr_stub = xr_stub(:,[2 1]); 
xr_rotating = xr_stub(1:2,:);
xr_stationary = xr_stub(end-1:end,:);

% Set exit boundary conditions
for bid = ts_get_bids(g)
    for pid = ts_get_pids(bid,g)
        if g{bid+1}.patch{pid+1}.attribute.kind == 1
            
            % Decide on exit boundary condition type
            if isfield(I,'m') == 0
                
                % Constant back pressure
                g{bid+1}.patch{pid+1}.pout = I.pout(1);
                g{bid+1}.patch{pid+1}.throttle_type = 0;
            else
                
                % PID throttle control
                g{bid+1}.patch{pid+1}.pout = I.pout(1);
                g{bid+1}.patch{pid+1}.throttle_type = 1;
                g{bid+1}.patch{pid+1}.throttle_target = I.m;
                g{bid+1}.patch{pid+1}.throttle_k0 = I.k;
                g{bid+1}.patch{pid+1}.throttle_k1 = I.k;
                g{bid+1}.patch{pid+1}.throttle_k2 = I.k;
            end
            
            % Set radial equlibrium control
            if numel(I.pout) == 2
                g{bid+1}.patch{pid+1}.ipout = I.pout(2); 
            else
                g{bid+1}.patch{pid+1}.ipout = 3;
            end
            
            % Record exit patch location
            bid_out = bid; pid_out = pid;
        end
    end
end

% Save patch data to file
ts_check_patches(g,2,[dr.ts file_prefix 'patch.dat'])

% Set rpms on blade rows and then cavities
g = ts_set_rpms(g,[1 0],xr_hub,xr_cas,xr_stationary,xr_rotating);

% Set stationary faces at bottom of shroud wells
for bid = ts_get_bids(g)
    x = sum(sum(sum(g{bid+1}.x))) / numel(g{bid+1}.x);
    r = sum(sum(sum(g{bid+1}.r))) / numel(g{bid+1}.r);
    if r < interp1(xr_hub(:,1),xr_hub(:,2),x)
        g{bid+1}.bv.rpmj1 = 0;
    end
end
    
% Display rotating surfaces
ts_check_rpms(g);

% Apply inlet boundary conditions and rotational speed
g = ts_apply_inlet(g,I.inlet, 0, -max(I.rpm));

% Apply initial settings
g = ts_create_initial(g,15000,I.poisson_nstep,I.guess);
g{1}.av.nomatch_int = 0;

% Write yplus values
g{1}.av.write_yplus = 1;

% Set initial velocity in cavities and slots to zero
for bid = ts_get_bids(g);
    x = mean(g{bid+1}.x(:)); r = mean(g{bid+1}.r(:));
    if r < interp1(xr_hub(:,1),xr_hub(:,2),x) || r > interp1(xr_cas(:,1),xr_cas(:,2),x)
        g{bid+1}.rovx = zeros(size(g{bid+1}.rovx));
        g{bid+1}.rovr = zeros(size(g{bid+1}.rovx));
        g{bid+1}.rorvt = g{bid+1}.rorvt + 5;
    end
end


%% Run file in TURBOSTREAM on Darwin
    
% Run poisson loop and coarse mixing length
% job.name = [file_prefix 'poisson'];
% job = ts_smonitor(ts_sbatch(job,g));
% 
% % Run fine mixing length solution
% if isnan(job.finished) == 0
%     g = ts_read_hdf5([job.directory job.name '.hdf5']);
%     g = ts_smooth_poisson(g);
%     g = ts_create_restart(g,15000,1);
% 
%     if isfield(I,'m') == 1
%         [~,~,dpout] = ts_plot_conv([job.directory job.logname],0);
%         g{bid_out+1}.patch{pid_out+1}.pout = g{bid_out+1}.patch{pid_out+1}.pout + dpout;
%     end
% 
%     job.name = [file_prefix 'ML'];
%     job = ts_smonitor(ts_sbatch(job,g)); 
% end

% Force run of SA only
job.finished = 1; job.name = [file_prefix 'ML']; job.logname = ['log_' file_prefix 'ML.txt'];

% Run SA restart
if isnan(job.finished) == 0
    g = ts_read_hdf5([job.directory job.name '.hdf5']);
    g = ts_create_restart(g,I.nstep,2);
     
    if isfield(I,'m') == 1
        [~,~,dpout] = ts_plot_conv([job.directory job.logname],0);
        g{bid_out+1}.patch{pid_out+1}.pout = g{bid_out+1}.patch{pid_out+1}.pout + dpout;
    end

    job.name = file_prefix(1:end-1);
    job = ts_smonitor(ts_sbatch(job,g));
    
end

% Plot convergence history
ts_plot_conv(job);

end


%% Run the characteristic
if run_chic == 1

% Specify start file in job structure
job.name = file_prefix(1:end-1);

% Coordinates of inlet and outlet planes
coords{1}.xr = [0.0066 0 ; 0.0066 1]; coords{1}.bids = 2;
coords{2}.xr = [0.105 0 ; 0.105 1]; coords{2}.bids = 9;

% Run the chic
job = ts_run_chic(job,coords,1);

end

%% Etract guess file for use in future runs
if extract_bc == 1 
    
    % Read in design point from datum geometry
    g = ts_read_hdf5([dr.ts 'datum_20_avg.hdf5']);
    
    %Find bounding x coordinates
    x_min = inf; x_max = -inf;
    for bid = ts_get_bids(g)
        x_min = min(min(min(min(g{bid+1}.x))), x_min);
        x_max = max(max(max(max(g{bid+1}.x))), x_max);
    end
    x = linspace(x_min,x_max,100);
    
    % Take cuts through domain and cut ld averages
    guess = zeros(length(x)-2,8); varnames = {'ro' 'rovx' 'rovr' 'rorvt' 'roe' 'trans_dyn_vis' 'phi'};
    for i = 2:length(x)-1
        c = ts_secondary(ts_unstructured_cut(g,x(i))); guess(i-1,1) = x(i);
        for v = 1:length(varnames)
            guess(i-1,v+1) = ts_mass_average(c,varnames{v},3);
        end
    end
    
    % Write 1D Guess
    dlmwrite([dr.ts 'Guess_p1.dat'],guess)
   
end


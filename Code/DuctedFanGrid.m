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
create_geom = 1;
redesign_geom = 0;
mesh_des = 0;
run_des = 0;
create_bc = 0;
run_chics = 0;
plot_des = 0;
plot_dev = 0;
plot_chics = 0;
write_ibl = 0;
pressure_load = 0;
plot_hubrat = 0;
plot_diff = 0;
plot_speed = 0;

% Current job variables
job.priority = 'MILLER-SL3-GPU';
job.cards = 1;
job.walltime = 4;

% Design to run
% des_name = 'CRF_3';
des_name = 'CRF_3b';
% des_name = 'NRF_3b';

% Load test matrix
test_name = 'NRF_SPD_';
design_name = 'NRF_design_';
% test_name = 'prelim';
% test_name = 'hubrat';
% test_name = 'diffusion';
% test_name = 'cube_2';
% load([dr.geom test_name '.mat'])
% load([dr.geom 'NRF_design_1.mat'])


%% Make blades to specified angle distributions from DuctedFanDesign
if create_geom == 1
    for i=1:25
        load([dr.geom test_name     num2str(i) '.mat'])
        load([dr.geom design_name   num2str(i) '.mat'])
        des_name = [design_name num2str(i)];
        
    % Duct sizes and gapping
    dx_gap = 70e-3; dx_duct = 0.3; x_in = 0; chord_min = 0.08;
    
    % Loop over all designs in test matrix
    for o = 1:numel(phi)
    
        % Load blade design from velocity triangles
%         des_name = ['CRF2_P' num2str(phi(o)*100) 'S' num2str(sig(o)*100) 'R' num2str(r_cas(o)*100)];
%         des_name = ['CRF2_P' num2str(phi(o)*100) 'S' num2str(sig(o)*100) 'R' num2str(r_cas(o)*100) ...
%             'DF' num2str(df(o)*100) 'HT' num2str(ht(o)*100)];
%         load([dr.geom des_name '.mat']);

        % Donor blade parameters
        b.sweep = [0 20e-3 20e-3 0]; b.lean = [0 -20e-3 -20e-3 0]; b.dcam_le = -1.1 * ones(1,4); 
        b.dcam_te = -0.4 * ones(1,4); b.qcam = 2 * ones(1,4); b.thick_te = 0.25 * ones(1,4); 
        b.s_thick_max = 0.3 * ones(1,4); b.rad_thick_max = 0.2 * ones(1,4); b.rad_le = 7 * ones(1,4); 
        b.wedge_te = 60 * ones(1,4); b.thick_max = [12 8 8 12] * 1e-3; 

        % Meridional parameters
        A = pi * (d.g.rcas^2 - d.g.rhub^2); r = 0.5 * (d.g.rcas + d.g.rhub);
        b.A_le = A; b.A_te = A; b.r_le = r; b.r_te = r; b.x_ref = 0;
        b.drdx_hub_le = 0; b.drdx_hub_te = 0; b.drdx_cas_le = 0; b.drdx_cas_te = 0;
        
        % Loop over both blades of fan design
        d.g.chi(:,:,1) = - d.g.chi(:,:,1);
        for n = 1:2
            
            % Set minimum chord
            d.g.c(:,n) = max(d.g.c(:,n),chord_min);

            % Fit with blade design from DuctedFanDesign
            r_nondim = (d.g.r - d.g.r(1)) / (d.g.r(end) - d.g.r(1));
            b.tchord = spap2(1,4,r_nondim,d.g.c(:,n));
            b.chi_le = spap2(1,4,r_nondim,d.g.chi(:,1,n));
            b.chi_te = spap2(1,4,r_nondim,d.g.chi(:,2,n));
            b.N = d.g.N(n);

            % Fit splines to numeric blade parameters
            B{n} = bl_spline_fit(b);

        end

        % Record blade speed
        I.rpm = [d.d.rpm1 d.d.rpm2];
        
        % Get leading and trailing edge coordinates
        for n = 1:length(B)
            B{n} = bl_construct_blade(B{n},0);
            i_te = round(size(B{n}.xrrt,1) / 2);
            x_le(:,n) = B{n}.xrrt(1,:,1); x_te(:,n) = B{n}.xrrt(i_te,:,1);
        end
        
        % Set first blade axial position from inlet plane
        dx_1 = -min(x_le(:,1)) + dx_duct; B{1}.x_ref = dx_1;
        
        % Set second blade axial position from gapping
        dx_2 = dx_gap - min(x_le(:,2) - x_te(:,1));
        B{2}.x_ref = dx_1 + dx_2;
        
        % Inlet and outlet coordinates
        B{1}.x_in = 0; B{2}.x_out = dx_1 + dx_2 + max(x_te(:,2)) + dx_duct;
        
        % Recreate geometry
        for n = 1:length(B); B{n} = bl_construct_blade(B{n},0); end;

        % Plot all blades
        figure(); hold on; grid on; box on; axis equal; cols = lines(length(B));
        for n = 1:length(B);
            plot3(B{n}.xrrt(:,2:2:end,1),-B{n}.xrrt(:,2:2:end,3),B{n}.xrrt(:,2:2:end,2),'k-');
            surf(B{n}.xrrt(:,:,1),-B{n}.xrrt(:,:,3),B{n}.xrrt(:,:,2),...
                repmat(reshape(cols(n,:),[1 1 3]),[601 43 1]),'linestyle','none');
        end
      
        % Save all blades
        save([dr.geom des_name '_splines.mat'],'B','I');
      
        % Store blade counts
        N(o,:) = [B{1}.N B{2}.N];
        
    end
    end
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


%% Create mesh for the design with minimum blade numbers
if mesh_des == 1
    
    % Loop over all designs in test matrix
    Jobs = cell(0,1);
    for o = 1:numel(phi)
    
        % Load design
        if exist('df','var') == 0
            des_name = ['CRF2_P' num2str(phi(o)*100) 'S' num2str(sig(o)*100) 'R' num2str(r_cas(o)*100)];
        elseif exist('ht','var') == 0
            des_name = ['CRF2_P' num2str(phi(o)*100) 'S' num2str(sig(o)*100) ...
                'R' num2str(r_cas(o)*100) 'DF' num2str(df(o) * 100)];
        else
            des_name = ['CRF2_P' num2str(phi(o)*100) 'S' num2str(sig(o)*100) ...
                'R' num2str(r_cas(o)*100) 'DF' num2str(df(o) * 100) 'HT' num2str(ht(o) * 100)];
        end
        load([dr.geom des_name '_splines.mat'])

        % Mesh parameters for normal or counter rotating fans
        if strcmp(des_name(1),'N') == 1
            M.stin = [1 0]; M.stout = [-1 0]; M.n_aid = [25 17]; M.n_aou = [25 17];
        else
            M.stin = [1 -1]; M.stout = [0 0]; M.n_bl = [81 49]; M.n_bu = [49 81]; M.n_so = 9;
        end
        
        % Turn off fillets
        F{1}.flag = [0 0]; F{2}.flag = [0 0]; 
        
        % Set minimum blade numbers for meshing
        for n = 1:length(B); if B{n}.N < 12; B{n}.N = 12; end; end;
        I.nstep = 10; I.poisson_nstep = 10;
%         B{2}.N = 32;
        % Mesh and run the design point in TURBOSTREAM
        if isempty(dir([dr.ts des_name '_avg.*'])) == 1
            nom_prefix = [des_name '_']; job.run_sa = 1;
            try
                Jobs{end+1,1} = bl_create_machine(B,dr,nom_prefix,job,I,M,F,D,1,1,1);
            end
        end
        
    end
    
    % Monitor jobs
    ts_smonitor(Jobs);

end


%% Run meshed cases in TURBOSTREAM
if run_des == 1
    
    % Mass flow increase to meet cruise conditions
    mass_fac = 1.2;
    
    % Block numbers for pitch extension
    bids{1} = [0 2 5]; bids{2} = [7 9 13];
    ks{1} = [53 53 17]; ks{2} = [1 1 17];
    bids_row{1} = 0:6; bids_row{2} = 7:13;
    
    % Pitchwise expansion ratio
    rk = 1.4; nk_rng = 5:100;
    
    % Process initial guess
%     g = ts_read_hdf5([dr.ts 'CRF_3b_avg.hdf5']);
%     g = ts_read_hdf5([dr.ts 'guess.hdf5']);
    g = ts_read_hdf5([dr.ts 'CRF2_P80S80R62DF38HT30_1_avg.hdf5']);
    guessnames = {'ro' 'rovx' 'rovr' 'rorvt' 'roe'};
    for v = 1:length(guessnames)
        for bid = ts_get_bids(g)
%             guess.(guessnames{v})(bid+1) = mean(g{bid+1}.(guessnames{v})(:));
            guess.(guessnames{v}){bid+1} = mean(g{bid+1}.(guessnames{v}),3);
%             guess.(guessnames{v}){bid+1} = g{bid+1}.(guessnames{v});
        end
    end
    
    % Loop over all designs in test matrix
    Jobs = cell(0,1);
    for o = 1:numel(phi)
        
        % Job name
        if exist('df','var') == 0
            des_name = ['CRF2_P' num2str(phi(o)*100) 'S' num2str(sig(o)*100) 'R' num2str(r_cas(o)*100)];
        elseif exist('ht','var') == 0
            des_name = ['CRF2_P' num2str(phi(o)*100) 'S' num2str(sig(o)*100) ...
                'R' num2str(r_cas(o)*100) 'DF' num2str(df(o) * 100)];
        else
            des_name = ['CRF2_P' num2str(phi(o)*100) 'S' num2str(sig(o)*100) ...
                'R' num2str(r_cas(o)*100) 'DF' num2str(df(o) * 100) 'HT' num2str(ht(o) * 100)];
        end
        
        % Set job name based on speed ratio
%         job.name = [des_name '_1']; rw = 0;
%         job.name = [des_name '_2']; rw = 0.15;
        job.name = [des_name '_3']; rw = -0.15;
        
        % Check if the job has already been run
        if isempty(dir([dr.ts job.name '_avg.*'])) == 1;
    
        % Load blade design from velocity triangles
        load([dr.geom des_name '_splines.mat'])
        g = ts_read_hdf5([dr.ts des_name '_avg.hdf5']);
        
        % Loop over both rows and correct pitch with extra points
        for n = 1:length(bids)
            if B{n}.N < double(g{bids{n}(1)+1}.bv.nblade)

                % Pitch difference between desired and meshed blade count
                dt = 2 * pi / B{n}.N - 2 * pi / double(g{bids{n}(1)+1}.bv.nblade);

                % Average cell size at boundary
                dt_cell = abs(diff(g{bids{n}(end)+1}.rt(:,:,end-1:end) ./ ...
                    g{bids{n}(end)+1}.r(:,:,end-1:end),1,3));
                dt_av = 0.5 * (max(dt_cell(:)) - min(dt_cell(:)));

                % Generate pitch coordinates with desired expansion ratio
                if dt_av < 0.4 * dt
                    rk_list = zeros(size(nk_rng)); t_list = cell(size(nk_rng));
                    for k = 1:length(nk_rng)
                        [t_list{k},rk_list(k)] = geo_bunch(nk_rng(k),0,dt,dt_av,dt_av,0.5);
                    end
                    k = find(rk_list < rk,1); t = t_list{k}; nk = nk_rng(k);
                else
                    t = linspace(0,dt,3); nk = 3;
                end
                
                % Record pitch coordinates in the right direction
                if n == 2; t = -t; end;
                t = reshape(t,[1 1 nk]);
                
                % Increase pitch of blocks
                for m = 1:length(bids{n})
                    
                    % Block id and index
                    bid = bids{n}(m); k = ks{n}(m);
                    
                    % Record block pitchwise coordinates
                    t_bk = g{bid+1}.rt ./ g{bid+1}.r;
                    
                    % Repeat other coordinates and flow
                    varnames = fieldnames(g{bid+1});
                    for v = 1:length(varnames)
                        if isnumeric(g{bid+1}.(varnames{v})) == 1 && strcmp(varnames{v},'rt') == 0
                            if k ~= 1
                                g{bid+1}.(varnames{v}) = cat(3,g{bid+1}.(varnames{v}),...
                                    repmat(g{bid+1}.(varnames{v})(:,:,k),[1 1 nk-1]));
                            else
                                g{bid+1}.(varnames{v}) = cat(3,repmat(g{bid+1}.(varnames{v})(:,:,k),...
                                    [1 1 nk-1]),g{bid+1}.(varnames{v}));                                
                            end
                        end
                    end
                    
                    % Pitchwise coordinates
                    ni = g{bid+1}.attribute.ni; nj = g{bid+1}.attribute.nj;
                    if k ~= 1
                        t_bk = cat(3,t_bk,repmat(t(2:end),[ni nj 1]) + repmat(t_bk(:,:,k),[1 1 nk-1]));
                        k_smth = k-3:size(t_bk,3);
                    else
                        t_bk = cat(3,repmat(flip(t(2:end),3),[ni nj 1]) + ...
                            repmat(t_bk(:,:,k),[1 1 nk-1]),t_bk);
                        k_smth = k:nk+3;
                    end
                    
                    % Smooth pitchwise coordinates
                    t_bk(:,:,k_smth(2:end-1)) = (t_bk(:,:,k_smth(1:end-2)) + ...
                        t_bk(:,:,k_smth(2:end-1)) + t_bk(:,:,k_smth(3:end))) / 3;
                    g{bid+1}.rt = t_bk .* g{bid+1}.r;
                    
                    % Coordinate list to update patch extents
                    nk_old = g{bid+1}.attribute.nk; nk_new = g{bid+1}.attribute.nk + nk - 1;
                    k_old = 0:nk_old; 
                    if k ~= 1; k_new = k_old; else; k_new = nk_new - nk_old:nk_new; end;
                    
                    % Increase number of points in attributes
                    g{bid+1}.attribute.nk = nk_new;
                    
                    % Update patch extents
                    for pid = ts_get_pids(g,bid)
                        if g{bid+1}.patch{pid+1}.attribute.kst == nk_old-1
                            g{bid+1}.patch{pid+1}.attribute.kst = nk_new - 1;
                        elseif g{bid+1}.patch{pid+1}.attribute.kst ~= 0
                            g{bid+1}.patch{pid+1}.attribute.kst = ...
                                k_new(k_old == g{bid+1}.patch{pid+1}.attribute.kst);
                        end
                        if g{bid+1}.patch{pid+1}.attribute.ken == nk_old
                            g{bid+1}.patch{pid+1}.attribute.ken = nk_new;
                        elseif g{bid+1}.patch{pid+1}.attribute.ken ~= 1
                            g{bid+1}.patch{pid+1}.attribute.ken = ...
                                k_new(k_old == g{bid+1}.patch{pid+1}.attribute.ken);
                        end
                    end
                    
                end
                
                % Update blade count
                for bid = bids_row{n}; g{bid+1}.bv.nblade = B{n}.N; g{bid+1}.bv.fblade = B{n}.N; end; 
                
            end
        end
        
        % Plot the new mesh
        figure(); hold on; grid on; box on; axis equal; j = round(M.n_j / 2);
        for bid = ts_get_bids(g)
            mesh(squeeze(g{bid+1}.x(:,j,:)),squeeze(g{bid+1}.rt(:,j,:)),squeeze(g{bid+1}.r(:,j,:)),...
                zeros(g{bid+1}.attribute.ni,g{bid+1}.attribute.nk,3))
        end
        
        % Calculate target mass flow and apply to boundary condition
        A = load([dr.ts 'inlet_cruise.mat']);
        r_mean = ((B{1}.xr_cas(1,2)^2 + B{1}.xr_hub(1,2)^2) / 2)^0.5;
        Vx_mean = phi(o) * r_mean * pi * (I.rpm(1) - I.rpm(2)) / 30;
        mass = A.ro * Vx_mean * B{1}.A_le * mass_fac; 
        for bid = ts_get_bids(g)
            for pid = ts_get_pids(g,bid)
                if g{bid+1}.patch{pid+1}.attribute.kind == 1
                    g{bid+1}.patch{pid+1}.throttle_type = 1;
                    g{bid+1}.patch{pid+1}.throttle_target = mass;
                    g{bid+1}.patch{pid+1}.throttle_k0 = 10;
                    g{bid+1}.patch{pid+1}.throttle_k1 = 10;
                    g{bid+1}.patch{pid+1}.throttle_k2 = 10;
                    g{bid+1}.patch{pid+1}.pout = 90e3;
                end
            end
        end
        
        % Apply initial guess
        for bid = ts_get_bids(g)
            for v = 1:length(guessnames)
%                 g{bid+1}.(guessnames{v})(:) = guess.(guessnames{v})(bid+1);
                g{bid+1}.(guessnames{v}) = repmat(guess.(guessnames{v}){bid+1},...
                    [1 1 g{bid+1}.attribute.nk]);
%                 g{bid+1}.(guessnames{v}) = guess.(guessnames{v}){bid+1};
            end
        end
        
        % Apply inlet boundary condition
        g = ts_apply_inlet(g,[dr.ts 'inlet_cruise_2percent.dat']);
        
        % Change the speed ratio
        for bid = ts_get_bids(g)
            rpm_ids = {'rpm' , 'rpmi1' , 'rpmi2' , 'rpmj1' , 'rpmj2' , 'rpmk1' , 'rpmk2'};
            for v = 1:length(rpm_ids)
                if g{bid+1}.bv.(rpm_ids{v}) < 0
                    g{bid+1}.bv.(rpm_ids{v}) = g{bid+1}.bv.(rpm_ids{v}) * (1 + rw);
                elseif g{bid+1}.bv.(rpm_ids{v}) > 0
                    g{bid+1}.bv.(rpm_ids{v}) = g{bid+1}.bv.(rpm_ids{v}) * (1 - rw);
                end
            end
        end
        
        % Update the run options
        g{1}.av.nstep = 100000; g{1}.av.nstep_save_start = 90000; g{1}.av.poisson_nstep = 10000;
%         g{1}.av.nstep = 1; g{1}.av.poisson_nstep = 2000;
%         g = ts_create_initial(g,100000,10000);

        % Submit the job
        Jobs{end+1,1} = ts_sbatch(job,g);
        
        else
            
            % Store job names for characteristics
            Jobs{end+1,1} = job;
        
        end
    end
    
    % Save all job names
    save([dr.ts 'CRF_jobs.mat'],'Jobs');
    
    % Remove empty jobs
    for n = length(Jobs):-1:1; if isfield(Jobs{n},'started') == 0; Jobs(n) = []; end; end;
    
    % Monitor all jobs
    ts_smonitor(Jobs);
    
end


%% Generate inlet boundary condition and guess files
if create_bc == 1
    
    % Flight condition
    T = 278.0957; P = 8.4438e+04; ro = 1.0566; V = 51.4; M = V / (ga * R * T)^0.5;
    flight_name = 'cruise';
    
    % Calculate stagnation quantities
    Comp_Tables; Po = P / interp1(Mach,P_Po,M); To = T / interp1(Mach,T_To,M); 
    save([dr.ts 'inlet_' flight_name '.mat'],'T','P','ro','V','M','Po','To')
    
    % Generate inlet boundary condition with 2% high boundary layers
    r = [0 0.02 0.98 1]; Po = [P Po Po P]; ni = 401; r_new = linspace(0,1,ni)';
    A = [r_new To * ones(ni,1) interp1(r,Po,r_new,'pchip') zeros(ni,3)];
    dlmwrite([dr.ts 'inlet_' flight_name '_2percent.dat'],A,'\t',1,0);
    
    
end


%% Run characteristics of all designs
if run_chics == 1
    
    % Load all jobs
%     load([dr.ts 'CRF_jobs.mat']);
    Jobs{1} = job; Jobs{1}.name = 'CRF2_P60S90R63DF35HT30_1';
    Jobs{2} = job; Jobs{2}.name = 'CRF2_P60S90R63DF35HT30_2';
    
    % Define cut planes for evaluation
    coords = [2 7 ; 7 23];
    
    % Run chics on all jobs
    ts_run_chic(Jobs,coords,1,[],0,0.02)
    
end


%% Plot output from fan designs
if plot_des == 1
    
    % Load test matrix
    load([dr.geom 'design_cube_2.mat']); 
    
    % Cruise velocity
    V_cruise = 51.4;
    
    % Process alldesigns in test matrix
    C = cell(size(phi)); B = cell(size(phi));
    for o = 1:numel(phi)
        
        % Job name
        des_name = ['CRF2_P' num2str(phi(o)*100) 'S' num2str(sig(o)*100) 'R' num2str(r_cas(o)*100) '_1'];
        job.name = [des_name '_1'];
    
        % Read the case and geometry
        g = ts_read_hdf5([dr.ts des_name '_avg.hdf5']);
        ts_export_paraview(g,[dr.pv des_name '.hdf5'],'HighSpeed');
%         ts_export_paraview(g,[dr.pv des_name '.hdf5'],'Track',[],1);
        B_temp = load([dr.geom des_name(1:end-2) '_splines.mat']); B{o} = B_temp.B;

        % Cut the rows
        C{o} = ts_cut_rows(g,[],[],0,0,[]);

        % Calculate fan 1D performance
        inlet = C{o}{1}.inlet; outlet = C{o}{2}.outlet; rpm = inlet.bv.rpm - outlet.bv.rpm;
        c = ts_calc_chic(inlet,outlet,rpm);

        % Calculate thrust and figure of merit
        [Vx,A] = ts_area_average(outlet,'Vx',3); [~,mdot] = ts_mass_average(outlet,'Vx',3);
        mdot = mdot * double(outlet.bv.nblade); A = A * double(outlet.bv.nblade);
        ro = ts_area_average(inlet,'ro',3);
        c.Wx = mdot * inlet.av.cp * (ts_mass_average(outlet,'To',3) - ts_mass_average(inlet,'To',3));
%         c.F = mdot * Vx; c.Mf = c.T^1.5 / (c.Wx * (2 * ro * A)^0.5); 
        c.F = mdot * (Vx / sig(o) - V_cruise);
        
        % Record data from design code
        load([dr.geom des_name(1:end-2) '.mat']); c.eta_p = d.d.eta_prop;
        
        % Calculate shaft speed and torque
        c.rpm = inlet.bv.rpm; c.w = pi * c.rpm / 30;
        c.T_1 = mdot * inlet.av.cp * (ts_mass_average(outlet,'To',3) - ...
            ts_mass_average(C{o}{1}.outlet,'To',3)) / c.w;
        c.T_2 = mdot * inlet.av.cp * (ts_mass_average(C{o}{1}.outlet,'To',3) - ...
            ts_mass_average(inlet,'To',3)) / c.w;
        
        % Flow mach numbers
        c.M_tip = max(max(ts_mass_average(C{o}{1}.inlet,'M_rel',2)),...
            max(ts_mass_average(C{o}{2}.inlet,'M_rel',2)));
        c.M_jet = max(ts_mass_average(C{o}{2}.outlet,'M',2)) / sig(o);
        
        % Blade interaction frequency
        c.bpf = double(inlet.bv.nblade) * abs(inlet.bv.rpm) / 60 + ...
            double(outlet.bv.nblade) * abs(outlet.bv.rpm) / 60;
        
        % Blade counts
        c.N_1 = double(inlet.bv.nblade); c.N_2 = double(outlet.bv.nblade); 
        
        % Record all data in arrays, preallocate on first case in loop
        scalnames = fieldnames(c);
        for v = 1:length(scalnames)
            if o == 1; p.(scalnames{v}) = zeros(size(phi)); end;
            p.(scalnames{v})(o) = c.(scalnames{v}); 
        end
        
    end
    
    % Increase resolution of cube by cubic interpolation
    p.phi = phi; p.sig = sig; p.D = r_cas * 2;
    varnames = fieldnames(p); 
    for v = 1:length(varnames)
        p.(varnames{v}) = interpn(p.(varnames{v}),3,'cubic'); 
    end
    
    % Save the efficiency data for visualisation
    pv_write_vtk(p,[dr.pv 'eta_cube_2.vtk'],{'phi','sig','D',},[],scalnames);
    save([dr.ts 'eta_cube_2.mat'],'p','B','C')
    
end


%% Plot deviation on all cases
if plot_dev == 1
    
    % Load data from file
    load([dr.ts 'eta_cube_2.mat'])
    
    % Plot whirl and blade angles at each location
    for o = 1:numel(C)
    
        % Create figure window
        figure(); r = linspace(0,1,100);
    
        % First rotor inlet plane
        subplot_tight(1,3,1,[0.1 0.05]); hold on; grid on; box on; title('Inlet');
        xlabel('Whirl Angle / ^\circ'); ylabel('Radius');
        plot(ts_mass_average(C{o}{1}.inlet,'Alpha',2),ts_mass_average(C{o}{1}.inlet,'r_nondim',2));
        plot(-fnval(B{o}{1}.chi_le,r),r);
        plot(ts_mass_average(C{o}{1}.inlet,'Alpha_rel',2),ts_mass_average(C{o}{1}.inlet,'r_nondim',2));
        legend('Absolute Whirl','B1 Metal Inlet','B1 Relative Whirl');

        % Second rotor inlet plane
        subplot_tight(1,3,2,[0.1 0.05]); hold on; grid on; box on; title('Betweeen Rows');
        xlabel('Whirl Angle / ^\circ'); set(gca,'yticklabel',[])
        plot(ts_mass_average(C{o}{1}.outlet,'Alpha_rel',2),ts_mass_average(C{o}{1}.outlet,'r_nondim',2));
        plot(-fnval(B{o}{1}.chi_te,r),r);
        plot(ts_mass_average(C{o}{2}.inlet,'Alpha_rel',2),ts_mass_average(C{o}{2}.inlet,'r_nondim',2));
        plot(fnval(B{o}{2}.chi_le,r),r);
        plot(ts_mass_average(C{o}{2}.inlet,'Alpha',2),ts_mass_average(C{o}{2}.inlet,'r_nondim',2));
        legend('B1 Relative Whirl','B1 Metal Outlet','B2 Relative Whirl',...
            'B2 Metal Inlet','Absolute Whirl');

        % Fan outlet
        subplot_tight(1,3,3,[0.1 0.05]); hold on; grid on; box on; title('Outlet');
        xlabel('Whirl Angle / ^\circ'); set(gca,'yticklabel',[])
        plot(ts_mass_average(C{o}{2}.outlet,'Alpha_rel',2),ts_mass_average(C{o}{2}.outlet,'r_nondim',2));
        plot(fnval(B{o}{2}.chi_te,r),r);
        plot(ts_mass_average(C{o}{2}.outlet,'Alpha',2),ts_mass_average(C{o}{2}.outlet,'r_nondim',2));
        legend('B2 Relative Whirl','B2 Metal Outlet','Absolute Whirl');
    
    end
    
    
end


%% Plot characteristics and calculate different fan operating points
if plot_chics == 1
    
    % Operating speeds and thrusts
    Vf = [0 48.9 55.4]; F = [6000 2700 1200]; alt = [0 0 5000] / 3.281; mar = {'x','o','s'};

    % Run number
    run_no = 3;
    
    % Define cut planes for evaluation
    coords = [2 7 ; 7 23];
    
    % Colour and marker lists
%     cols = repmat(reshape(repmat(reshape(lines(3),[1 3 3]),[3 1]),[],3),[3 1]); % Colour by nozzle area
%     cols = reshape(repmat(reshape(lines(3),[1 3 3]),[9 1]),[],3); % Colour by diameter
    %mar = reshape(repmat({'x','o','s'},[9 1]),[],1);
    
    % Plot middle diameter only
%     phi_des = phi_des(:,:,1); sig = sig(:,:,1); r_cas = r_cas(:,:,1); 
    
    % Plot specific cases
%     phi_des = phi_des(3,1,2); sig = sig(3,1,2); r_cas = r_cas(3,1,2); 
%     phi_des = phi_des(2,2,2); sig = sig(2,2,2); r_cas = r_cas(2,2,2); 
%     phi_des = phi_des(1,3,1); sig = sig(1,3,1); r_cas = r_cas(1,3,1); 
%     phi_des = phi_des(2,2,1); sig = sig(2,2,1); r_cas = r_cas(2,2,1);   
    
    % Plot size trend
%     phi_des = phi_des(2,2,:); sig = sig(2,2,:); r_cas = r_cas(2,2,:);   

    % Nine colours
    cols = [lines(7) ; 0.6 0.6 0.6 ; 0 0 0];

    % Colour by diameter for diffusion cases
%     cols = reshape(repmat(reshape(lines(4),[1 4 3]),[3 1 1]),[],3);

    % Plot different diffusion factors at current diameter
%     cols = lines(3); phi = phi(:,3); sig = sig(:,3); r_cas = r_cas(:,3); df = df(:,3); ht = ht(:,3); 
    
    % Figure windows
    h.chic = figure(); hold on; grid on; box on; 
    xlabel('Flow Coefficient'); ylabel('Total-Total Pressure Rise Coefficient');
    h.eta = figure(); hold on; grid on; box on; 
    xlabel('Flow Coefficient'); ylabel('Total-Total Fan Efficiency');
    
    % Process all designs in test matrix
    phi_des = phi; nc = zeros(size(phi_des)); leg = cell(numel(phi_des),1);
    for o = 1:numel(phi_des)
        
        % Process the chic
        if exist('df','var') == 0
            des_name = ['CRF2_P' num2str(phi_des(o)*100) 'S' num2str(sig(o)*100) 'R' num2str(r_cas(o)*100)];
        elseif exist('ht','var') == 0
            des_name = ['CRF2_P' num2str(phi_des(o)*100) 'S' num2str(sig(o)*100) ...
                'R' num2str(r_cas(o)*100) 'DF' num2str(df(o) * 100)];
        else
            des_name = ['CRF2_P' num2str(phi_des(o)*100) 'S' num2str(sig(o)*100) ...
                'R' num2str(r_cas(o)*100) 'DF' num2str(df(o) * 100) 'HT' num2str(ht(o) * 100)];
        end
        [inlets,outlets,c,h_temp] = ts_plot_chic(dr.ts,[des_name '_' num2str(run_no) '_M'],[],1,coords); 
        close(h_temp);
        
        % Legend
%         leg{o} = sprintf('Phi = %.1f Sigma = %.1f D = %.1f',[phi_des(o) sig(o) r_cas(o)*2]);
%         leg{o} = sprintf('DF = %.1f D = %.1f',[df(o) r_cas(o)*2]);
        leg{o} = sprintf('DF = %.1f D = %.1f',[df(o) r_cas(o)*2]);
        
        % Plot the chic
        figure(h.chic); plot(c.phi,c.Cptt,'.-','color',cols(o,:));
        figure(h.eta); plot(c.phi,c.eta_poly,'.-','color',cols(o,:));
        
        % Plot exit contours for selected designs
        if numel(phi_des) < 4
            figure; is = unique(round(linspace(1,length(outlets),4)));
            for i = 1:length(is)
                a = subplot_tight(1,4,i); hold on; axis equal; set(a,'visible','off');
                cut = outlets{is(i)}; v = linspace(0.8,0.95,20);
                contourf(cut.z,cut.y,cut.sfunc,v); caxis(v([1 end]));
            end
        end
        
        % Check the length of the chic
        nc(o) = length(inlets);
        if nc(o) > 1
            
            % Calculate tip Mach number for all points
            M_chic = zeros(size(inlets)); T_chic = zeros(size(inlets)); 
            omega_chic = pi * inlets{1}.bv.rpm / 30;
            for n = 1:length(inlets)
                M_chic(n) = max(ts_mass_average(inlets{n},'M_rel',2));
                T_chic(n) = ts_mass_average(inlets{n},'T',3);
            end
            M_chic = M_chic([1 1:end end]); T_chic = T_chic([1 1:end end]); 
            
            % Extend the chic to the axes
            phi_chic = c.phi; Cptt_chic = c.Cptt; eta_chic = [0 ; c.eta_poly ; 0];
            phi_chic = [interp1(Cptt_chic(1:2),phi_chic(1:2),0,'linear','extrap') ; phi_chic ; 0];
            Cptt_chic = [0 ; Cptt_chic ; 0];

            % Plot the extended chic
%             figure(h.chic); plot(phi_chic,Cptt_chic,'--','color',cols(o,:))
        
            % Fan disc area
            A = pi * (inlets{1}.r(end,1)^2 - inlets{1}.r(1,1)^2);

            % Calculate fan operating point at all vehicle operating points
            eta_f = zeros(size(Vf)); eta_p = zeros(size(Vf)); M_tip = zeros(size(Vf)); 
            phi = zeros(size(Vf)); Cptt = zeros(size(Vf)); omega = zeros(size(Vf)); T = zeros(size(Vf));
            for n = 1:length(Vf)

                % Air conditions at altitude
                [Ta,Pa,ro,~] = Altitude(alt(n));

                % Calculate axial velocity from thrust quadratic formula 
                Vx = (ro * A * Vf(n) + (ro^2 * A^2 * Vf(n)^2 + 4 * ro * A * F(n) / sig(o))^0.5) / ...
                    (2 * ro * A / sig(o));
                Vj = Vx / sig(o);
                
                % Bisect to find non-dimensional operating point
                phi_min = phi_chic(end); phi_max = phi_chic(1);
                for m = 1:15
                    
                    % Current flow coefficient and blade speed
                    phi(n) = 0.5 * (phi_min + phi_max);
                    U_1 = Vx / phi(n);
                    
                    % Pressure rise coefficient and 
                    Cptt(n) = interp1(phi_chic,Cptt_chic,phi(n),'linear','extrap');
                    U_2 = ((Vj^2 - Vf(n)^2) / Cptt(n))^0.5;
                    
                    % Update boundaries on flow coefficient
                    if U_1 > U_2
                        phi_min = phi(n);
                    else
                        phi_max = phi(n);
                    end

                end
                
                % Calculate fan and propulsive efficiencies
                eta_f(n) = interp1(phi_chic,eta_chic,phi(n),'linear','extrap');
                if Vf(n) > 0
                    eta_p(n) = 2 * Vf(n) / (Vf(n) + Vj);
                end
                
                % Compressor inlet and outlet conditions
                Po_1 = Pa + 0.5 * ro * Vf(n)^2; Po_2 = Pa + 0.5 * ro * Vj^2;
                To_1 = Ta; To_2 = To_1 * (Po_2 / Po_1)^((ga-1)/(ga * eta_f(n)));
                
                % Calculate power requirements, shaft speed and torque
                Wx = ro * A * Vx * cp * (To_2 - To_1);
                r = ts_area_average(inlets{1},'r',3); 
                omega(n) = U_1 / r; T(n) = Wx / (2 * omega(n));
                
                % Calculate figure of merit
                if Vf(n) == 0
                    eta_p(n) = F(n)^1.5 / (Wx * (2 * ro * A)^0.5); 
                end
                
                % Plot the operating point on the chic
                figure(h.chic); plot(phi(n),Cptt(n),'color',cols(o,:),'marker',mar{n})
                figure(h.eta); plot(phi(n),eta_f(n),'color',cols(o,:),'marker',mar{n})
                
                % Print motor requirements
                fprintf('Design %d\tPoint %d\tPower = %.1f kW\tShaft Speed = %.1f rad/s\n',...
                    [o n Wx / 2000 omega(n) / 2]); 
                
                % Scale Mach number to actual value at the correct conditions
                M_tip(n) = interp1(phi_chic,M_chic,phi(n),'linear','extrap') * omega(n) * ...
                    interp1(phi_chic,T_chic,phi(n),'linear','extrap')^0.5 / (omega_chic * Ta^0.5);
                
                % Print overall efficiency for cruise and climb cases
                fprintf('Design %d\tPoint %d\tFan Eta = %.1f\tProp. Eta = %.1f\tOv. Eta = %.1f\n',...
                    [o n eta_f(n)*100 eta_p(n)*100 eta_f(n)*eta_p(n)*100]);   
                
                % Print non-dimensionals for Dave's chart
%                 disp([F(n) / (A * ro * Vf(n)^2) F(n) * Vf(n) / Wx]);
                
            end
            
            % Save the chic
            save([dr.ts des_name '_' num2str(run_no) '_chic.mat'],'c','phi','Cptt','eta_f')
            
            % Calculate stall margin at static conditions
            sm = (phi(1) - phi_chic(end-1)) / (phi(end) - phi_chic(end-1));  
            
            % Record all variables
            clear c; c.eta_f = eta_f(end); c.eta_p = eta_p(end); c.eta_ov = c.eta_f * c.eta_p;
            c.eta_av = mean(eta_f); c.M_tip = M_tip(1); c.sm = sm; c.eta_cl = eta_f(2) * eta_p(2);
            c.omega_cr = omega(end); c.omega_st = omega(1); c.Mf = eta_p(1);
            c.T_cr = T(end); c.T_st = T(1);
            
            % Record all data in arrays, preallocate on first case in loop
            scalnames = fieldnames(c);
            for v = 1:length(scalnames)
                if o == 1; p.(scalnames{v}) = zeros(size(phi_des)); end;
                p.(scalnames{v})(o) = c.(scalnames{v}); 
            end
            
            % Print blade number
            disp([inlets{1}.bv.nblade outlets{1}.bv.nblade]);
            
        end        
    end    
    
    % Add the legend
    varnames = fieldnames(h); loc = {'northeast','southeast'};
    for v = 1:length(varnames)
        figure(h.(varnames{v}));
        l = flip(get(gca,'children'),1); 
        for a = length(l):-1:1; if strcmp(l(a).Marker,'.') == 0; l(a) = []; end; end;
        legend(l,leg,'location',loc{v});
    end
    
    % Save operating point data to file
    save([dr.ts 'chic_' test_name '_' num2str(run_no) '.mat'],'p')
    
    % Save data for PARAVIEW
    if numel(phi_des) == 27
        
        % Interpolate missing values
        varnames = fieldnames(p); 
        for v = 1:length(varnames)
            p.(varnames{v})(1,3,2) = 0.5 * (p.(varnames{v})(1,3,1) + p.(varnames{v})(1,3,3));
        end        
    
        % Increase resolution of cube by interpolation
        p.phi = phi_des; p.sig = sig; p.D = r_cas * 2;
        varnames = fieldnames(p); 
        for v = 1:length(varnames)
            p.(varnames{v}) = interpn(p.(varnames{v}),3,'linear'); 
            p.(varnames{v}) = smoothn(p.(varnames{v}),0.1);
        end

        % Save the efficiency data for visualisation
        pv_write_vtk(p,[dr.pv 'chic_cube_linear.vtk'],{'phi','sig','D',},[],scalnames);
        
    end
    
end


%% Create IBL files for drawing in Creo
if write_ibl == 1
    
    % Load design
    des_name = 'CRF2_P70S100R70';
    load([dr.geom des_name '_splines.mat'])
    
    % Create IBL files for all rows
    for n = 1:length(B)
        bl_write_ibl([dr.cad des_name '_' num2str(n) '.ibl'],B{n},n-2);
    end
    
    
end


%% Extract pressure loading data for specific fan cases
if pressure_load == 1
    
    % Load the case
    des_name = 'CRF2_P70S100R70_1_M9600';
    g = ts_read_hdf5([dr.ts des_name '_avg.hdf5']);
    
    % Take cuts and extract pressure and coordinates
    C = ts_cut_rows(g,[],[],0,0,[]);
    B = cell(size(C)); varnames = {'x','y','z','P'};
    for n = 1:length(B)
        c = ts_secondary(C{n}.blade);
        for v = 1:length(varnames)
            B{n}.(varnames{v}) = c.(varnames{v});
        end
    end
    
    % Save to file
    save([dr.cad des_name '_pressure.mat'],'B');
    
end


%% Plot efficiency trend with hub to tip ratio
if plot_hubrat == 1
    
    % Load data
    load([dr.geom 'design_hubrat.mat'])
    load([dr.ts 'chic_hubrat.mat']);
    
    % Open figure window
    figure(); hold on; grid on; box on; xlabel('Hub - Tip Ratio'); ylabel('Overall Efficiency');
    axis([ht(1) ht(end) 0.6 0.78]);
    
    % Correct lowest radius data
    p.eta_ov(1) = 0.71; p.eta_cl(1) = 0.61; 
    
    % Plot performance trends
    plot(ht,p.eta_ov,'.-')
    plot(ht,p.eta_cl,'.-')
    legend('Cruise','Climb','location','east');

    % Fit curve through data points
%     p = polyfit(ht(2:end),p.eta_ov(2:end),3);
%     x = linspace(ht(1),ht(end),100);
%     plot(x,polyval(p,x),'--')
%     plot(x,interp1(ht,p.eta_ov,x,'pchip'),'--')
    
    % Save the figure
    export_pres(gcf,[dr.fig 'hubtip_efficiency.pdf'],'full')

end

    
%% Plot efficiency trends with fan diameter and diffusion factor
if plot_diff == 1
    
    % Load data
    load([dr.geom 'design_diffusion.mat'])
    load([dr.ts 'chic_diffusion.mat']);
    
    % Increase resolution of test matrix by interpolation
    p.phi = phi; p.sig = sig; p.D = r_cas * 2; p.df = df;
    varnames = fieldnames(p); 
    for v = 1:length(varnames)
        p.(varnames{v}) = interpn(p.(varnames{v}),3,'linear'); 
        if any(strcmp(varnames{v},{'D','df'})) ~= 1
            p.(varnames{v}) = smoothn(p.(varnames{v}),0.2);
        end
    end
    
    % Colours
    cols = lines(7);
    
    % Create all figure windows
    varnames = {'stall','torque','cruise','climb','static'};
    for v = 1:length(varnames)
        h.(varnames{v}) = figure(); hold on; grid on; box on; 
        xlabel('Diffusion Factor'); ylabel('Diameter / m');
    end
    
    % Plot stall margin 
    figure(h.stall); title('Stall Margin'); v = linspace(-0.6,0.3,10);
    contourf(p.df,p.D,p.sm,v); colorbar; caxis(v([1 end]));
    contour(p.df,p.D,p.sm,[0 0],'linewidth',2,'linecolor',cols(2,:))

    % Specify existing motor designs
    Ct_mot = [0.058 0.192]; v = diff(Ct_mot) * linspace(-1,2,11) + Ct_mot(1);  
    
    % Motor torque coefficient
    figure(h.torque); title('Torque Coefficient');
    Ct = p.T_st ./ p.omega_st.^2;
    contourf(p.df,p.D,Ct,v); colorbar; caxis([min(Ct(:)) max(Ct(:))]);
    contour(p.df,p.D,Ct,Ct_mot,'linewidth',2,'linecolor',[0.8 0.8 0.8])
    contour(p.df,p.D,p.sm,[0 0],'linewidth',2,'linecolor',cols(2,:))
    
    % Efficiency across the different operating points
    varnames = {'eta_ov','eta_cl','Mf'}; fignames = {'cruise','climb','static'};
    titlenames = {'Cruise Overall Efficiency','Climb Overall Efficiency','Static Figure of Merit'};
    for v = 1:length(varnames)
        figure(h.(fignames{v})); title(titlenames{v}); colorbar;
        contourf(p.df,p.D,p.(varnames{v})); caxis([min(p.(varnames{v})(:)) max(p.(varnames{v})(:))]); 
        contour(p.df,p.D,Ct,Ct_mot,'linewidth',2,'linecolor',[0.8 0.8 0.8])
        contour(p.df,p.D,p.sm,[0 0],'linewidth',2,'linecolor',cols(2,:))
    end

    % Save all figures and add boxes
    varnames = fieldnames(h);
    for v = 1:length(varnames)
        export_pres(h.(varnames{v}),[dr.fig 'temp.pdf'],'full'); add_box(h.(varnames{v}));
        export_pres(h.(varnames{v}),[dr.fig 'diffusion_' varnames{v} '.pdf'],'full');
        export_pres(h.(varnames{v}),[dr.fig 'diffusion_' varnames{v} '.png'],'full');
    end
    
end


%% Plot characteristics and performance of design run with different speed ratios
if plot_speed == 1
    
    % Figure windows
    h.chic = figure(); hold on; grid on; box on; axis([0.4 0.7 0.12 0.32]);
    xlabel('Flow Coefficient'); ylabel('Total-Total Pressure Rise Coefficient');
    h.eta = figure(); hold on; grid on; box on; axis([0.4 0.7 0.6 0.9]);
    xlabel('Flow Coefficient'); ylabel('Total-Total Fan Efficiency');
    
    % Colours, markers and legend
    cols = lines(3); leg = {'1.00 : 1.00','1.15 : 0.85','0.85 : 1.15'}; mar = {'x','o','s'};
    
    % Plot the performance chics
    for o = 1:3
        
        % Load chic data
        load([dr.ts 'CRF2_P60S90R63DF35HT30' '_' num2str(o) '_chic.mat'])
        
        % Plot the chic
        figure(h.chic); l.chic(o) = plot(c.phi,c.Cptt,'.-','color',cols(o,:));
        figure(h.eta); l.eta(o) = plot(c.phi,c.eta_poly,'.-','color',cols(o,:));
        
        % Plot the operating points on the chic
        for n = 1:3     
            figure(h.chic); ob = plot(phi(n),Cptt(n),'color',cols(o,:),'marker',mar{n},...
                'markerfacecolor',cols(o,:));
            figure(h.eta); ob = plot(phi(n),eta_f(n),'color',cols(o,:),'marker',mar{n},...
                'markerfacecolor',cols(o,:));
            if n == 1; set(ob,'markersize',10); end;
        end
        
    end
    
    % Add the legends
    varnames = fieldnames(h); loc = {'northeast','southeast'};
    for v = 1:length(varnames)
        legend(l.(varnames{v}),leg,'location',loc{v});
    end
    
    % Save the figures
    export_pres(h.chic,[dr.fig 'speed_cptt.pdf'],'half');
    export_pres(h.eta,[dr.fig 'speed_eta.pdf'],'half');
    
end




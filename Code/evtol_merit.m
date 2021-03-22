function fig_merit = evtol_merit(x)
    close all
    eVTOL
    global run_case ulim ub_actual lb_actual ub_scaled lb_scaled
    
    
    n_base = 10; n_num = 1; % what base and how many numbers/ letters per parameter change this to change how big the sweeps are
    
    normalisedb = (x - lb_scaled) ./ (ub_scaled - lb_scaled); % normalise 0<x<1
    x_scale = (n_base^n_num - 1) * normalisedb'; % normalise 0<=x<n_base(00<x<ZZ, if n_base =36 and n_num =2)
    iter_id = reshape(dec2base(round(x_scale),n_base,n_num)',1,[]); % (...,1,[]) gives 1 x whatever

    chi_le = x(1)/ulim * (ub_actual(1) - lb_actual(1)) + lb_actual(1);
    chi_te = x(2)/ulim * (ub_actual(2) - lb_actual(2)) + lb_actual(2);
    sweep = x(3)/ulim * (ub_actual(3) - lb_actual(3)) + lb_actual(3);
    N = x(4);
    lean = x(5)/ulim * (ub_actual(5) - lb_actual(5)) + lb_actual(5);
    
    
    % Current job variables
    job.priority = 'MILLER-SL2-GPU';
    job.cards = 1;
    job.walltime = 0.50;

    % Options for using interactive node and command to activate on Wilkes
    % sintr -p pascal -t 8:00:00 -A MILLER-SL2-GPU --gres=gpu:1 -N 1
    job.node = 'gpu-e-52';
    
    % Mesh parameters
    dx_duct = 0.2; tol = 2e-4;
    ni_in = 17; ni_out = 25; ni_bl = 81; nj = 57; nk = 41;
    e_b = 5e-5; e_j = 5e-5; e_le = 6e-4; e_te = 8e-4;

    % Boundary conditions
    Po_in = 101325; To_in = 288.15; Alpha_in = 23; Vx_in = 80; d_hub = 0.1; d_cas = 0.1;

    % Number of iterations to run
    nstep = 20000; nstep_poisson = 2000; nstep_avg = 5000;
    
    %% Set job name
    job.name = [run_case '_' num2str(iter_id)];
    
    %Check if already run
    if exist([dr.geom job.name '.mat'])  ==2   
        disp(['Found ' job.name ' skipping TURBOSTREAM'])
        load([dr.fig run_case],'designstruct'); % load table
        load([dr.geom job.name '.mat'],'designvars'); % load the run

        designstruct.name(end+1,1) = designvars.name;
        designstruct.chi_le(end+1,1) = designvars.chi_le;
        designstruct.chi_te(end+1,1) = designvars.chi_te;
        designstruct.sweep(end+1,1) = designvars.sweep;
        designstruct.N(end+1,1) = designvars.N;
        designstruct.lean(end+1,1) = designvars.lean;
        designstruct.merit(end+1,1) = designvars.merit;
        
        fig_merit = designvars.merit;
        
        save([dr.fig run_case],'designstruct');
        save([dr.geom job.name '.mat'],'designvars'); 
    else
    %% Create the geometry

    % Load stator example
    load([dr.geom 'NRF_3b_splines.mat']); b = B{2};

    % Increase 3D stacking
    b.sweep.coefs = [0 sweep sweep sweep 0];
    b.lean.coefs = [0 lean lean lean 0];
    % Modify hub and casing lines
    b.x_in = 0.1; b.x_out = 0.35; b.A_te = 1.1;
    b.drdx_hub_le = -0.15; b.drdx_hub_te = -0.15;
    b.drdx_cas_le = 0.15; b.drdx_cas_te = 0.15;

    % Apply blade angles and count
    b.chi_le = spap2(1,2,[0 1],[1 1] * chi_le);
    b.chi_te = spap2(1,2,[0 1],[1 1] * chi_te);
    b.N = N; dt = 2 * pi / N;
    % Create and plot the blade
    b = bl_construct_blade(b,1);

    % Create the meridional lines and extend with constant area ducts
    m = bl_construct_mer({b},tol); xr_hub = m.xr_hub; xr_cas = m.xr_cas;
    x_in = linspace(b.x_in-dx_duct,b.x_in,dx_duct/tol)'; 
    x_out = linspace(b.x_out,b.x_out+dx_duct,dx_duct/tol)';
    xr_hub = [x_in xr_hub(1,2) * ones(size(x_in)) ; xr_hub(2:end-1,:) ; ...
        x_out xr_hub(end,2) * ones(size(x_out))];
    xr_cas = [x_in xr_cas(1,2) * ones(size(x_in)) ; xr_cas(2:end-1,:) ; ...
        x_out xr_cas(end,2) * ones(size(x_out))];

    % Plot the meridional lines
    plot3(xr_hub(:,1),zeros(size(xr_hub,1),1),xr_hub(:,2),'k-')
    plot3(xr_cas(:,1),zeros(size(xr_cas,1),1),xr_cas(:,2),'k-')
    
    %% Generate the mesh

    % Correct number of streamwise points to include inlet and outlet
    ni = ni_bl + ni_in + ni_out - 2;

    % Split the blade based upon minimum and maximum axial coordinates and reinterpolate along chord
    xrrti = cell(2,1); 
    for j = 1:size(b.xrrt,2)

        % Find leading and trailing edge indices
        [~,i_le] = min(b.xrrt(:,j,1)); [~,i_te] = max(b.xrrt(:,j,1)); 

        % Split into two sides
        [xrrt{1},xrrt{2}] = bl_split(squeeze(b.xrrt(:,j,:)),i_le,i_te);

        % Offset second side by a pitch
        xrrt{2}(:,3) = xrrt{2}(:,3) + xrrt{2}(:,2) * dt;

        % Interpolate to correct spacing
        for n = 1:2

            % Extend the surface to inlet and outlet planes
            xrrt{n} = xrrt{n}([1 1:end end],:);
            xrrt{n}([1 end],1) = xr_hub([1 end],1);

            % Pinch the leading and trailing edges to prevent negative volumes
            xrrt{n}(2,1) = xrrt{n}(2,1) - 0.5 * e_le;
            xrrt{n}(end-1,1) = xrrt{n}(end-1,1) + 0.5 * e_te;

            % Generate new spacings
            s = dist_2d(xrrt{n}(:,[1 2]),0);
            si = 0.5 * (geo_bunch(s(2),s(end-1),e_le,e_te,ni_bl) + ...
                poly_bunch(s(2),s(end-1),e_le,e_te,ni_bl));

            % Add on inlet and outlet spacings
            si_in = flip(s(2) - geo_bunch(s(1),s(2),e_le,[],ni_in),2);
            si_out = geo_bunch(s(end-1),s(end),e_te,[],ni_out);
            si = [si_in si(2:end-1) si_out];

            % Interpolate new coordinates
            i = [1:ni_in (ni-ni_out+1):ni];
            xrrti{n}(i,j,:) = reshape(interp1(s,xrrt{n},si(i),'linear','extrap'),[length(i) 1 3]);
            i = ni_in:(ni-ni_out+1);
            xrrti{n}(i,j,:) = reshape(interp1(s,xrrt{n},si(i),'pchip'),[ni_bl 1 3]);      

        end
    end

    % Interpolate to correct radial spacing
    xrrtj = cell(2,1); 
    for i = 1:ni
        for n = 1:2

            % Calculate non-dimensional radial coordinates
            xrrt = squeeze(xrrti{n}(i,:,:));
            r_hub = interp1(xr_hub(:,1),xr_hub(:,2),xrrt(:,1),'pchip');
            r_cas = interp1(xr_cas(:,1),xr_cas(:,2),xrrt(:,1),'pchip');
            r_nondim = (xrrt(:,2) - r_hub) ./ (r_cas - r_hub);

            % Calculate new spanwise spacing
            s = dist_2d(xrrt(:,1:2),0); s_mer = interp1(r_nondim,s,[0 1],'pchip');
            sj = 0.8 * geo_bunch(s_mer(1),s_mer(2),e_j,e_j,nj) + ...
                0.2 * poly_bunch(s_mer(1),s_mer(2),e_j,e_j,nj);

            % Interpolate to correct spacing
            xrrtj{n}(i,:,:) = reshape(interp1(s,xrrt,sj,'spline'),[1 nj 3]);        
        end
    end

    % Plot the blade surfaces
    figure(); hold on; grid on; box on; axis equal; cols = lines(7);
    i = ni_in:(ni_in+ni_bl-1);
    mesh(xrrtj{1}(i,:,1),-xrrtj{1}(i,:,3),xrrtj{1}(i,:,2),zeros(ni_bl,nj,3))
    mesh(xrrtj{2}(i,:,1),-xrrtj{2}(i,:,3),xrrtj{2}(i,:,2),zeros(ni_bl,nj,3))

    % Generate passage block by applying pitchwise spacing
    xrrt_pass = zeros(ni,nj,nk,3);
    for i = 1:ni
        for j = 1:nj

            % Pitchwise spapcing at current location
            s = xrrtj{2}(i,j,3) - xrrtj{1}(i,j,3);

            % Interpolate spacing factor to blend between geometric and linear
            f = interp1(xrrtj{1}([1 ni_in ni_in+ni_bl-1 end],j,1),[0 1 1 0],xrrtj{1}(i,j,1),'linear');

            % New spacing
            sk = f * (0.7 * geo_bunch(0,s,e_b,e_b,nk) + 0.3 * poly_bunch(0,s,e_b,e_b,nk)) ...
                + (1-f) * linspace(0,s,nk); 
            sk = repmat(reshape(sk / sk(end),[1 1 nk 1]),[1 1 1 3]);

            % Block coordinates
            xrrt_pass(i,j,:,:) = repmat(reshape(xrrtj{1}(i,j,:),[1 1 1 3]),[1 1 nk 1]) .* (1 - sk) + ...
                repmat(reshape(xrrtj{2}(i,j,:),[1 1 1 3]),[1 1 nk 1]) .* sk;

        end
    end

    % Plot the endwall surfaces
    mesh(squeeze(xrrt_pass(:,1,:,1)),-squeeze(xrrt_pass(:,1,:,3)),squeeze(xrrt_pass(:,1,:,2)),...
        repmat(reshape(cols(1,:),[1 1 3]),[ni nk 1]));
    mesh(squeeze(xrrt_pass(:,end,:,1)),-squeeze(xrrt_pass(:,end,:,3)),squeeze(xrrt_pass(:,end,:,2)),...
        repmat(reshape(cols(1,:),[1 1 3]),[ni nk 1]));

    % Plot periodic patches
    i = cell(2,1); i{1} = 1:ni_in; i{2} = (ni_in+ni_bl-1):ni;
    for k = [1 nk]
        for n = 1:length(i)
            mesh(squeeze(xrrt_pass(i{n},:,k,1)),-squeeze(xrrt_pass(i{n},:,k,3)),...
                squeeze(xrrt_pass(i{n},:,k,2)),repmat(reshape(cols(2,:),[1 1 3]),[length(i{n}) nj 1]),...
                'facecolor','none');
        end
    end

    % Check for negative axial spacings
    q = diff(xrrt_pass(:,:,:,1),1,1) < 0;
    if any(q); error([num2str(sum(q(:) < 0)) ' Negative axial spacings']); end;


    %% Patch together and preprocess for TURBOSTREAM

    % Add mesh to the first block in cell array
    g = cell(1); varnames = {'x' 'r' 'rt'}; 
    for v = 1:length(varnames); g{1}.(varnames{v}) = xrrt_pass(:,:,:,v); end;

    % Block attributes
    a.bid = 0; a.np = 2; a.ni = ni; a.nj = nj; a.nk = nk; a.nc = 0; a.ncl = 0; a.procid = 0; a.threadid = 0;
    g{1}.attribute = a; g{1}.patch = cell(2,1); g{1}.bv.nblade = N; g{1}.bv.fblade = N;

    % Approximate dynamic head
    ro = (Po_in / (R * To_in));
    dP = 0.5 * ro * (Vx_in / cosd(Alpha_in))^2;

    % Calculate mass flow 
    mass = pi * (xr_cas(1,2)^2 - xr_hub(1,2)^2) * ro * Vx_in;

    % Default patch attributes
    p.idir = 0; p.jdir = 1; p.kdir = 2; p.nface = 0; p.nt = 1; p.nxbid = 0; p.nxpid = 0;

    % Apply exit boundary condition
    p.ist = ni-1; p.ien = ni; p.jst = 0; p.jen = nj; p.kst = 0; p.ken = nk;
    p.bid = 0; p.pid = 0; p.kind = 1;
    g{1}.patch{1}.attribute = p;
    g{1}.patch{1}.pout = Pref - dP; g{1}.patch{1}.throttle_type = 1; g{1}.patch{1}.throttle_target = mass;
    g{1}.patch{1}.throttle_k0 = 50; g{1}.patch{1}.throttle_k1 = 50; g{1}.patch{1}.throttle_k2 = 50;
    g{1}.patch{1}.ipout = 3;

    % Turn off mass flow control
    g{1}.patch{1}.throttle_type = 0;

    % Inlet stagnation pressure profile
    inlet.r_nondim = (g{1}.r(1,:,1) - g{1}.r(1,1,1)) / diff(g{1}.r(1,[1 end],1));
    inlet.Po = interp1([0 d_hub 1-d_cas 1],[-dP 0 0 -dP] + Po_in,inlet.r_nondim,'pchip');

    % Remaining inlet boundary condition application
    inlet.Alpha = ones(1,nj) * Alpha_in; inlet.Beta = zeros(1,nj); inlet.To = ones(1,nj) * To_in;
    p.ist = 0; p.ien = 1; p.pid = 1; p.kind = 0; g{1}.patch{2}.attribute = p;
    g = ts_apply_inlet(g,inlet,1);

    % Periodic patch across pitch upstream of the blade
    g = ts_add_periodic_patch(0,0,[0 ni_in ; 0 nj ; 0 1],[0 ni_in ; 0 nj ; nk-1 nk],[0 1 2],[0 1 2],g,5);

    % Periodic patch across pitch downstream of the blade
    g = ts_add_periodic_patch(0,0,[ni-ni_out ni ; 0 nj ; 0 1],[ni-ni_out ni ; 0 nj ; nk-1 nk],...
        [0 1 2],[0 1 2],g,5);

    % Create default TURBOSTREAM settings
    g = ts_create_initial(g,nstep,nstep_poisson,[ro Vx_in 0 0 2e5 0.0002]);

    % Reduce smoothing but run with mixng length turbulence model
    g = ts_create_restart(g,nstep,1); g{1}.av.nstep_save_start = nstep - nstep_avg;
    g{1}.av.poisson_nstep = nstep_poisson; g{1}.av.poisson_restart = 0;

    % Check patches are valid
    %ts_check_patches(g,2);
    
    %% adding values to structure
    
    load([dr.fig run_case],'designstruct');
    
    designstruct.name(end+1,1) = {num2str(iter_id)};
    designstruct.chi_le(end+1,1) = chi_le;
    designstruct.chi_te(end+1,1) = chi_te;
    designstruct.sweep(end+1,1) = sweep;
    designstruct.N(end+1,1) = N;
    designstruct.lean(end+1,1) = lean;

    %% Submit and monitor the job

    % Set the job name
    %job.name = [run_case '_' num2str(iter_id)];

    % Submit the job to Wilkes
    job = ts_sbatch(job,g);
    
    ts_smonitor(job);

    % Plot the job convergence history
    ts_plot_conv(job);

    % Export for viewing in PARAVIEW
    g = ts_read_hdf5([dr.ts job.name '_avg.hdf5']);
    ts_export_paraview(g,[dr.pv job.name '.hdf5'],'Track',[],1);


    %% Post-process the job
    
    % Record for each run the paramter and then add them into a table
    
    
    designvars.name = {num2str(iter_id)};
    designvars.chi_le = chi_le;
    designvars.chi_te = chi_te;
    designvars.sweep = sweep;
    designvars.N = N;
    designvars.lean = lean;
   

     %Take cuts at inlet and outlet
    inlet = ts_structunstruct_cut(g,[b.x_in 0 ; b.x_in 10],[],2,1); % inlet plane
    outlet = ts_structunstruct_cut(g,[b.x_out 0 ; b.x_out 10],[],2,1); % outlet plane

    
    inlet_area = mass/ ts_mass_average(inlet,'rovx',3); % effective inlet area
    outlet_area = mass/ ts_mass_average(outlet,'rovx',3); % effective outlet area
   
    fig_merit = sqrt(2*outlet_area/inlet_area); % figure of merit 
    
    designstruct.merit(end+1,1) = fig_merit;
    save([dr.fig run_case],'designstruct');
    
    designvars.merit = fig_merit;
    save([dr.geom job.name],'designvars')
    disp(['figure of merit: ' num2str(fig_merit)])
   

end
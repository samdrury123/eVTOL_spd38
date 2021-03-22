function g = ts_create_initial(g,nstep,poisson_nstep,initial_value)
% TS_CREATE_INITIAL  Apply default options to a TURBOSTREAM mesh
%
%   g = TS_CREATE_INITIAL(g,nstep,poisson_nstep,initial_value)
%
%   g - cell array of TURBOSTREAM solution
%   nstep - number of main loop iterations to use
%   poisson_nstep - number of poisson loop iterations to use
%   initial_value - vector of primary flow variables or string pointing to a text file for initial guess

% Set poisson to restart if no steps are requested
if poisson_nstep == 0
    poisson_restart = 1;
else
    poisson_restart = 0;
end

% Create default application variables
av.cp = 1005.0;
av.ga = 1.4;
av.cfl = 0.4;
av.prandtl = 1.0;
av.viscosity = 1.79E-5;
av.sfin = 0.01;
av.dampin = 10.0;
av.rfvis = 0.2;
av.poisson_sfin = 0.05;
av.poisson_cfl = 0.7;
av.restart = 1;
av.ilos = 1;
av.nlos = 5;
av.dts = 0;
av.ncycle = 0;
av.frequency = 0.0;
av.nstep_cycle = 0;
av.dts_conv = 0.0;
av.facsafe = 0.0;
av.nstep_inner = 0;
av.nstep_save = 0;
av.nstep_save_start = 100000;
av.nstep_save_probe = 0;
av.nstep_save_start_probe = 0;
av.rfmix = 0.2;
av.facsecin = 0.8;
av.viscosity_law = 1;
av.wall_law = 0;
av.write_yplus = 0;
av.write_force = 0;
av.adaptive_smoothing = 1;
av.ifgas = 0;
av.tref = 300.0;
av.pref = 1e5;
av.rg_cp0 = 1005.0;
av.rg_cp1 = 0.0;
av.rg_cp2 = 0.0;
av.rg_cp3 = 0.0;
av.rg_cp4 = 0.0;
av.rg_cp5 = 0.0;
av.rg_rgas = 287.15;
av.cp0_0 = 1005.0;
av.cp1_0 = 0.0;
av.cp2_0 = 0.0;
av.cp3_0 = 0.0;
av.cp4_0 = 0.0;
av.cp5_0 = 0.0;
av.rgas_0 = 287.15;
av.cp0_1 = 1005.0;
av.cp1_1 = 0.0;
av.cp2_1 = 0.0;
av.cp3_1 = 0.0;
av.cp4_1 = 0.0;
av.cp5_1 = 0.0;
av.rgas_1 = 287.15;
av.nspecies = 1;
av.sf_scalar = 0.05;
av.schmidt_0 = 1.0;
av.schmidt_1 = 1.0;
av.fac_stmix = 0.0;
av.fac_st0 = 1.0;
av.fac_st1 = 1.0;
av.fac_st2 = 1.0;
av.fac_st3 = 1.0;
av.fac_sa_step = 1.0;
av.fac_sa_smth = 4.0;
av.turb_vis_damp = 1.0;
av.fac_wall = 1.0;
av.if_no_mg = 0;
av.turbvis_lim = 3000.0;
av.nomatch_int = 1;
av.poisson_nsmooth = 10;
av.write_tdamp = 0;
av.write_egen = 0;
av.viscosity_a1 = 0.0;
av.viscosity_a2 = 0.0;
av.viscosity_a3 = 0.0;
av.viscosity_a4 = 0.0;
av.viscosity_a5 = 0.0;
av.ifsuperfac = 0;
av.poisson_limit = 0;
av.cfl_ko = 0.4;
av.cfl_st_ko = 0.01;
av.cfl_en_ko = 0.4;
av.sfin_sa = 0.05;
av.sfin_ko = 0.05;
av.ko_dist = 1e-4;
av.ko_restart = 0;
av.fac_st0_option = 0;
av.sa_helicity_option = 0;
av.sa_ch1 = 1.0;
av.sa_ch2 = 1.0;
      
% Apply default application variables
varnames = fieldnames(av);
for v = 1:length(varnames)
    g{1}.av.(varnames{v}) = av.(varnames{v});
end

% Create default block variables
bv.fmgrid = 0.4;
bv.poisson_fmgrid = 0.1;
bv.xllim = 0.0;
bv.vgridin = 100.0;
bv.pstatin = 1e5;
bv.tstagin = 300.0;
bv.vgridout = 100.0;
bv.pstatout = 1e5;
bv.tstagout = 300.0;
bv.ftype = 0;
bv.fracann = 1.0;
bv.sfin_mul = 1.0;
bv.facsecin_mul = 1.0;
bv.dampin_mul = 1.0;
bv.fsturb = 1.0;
bv.nimixl = 0;
bv.srough_i0 = 0.0;
bv.srough_i1 = 0.0;
bv.srough_j0 = 0.0;
bv.srough_j1 = 0.0;
bv.srough_k0 = 0.0;
bv.srough_k1 = 0.0;

bv.rpm = 0.0;
bv.rpmi1 = 0.0;
bv.rpmi2 = 0.0;
bv.rpmj1 = 0.0;
bv.rpmj2 = 0.0;
bv.rpmk1 = 0.0;
bv.rpmk2 = 0.0;

bv.itrans = 0;
bv.itrans_j1_st = 0;
bv.itrans_j2_st = 0;
bv.itrans_k1_st = 0;
bv.itrans_k2_st = 0;
bv.itrans_j1_en = 0;
bv.itrans_j2_en = 0;
bv.itrans_k1_en = 0;
bv.itrans_k2_en = 0;
bv.itrans_j1_frac = 0.0;
bv.itrans_j2_frac = 0.0;
bv.itrans_k1_frac = 0.0;
bv.itrans_k2_frac = 0.0;

bv.jtrans = 0;
bv.jtrans_i1_st = 0;
bv.jtrans_i2_st = 0;
bv.jtrans_k1_st = 0;
bv.jtrans_k2_st = 0;
bv.jtrans_i1_en = 0;
bv.jtrans_i2_en = 0;
bv.jtrans_k1_en = 0;
bv.jtrans_k2_en = 0;
bv.jtrans_i1_frac = 0.0;
bv.jtrans_i2_frac = 0.0;
bv.jtrans_k1_frac = 0.0;
bv.jtrans_k2_frac = 0.0;

bv.ktrans = 0;
bv.ktrans_i1_st = 0;
bv.ktrans_i2_st = 0;
bv.ktrans_j1_st = 0;
bv.ktrans_j2_st = 0;
bv.ktrans_i1_en = 0;
bv.ktrans_i2_en = 0;
bv.ktrans_j1_en = 0;
bv.ktrans_j2_en = 0;
bv.ktrans_i1_frac = 0.0;
bv.ktrans_i2_frac = 0.0;
bv.ktrans_j1_frac = 0.0;
bv.ktrans_j2_frac = 0.0;

bv.dro_avg = 0;
bv.dro_superfac = 0;
bv.droe_avg = 0;
bv.droe_superfac = 0;
bv.drorvt_avg = 0;
bv.drorvt_superfac = 0;
bv.drovr_avg = 0;
bv.drovr_superfac = 0;
bv.drovx_avg = 0;
bv.drovx_superfac = 0;
bv.dtrans_dyn_vis = 0;

bv.superstep = 0.0;
bv.superfac = 0.0;
bv.fac_st0 = 1.5;

bv.ndup_phaselag = 1;
bv.turb_intensity = 5.0;

% Apply default block variables
varnames = fieldnames(bv);
for bid = ts_get_bids(g)
    for v = 1:length(varnames)
        g{bid+1}.bv.(varnames{v}) = bv.(varnames{v});
    end
end

% Tune variables appropriate to an initial run

% Slower and more stable poisson loop
g{1}.av.poisson_nstep = poisson_nstep;
g{1}.av.poisson_restart = poisson_restart;
g{1}.av.poisson_cfl = 0.5;

% More stable main loop
g{1}.av.nstep = nstep;
g{1}.av.nstep_save_start = nstep+1;
g{1}.av.nchange = 10000;
g{1}.av.facsecin = 0.02;
g{1}.av.sfin = 2.5;
g{1}.av.cfl = 0.3;
g{1}.av.dampin = 3;
g{1}.av.rfmix = 0;
g{1}.av.ilos = 1;
g{1}.av.restart = 1;

% Set block variables
for bid = ts_get_bids(g)
    % Multigrid
    g{bid+1}.bv.poisson_fmgrid = 0.08;
    g{bid+1}.bv.fmgrid = 0.25;

    % Freestream mixing
    g{bid+1}.bv.free_turb = 0.05;
    g{bid+1}.bv.fsturb = 2;
    
    % Mixing length
    g{bid+1}.bv.xllim_free = 0.1;
    r_nondim = (g{bid+1}.r(1,:,1) - min(g{bid+1}.r(1,:,1))) / ...
        (max(g{bid+1}.r(1,:,1)) - min(g{bid+1}.r(1,:,1)));
    [~, J] = min(abs(r_nondim-0.5));
    nblade = double(g{bid+1}.bv.nblade);
    pitch = g{bid+1}.r(1,J,1) * 2 * pi / nblade;
    g{bid+1}.bv.xllim = 0.03 * pitch;  
    
    % No multigrid at mixing planes
    for pid = ts_get_pids(bid,g)
        if g{bid+1}.patch{pid+1}.attribute.kind == 2
            g{bid+1}.bv.fmgrid = 0;
        end
    end
end

% Create initial guess if necessary
if exist('initial_value','var') ~= 0
    
    % Choose which initial guess to use
    if ischar(initial_value) == 1 && isempty(strfind(initial_value,'.hdf5')) == 0
        
        % Copy an initial guess from a converged solution point for point
        g_guess = ts_read_hdf5(initial_value);
        varnames = {'ro' 'rovx' 'rovr' 'rorvt' 'roe' 'trans_dyn_vis' 'phi'};
        for bid = ts_get_bids(g)
            for v = 1:length(varnames)
                g{bid+1}.(varnames{v}) = g_guess{bid+1}.(varnames{v});
            end
        end
        
    elseif ischar(initial_value) == 1
        
        % Read initial guess through the machine from file
        guess = dlmread(initial_value);
        
        % Determine type of guess to use, derived or primary flow variables
        if size(guess,2) == 7
        
            % Extend initial guess
            guess = [guess(1,:) ; guess ; guess(end,:)];
            guess(1,2) = -1e6; guess(end,2) = 1e6;

            % Extract guess parameters from file
            x_guess = guess(:,2); Vx_guess = guess(:,3); Vr_guess = guess(:,4);
            Vt_guess = guess(:,5); T_guess = guess(:,6); P_guess = guess(:,7);
            
            % Initialise other values
            trans_dyn_vis = 0.0002;
            
            % Loop over every block
            for bid = ts_get_bids(g)
                
                % Interpolate derived properties
                Vx = interp1(x_guess,Vx_guess,g{bid+1}.x);
                Vr = interp1(x_guess,Vr_guess,g{bid+1}.x);
                Vt = interp1(x_guess,Vt_guess,g{bid+1}.x);
                V = (Vx.^2 + Vr.^2 + Vt.^2).^0.5;
                T = interp1(x_guess,T_guess,g{bid+1}.x);
                P = interp1(x_guess,P_guess,g{bid+1}.x);

                % Calculate primary flow variables
                g{bid+1}.ro = P ./ (287 * T);
                g{bid+1}.rovx = g{bid+1}.ro .* Vx;
                g{bid+1}.rovr = g{bid+1}.ro .* Vr;
                g{bid+1}.rorvt = g{bid+1}.ro .* g{bid+1}.r .* Vt;
                g{bid+1}.roe = g{bid+1}.ro .* (T * 718 + 0.5 * V.^2);

                g{bid+1}.trans_dyn_vis = ones(size(g{bid+1}.ro)) * trans_dyn_vis;
                
            end
            
        elseif size(guess,2) == 8
            
            % Apply primary flow variables directly
            varnames = {'ro' 'rovx' 'rovr' 'rorvt' 'roe' 'trans_dyn_vis' 'phi'};
            for bid = ts_get_bids(g)
                for v = 1:length(varnames)
                    g{bid+1}.(varnames{v}) = interp1(guess(:,1),guess(:,v+1),g{bid+1}.x,...
                        'linear','extrap');
                end
%                 ni = g{bid+1}.attribute.ni; nj = g{bid+1}.attribute.nj; nk = g{bid+1}.attribute.nk; 
%                 x = sum(sum(sum(g{bid+1}.x))) / numel(g{bid+1}.x);
%                 g{bid+1}.ro = ones(ni,nj,nk) * interp1(guess(:,1),guess(:,2),x,'linear','extrap');
%                 g{bid+1}.rovx = ones(ni,nj,nk) * interp1(guess(:,1),guess(:,3),x,'linear','extrap');
%                 g{bid+1}.rovr = ones(ni,nj,nk) * interp1(guess(:,1),guess(:,4),x,'linear','extrap');
%                 g{bid+1}.rorvt = ones(ni,nj,nk) * interp1(guess(:,1),guess(:,5),x,'linear','extrap');
%                 g{bid+1}.roe = ones(ni,nj,nk) * interp1(guess(:,1),guess(:,6),x,'linear','extrap');
%                 g{bid+1}.trans_dyn_vis = ones(ni,nj,nk) * interp1(guess(:,1),guess(:,7),x,...
%                     'linear','extrap');
%                 g{bid+1}.phi = ones(ni,nj,nk) * interp1(guess(:,1),guess(:,8),x,'linear','extrap');
            end
            
        elseif size(guess,2) == 9
            
            % Extend initial guess
            guess = [guess(1,:) ; guess ; guess(end,:)];
            guess(1,1) = -1e6; guess(end,1) = 1e6;

            % Extract guess parameters from file
            x_guess = guess(:,1); ro_guess = guess(:,3); rovx_guess = guess(:,4);
            rovr_guess = guess(:,5); rorvt_guess = guess(:,6); roe_guess = guess(:,7);
            trans_dyn_vis_guess = guess(:,8); phi_guess = guess(:,9);
            
            % Loop over every block
            for bid = ts_get_bids(g)
                
                % Interpolate primary flow properties
                g{bid+1}.ro = interp1(x_guess,ro_guess,g{bid+1}.x,'linear','extrap');
                g{bid+1}.rovx = interp1(x_guess,rovx_guess,g{bid+1}.x,'linear','extrap');
                g{bid+1}.rovr = interp1(x_guess,rovr_guess,g{bid+1}.x,'linear','extrap');
                g{bid+1}.rorvt = interp1(x_guess,rorvt_guess,g{bid+1}.x,'linear','extrap');
                g{bid+1}.roe = interp1(x_guess,roe_guess,g{bid+1}.x,'linear','extrap');
                g{bid+1}.trans_dyn_vis = interp1(x_guess,trans_dyn_vis_guess,g{bid+1}.x,...
                    'linear','extrap');
                g{bid+1}.phi = interp1(x_guess,phi_guess,g{bid+1}.x,'linear','extrap');                

            end
            
            % Restart poisson loop
            g{bid+1}.av.poisson_restart = 1;
        end
        
    elseif isnumeric(initial_value) == 1 && size(initial_value,2) == 6
        
        % Apply a uniform guess throughout mesh
        varnames = {'ro', 'rovx', 'rovr', 'rorvt', 'roe', 'trans_dyn_vis'}; 
        initial_value = repmat(initial_value,size(g,1),1);
        for bid = ts_get_bids(g)
            for v = 1:length(varnames)
                varname = varnames{v};
                g{bid+1}.(varname) = initial_value(bid+1,v) * ones(g{bid+1}.attribute.ni,...
                    g{bid+1}.attribute.nj, g{bid+1}.attribute.nk);
            end
        end
    else
        
        % Calculate and apply 1D guess to mesh
        g = ts_calc_guess(g);
        
    end

end

end
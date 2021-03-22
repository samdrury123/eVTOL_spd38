function job = bl_create_machine(B,dr,file_prefix,job,I,M,F,D,plot_stuff,monitor_job,run_ag,run_ts,dev)
% BL_CREATE_MACHINE  Create a fully running TURBOSTREAM mesh from a set of blade spline definitions

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Default to monitor jobs
if exist('monitor_job','var') == 0
    monitor_job = 1;
end

% Default to run AUTOGRID
if exist('run_ag','var') == 0
    run_ag = 1;
end

% Default to run TURBOSTREAM
if exist('run_ts','var') == 0
    run_ts = 1;
end

% Default to no deviations
if exist('dev','var') == 0
    dev = cell(length(B),1);
end

% Default to run SA model straight away
if isfield(job,'run_sa') == 0
    job.run_sa = 1;
end

% Open figure windows
if plot_stuff == 1
    h.mer = figure('position',[1 26 1280 895]); hold on; grid on; box on; axis equal; view(2);
    xlabel('Axial Coordinate'); ylabel('Radial Coordinate'); 
end

% Number of Autogrid attempts
n_ag_lim = 1;

% Number of rows
nn = length(B);

% Define colours
C = lines(nn);

% Meridional lines to use
use_F0 = 1;

% Repeat mesh parameters if only specified once
varnames = fieldnames(M); 
for v = 1:length(varnames)
    if numel(M.(varnames{v})) == 1
        M.(varnames{v}) = ones(1,nn) * M.(varnames{v}); 
    end
end


%% Create blade geometries depending on specification

% Create blades and meridional lines if spline parameters exist
for n = 1:nn
    if isfield(B{n},'chi_le') == 1
        B{n} = bl_construct_blade(B{n},0,[],dev{n});
    elseif isfield(B{n},'A_le') == 1
        B{n} = bl_construct_mer(B{n});
    end
end

% Flip blades to run in the positive direction
for n = 1:nn
    
    % Get important indices
    i_te = round(size(B{n}.xrrt,1)/2); j_mid = round(size(B{n}.xrrt,2)/2);
    
    % Check orientation
    if I.rpm(n) > 0 && B{n}.xrrt(1,j_mid,3) < B{n}.xrrt(i_te,j_mid,3)
        
        % Flip rotor direction
        B{n}.xrrt(:,:,3) = - B{n}.xrrt(:,:,3);
        
    elseif I.rpm(n) == 0 && B{n}.xrrt(1,j_mid,3) > B{n}.xrrt(i_te,j_mid,3)
        
        % Flip stator direction
        B{n}.xrrt(:,:,3) = - B{n}.xrrt(:,:,3);

    end
end

% Plot top and bottom blade sections
if plot_stuff == 1
    figure(h.mer)
    for n = 1:nn
        plot3(B{n}.xrrt(:,[1 end],1),B{n}.xrrt(:,[1 end],2),B{n}.xrrt(:,[1 end],3),'-','color',C(n,:))
%         plot3(B{n}.xrrt(:,:,1),B{n}.xrrt(:,:,2),B{n}.xrrt(:,:,3),'-','color',C(n,:))
    end
end

% Generate meridional lines for the whole machine
d_tol = (B{end}.x_out - B{1}.x_in) / 3000; d_tol = max(d_tol,2e-4);
b = bl_construct_mer(B,d_tol); xr_hub = b.xr_hub; xr_cas = b.xr_cas;

% Plot meridional lines
if plot_stuff == 1
    figure(h.mer)
    t = bl_guess_camber(B{1}.xrrt(:,1,[1 3])); z_plot = 0.5*t.chord;
    plot3(xr_hub(:,1),xr_hub(:,2),z_plot*ones(size(xr_hub,1),1),'k-')
    plot3(xr_cas(:,1),xr_cas(:,2),z_plot*ones(size(xr_cas,1),1),'k-')
end

% Add fillets on to blades
if isempty(F) == 0
    for n = 1:nn
        B{n}.xrrt_F0 = B{n}.xrrt;
        if isempty(F{n}) == 0 && any(F{n}.flag) == 1
            B{n}.xrrt = bl_construct_fillet(B{n}.xrrt,F{n},xr_hub,xr_cas);
        end
    end
end

% Plot blade section definitions
if plot_stuff == 1
    figure(h.mer)
    for n = 1:nn
        mesh(B{n}.xrrt(:,:,1),B{n}.xrrt(:,:,2),B{n}.xrrt(:,:,3),...
            repmat(reshape(C(n,:),[1 1 3]),[size(B{n}.xrrt,1) size(B{n}.xrrt,2) 1]));
    end
end

% Generate shroud or bleed stubs if required
xr_stb = cell(nn,1);
for n = 1:nn
    
    % Shrouds at hub LE
    if isfield(D{n},'shroud') == 1 && sum(D{n}.shroud) > 0
        
        % Generate axial coordinates by offset or absolute values
        xr_stb{n} = zeros(2,2,2);
        if isfield(D{n},'dx_shroud') == 1
            xr_stb{n}(:,1,1) = min(B{n}.xrrt(:,1,1)) - D{n}.dx_shroud - D{n}.dx_gap;
            xr_stb{n}(:,2,1) = min(B{n}.xrrt(:,1,1)) - D{n}.dx_shroud;
        else
            xr_stb{n}(:,1,1) = D{n}.x_shroud(1); xr_stb{n}(:,2,1) = D{n}.x_shroud(2); 
        end
        
        % Generate radial coordinates
        xr_stb{n}(1,:,2) = interp1(xr_hub(:,1),xr_hub(:,2),xr_stb{n}(1,:,1),'pchip');
        xr_stb{n}(2,:,2) = xr_stb{n}(1,:,2) - D{n}.dr_shroud;
        
    end
    
    % Bleed slots at casing LE
    if isfield(D{n},'bleed') == 1
        xr_stb{n} = zeros(2,2,2);
        xr_stb{n}(:,1,1) = min(B{n}.xrrt(:,end,1)) - D{n}.dx_cas - D{n}.dx_gap;
        xr_stb{n}(:,2,1) = min(B{n}.xrrt(:,end,1)) - D{n}.dx_cas;
        xr_stb{n}(2,2,1) = xr_stb{n}(2,2,1) + D{n}.dr_bleed * tand(D{n}.bl_angle);
        xr_stb{n}(1,:,2) = interp1(xr_cas(:,1),xr_cas(:,2),xr_stb{n}(1,:,1),'pchip');
        xr_stb{n}(2,:,2) = xr_stb{n}(1,:,2) + D{n}.dr_bleed;
        xr_stb{n}(2,1,1) = xr_stb{n}(2,2,1) - D{n}.A_bleed * D{n}.dx_gap * ...
            xr_stb{n}(1,1,2) / xr_stb{n}(2,1,2);
    end
end

% Plot shroud stubs
if plot_stuff == 1
    figure(h.mer)
    for n = 1:nn
        if isempty(xr_stb{n}) == 0
            plot3(xr_stb{n}(:,:,1),xr_stb{n}(:,:,2),z_plot*ones(size(xr_stb{n}(:,:,1))),'k-')
        end
    end
end


%% Calculate all boundary conditions and interstage lines

% Resolution of axially constant lines
nj = 200; r_nondim = linspace(0,1,nj);

% Preallocate boundary condition line arrays
xr_bc = zeros(nn+1,nj,2);

% Inlet plane
xr_bc(1,:,:) = cat(3,ones(1,nj) * B{1}.x_in,r_nondim * (xr_cas(1,2) - xr_hub(1,2)) + xr_hub(1,2));

% Outlet plane
xr_bc(end,:,:) = cat(3,ones(1,nj) * B{end}.x_out,...
    r_nondim * (xr_cas(end,2) - xr_hub(end,2)) + xr_hub(end,2));

% Choose type of leading and trailing edge to use
if use_F0 == 0; xr_type = 'xrrt'; else xr_type = 'xrrt_F0'; end;

% Find blade leading and trailing edge coordinates
xr_le{n} = cell(nn,1); xr_te{n} = cell(nn,1);
for n = 1:nn
    xr_le{n} = zeros(size(B{n}.(xr_type),2),2); xr_te{n} = zeros(size(B{n}.(xr_type),2),2);
    for j = 1:size(B{n}.(xr_type),2)
        [~,i_le] = min(B{n}.(xr_type)(:,j,1)); [~,i_te] = max(B{n}.(xr_type)(:,j,1)); 
        xr_le{n}(j,:) = B{n}.(xr_type)(i_le,j,1:2); xr_te{n}(j,:) = B{n}.(xr_type)(i_te,j,1:2);
    end
end

% Non-dimensional leading and trailing edge coordinates
x_le = zeros(nn,nj); x_te = zeros(nn,nj);
for n = 1:nn
    
    % Initial coarse interpolation
    r_coarse = linspace(0.02,0.98,9);
    x_le_temp = interp1((xr_le{n}(:,2) - xr_le{n}(1,2)) / (xr_le{n}(end,2) - xr_le{n}(1,2)),...
        xr_le{n}(:,1),r_coarse,'pchip');
    x_te_temp = interp1((xr_te{n}(:,2) - xr_te{n}(1,2)) / (xr_te{n}(end,2) - xr_te{n}(1,2)),...
        xr_te{n}(:,1),r_coarse,'pchip');
    
    % Fine interpolation to smooth effect of fillets
    x_le(n,:) = interp1(r_coarse,x_le_temp,r_nondim,'pchip','extrap');
    x_te(n,:) = interp1(r_coarse,x_te_temp,r_nondim,'pchip','extrap');
    
end

% Mixing planes midway between leading and trailing edges
xr_bc(2:end-1,:,1) = 0.5 * (x_le(2:end,:) + x_te(1:end-1,:));

% Move mixing planes if constant coordinates are specified
for n = 1:nn-1
    if isfield(B{n},'x_mp') == 1
        xr_bc(n+1,:,1) = interp1(linspace(0,1,length(B{n}.x_mp)),B{n}.x_mp,r_nondim,'spline');
    end
end

% Move mixing planes in the case of shrouded or bleed slot rows
for n = 2:nn
    if isfield(D{n},'shroud') == 1 && sum(D{n}.shroud) > 0
        if D{n}.shroud(1) ~= 0 && xr_bc(n,1,1) > xr_stb{n}(1,1,1)
            dx = interp1([0 0.5 1],[xr_bc(n,1,1) - xr_stb{n}(1,1,1) + D{n}.dx_mp 0 0],...
                r_nondim,'pchip');
            xr_bc(n,:,1) = xr_bc(n,:,1) - dx;
        end
    end
    if isfield(D{n},'bleed') == 1
        if xr_bc(n,end,1) > xr_stb{n}(1,1,1)
            dx = interp1([0 0.5 1],[0 0 xr_bc(n,end,1) - xr_stb{n}(1,1,1) + D{n}.dx_mp],...
                r_nondim,'pchip');
            xr_bc(n,:,1) = xr_bc(n,:,1) - dx;
        end
    end    
end

% Row leading and trailing edge z-constant lines by offsetting blade edges
xr_row = cell(nn,1);
for n = 1:nn
    xr_row{n} = zeros(2,nj,2);
    xr_row{n}(1,:,1) = x_le(n,:) - M.dx_row(n);
    xr_row{n}(2,:,1) = x_te(n,:) + M.dx_row(n);
end

% Fit trailing edge z-constant lines within the domain
for n = 1:nn-1
    if sum(xr_row{n}(2,:,1) > xr_bc(n+1,:,1)) > 0
        p = (xr_row{n}(2,:,1) - x_te(n,:)) ./ (xr_bc(n+1,:,1) - x_te(n,:));
        p = mean(p);
        xr_row{n}(2,:,1) = p * (xr_bc(n+1,:,1) - x_te(n,:)) + x_te(n,:);
    end
end

% Move and create leading edge z-constant lines in the case of shrouded or bleed slot rows
for n = 1:nn
    if isfield(D{n},'shroud') == 1 && sum(D{n}.shroud) > 0
        
        % Move leading edge line if overlapping the mixing plane
        dx = xr_row{n}(1,:,1) - xr_bc(n,:,1) - D{n}.dx_mp;
        dx(dx > 0) = 0; dx = smooth(dx)';
        xr_row{n}(1,:,1) = xr_row{n}(1,:,1) - dx;
        
        % Match the leading edge line
        dx = xr_row{n}(1,1,1) - xr_stb{n}(1,2,1);
        dx = interp1([0 0.2 1],[dx 0 0],r_nondim,'pchip');
        xr_row{n}(1,:,1) = xr_row{n}(1,:,1) - dx;
        
        % Create a new line on upstream edge of stub
        p = (xr_stb{n}(1,1,1) - xr_stb{n}(1,2,1)) / (xr_bc(n,1,1) - xr_stb{n}(1,2,1));
        x = p * (xr_bc(n,:,1) - xr_row{n}(1,:,1)) + xr_row{n}(1,:,1);
        xr_row{n} = [cat(3,x,zeros(1,nj,1)) ; xr_row{n}];
    end
    if isfield(D{n},'bleed') == 1
                
        % Move leading edge line if overlapping the mixing plane
        dx = xr_row{n}(1,:,1) - xr_bc(n,:,1) - D{n}.dx_mp;
        dx(dx > 0) = 0; dx = smooth(dx)';
        xr_row{n}(1,:,1) = xr_row{n}(1,:,1) - dx;
        
        % Match the leading edge line
        dx = xr_row{n}(1,end,1) - xr_stb{n}(1,2,1);
        dx = interp1([0 0.8 1],[0 0 dx],r_nondim,'pchip');
        xr_row{n}(1,:,1) = xr_row{n}(1,:,1) - dx;
        
        % Create a new line on upstream edge of stub
        p = (xr_stb{n}(1,1,1) - xr_stb{n}(1,2,1)) / (xr_bc(n,end,1) - xr_stb{n}(1,2,1));
        x = p * (xr_bc(n,:,1) - xr_row{n}(1,:,1)) + xr_row{n}(1,:,1);
        xr_row{n} = [cat(3,x,zeros(1,nj,1)) ; xr_row{n}];
        
    end
end
    
% Interpolate radial mixing plane coordinates
xr_bc(:,1,2) = interp1(xr_hub(:,1),xr_hub(:,2),xr_bc(:,1,1),'pchip');
xr_bc(:,end,2) = interp1(xr_cas(:,1),xr_cas(:,2),xr_bc(:,end,1),'pchip');
xr_bc(:,:,2) = repmat(r_nondim,[nn+1 1]) .* repmat((xr_bc(:,end,2) - xr_bc(:,1,2)),[1 nj]) +...
    repmat(xr_bc(:,1,2),[1 nj]);
    
% Interpolate radial z-constant line coordinates
for n = 1:nn
    xr_row{n}(:,1,2) = interp1(xr_hub(:,1),xr_hub(:,2),xr_row{n}(:,1,1),'pchip');
    xr_row{n}(:,end,2) = interp1(xr_cas(:,1),xr_cas(:,2),xr_row{n}(:,end,1),'pchip');
    xr_row{n}(:,:,2) = repmat(r_nondim,[size(xr_row{n},1) 1]) .* repmat((xr_row{n}(:,end,2) - ...
        xr_row{n}(:,1,2)),[1 nj]) + repmat(xr_row{n}(:,1,2),[1 nj]);
end

% Plot all axially constant lines
if plot_stuff == 1
    figure(h.mer)
    plot3(xr_bc(:,:,1)',xr_bc(:,:,2)',z_plot*ones(nj,nn+1,1),'-','color',[0.3 0.3 0.3])
    for n = 1:nn
        plot3(xr_row{n}(:,:,1)',xr_row{n}(:,:,2)',...
            z_plot*ones(nj,size(xr_row{n},1),1),'-','color',[0.5 0.5 0.5])
    end
    
    % Update axis limits
    v = axis; axis([B{1}.x_in - M.dx_row(1) B{end}.x_out + M.dx_row(1) v(3:4)]);
end


%% Write input files to AUTOGRID

% Write blade coordinate files
for n = 1:nn
    b = B{n}; xyz = bl_pol2cart(b.xrrt);
    ag_write_dat(dr.ag,[file_prefix 'B' num2str(n) '.dat'],xyz(:,:,[2 3 1]));
    save([dr.geom file_prefix 'B' num2str(n) '_spline.mat'],'b')
end

% Write meridional curves
ag_write_dat(dr.ag,[file_prefix 'xr_hub.dat'],xr_hub(:,[2 1]))
ag_write_dat(dr.ag,[file_prefix 'xr_cas.dat'],xr_cas(:,[2 1]))
save([dr.geom file_prefix 'xr_hub.mat'],'xr_hub')
save([dr.geom file_prefix 'xr_cas.mat'],'xr_cas')

% Write boundary condition planes
for n = 1:nn+1
    ag_write_dat(dr.ag,[file_prefix 'xr_bc_' num2str(n) '.dat'],squeeze(xr_bc(n,:,[2 1])))
end
save([dr.geom file_prefix 'xr_bc.mat'],'xr_bc')

% Write other axially constant curves
for n = 1:nn
    for m = 1:size(xr_row{n},1)
        ag_write_dat(dr.ag,[file_prefix 'xr_row_' num2str(n) '_' num2str(m) '.dat'],...
            squeeze(xr_row{n}(m,:,[2 1])))
    end
end
save([dr.geom file_prefix 'xr_row.mat'],'xr_row')

% Write shroud stubs
for n = 1:nn
    if isempty(xr_stb{n}) == 0
        ag_write_dat(dr.ag,[file_prefix 'xr_stb_B' num2str(n) '_1.dat'],squeeze(xr_stb{n}(:,1,[2 1])));
        ag_write_dat(dr.ag,[file_prefix 'xr_stb_B' num2str(n) '_2.dat'],squeeze(xr_stb{n}(:,2,[2 1])));
    end
end
save([dr.geom file_prefix 'xr_stb.mat'],'xr_stb')

% Assemble configuration variables into a single structure
c.s_hub = zeros(1,nn); c.s_cas = zeros(1,nn); c.d_hub = zeros(1,nn); c.d_cas = zeros(1,nn);
c.f_hub = zeros(1,nn); c.f_cas = zeros(1,nn); 
for n = 1:nn; c.N(n) = B{n}.N; end;
for n = 1:nn; if isfield(D{n},'hub') == 1 && D{n}.flag(1) ~= 0; c.d_hub(n) = D{n}.hub; end; end;
for n = 1:nn; if isfield(D{n},'cas') == 1 && D{n}.flag(2) ~= 0; c.d_cas(n) = D{n}.cas; end; end;
for n = 1:nn
    if isfield(D{n},'shroud') == 1
        c.s_hub(n) = D{n}.shroud(1); c.s_cas(n) = D{n}.shroud(2); 
    end
    if isfield(D{n},'bleed') == 1; c.s_cas(n) = 1; end
end
for n = 1:nn; if isfield(F{n},'hub') == 1 && F{n}.flag(1) ~= 0; 
        c.f_hub(n) = min(F{n}.hub.rad1); end; end;
for n = 1:nn; if isfield(F{n},'cas') == 1 && F{n}.flag(2) ~= 0; 
        c.f_cas(n) = min(F{n}.cas.rad1); end; end;
c.rpm = I.rpm;

% Calculate z-constant spacings for flow path control
for n = 1:nn
    
    % Default spacings
    c.o_le_hub(n) = 0.02; c.o_le_cas(n) = 0.98; c.o_te_hub(n) = 0.02; c.o_te_cas(n) = 0.98;
    c.e_le_hub(n) = M.e_f(n); c.e_le_cas(n) = M.e_f(n); 
    c.e_te_hub(n) = M.e_f(n); c.e_te_cas(n) = M.e_f(n);
    
    % Leading and trailing edge spacings for fillets
    if F{n}.flag(1) ~= 0
        c.o_le_hub(n) = c.f_hub(n) / (diff(xr_le{n}([1 end],2))); 
        c.o_te_hub(n) = c.f_hub(n) / (diff(xr_te{n}([1 end],2)));
    end
    if F{n}.flag(2) ~= 0
        c.o_le_cas(n) = c.f_cas(n) / (diff(xr_le{n}([1 end],2))); 
        c.o_te_cas(n) = c.f_cas(n) / (diff(xr_te{n}([1 end],2))); 
    end

    % Leading and trailing edge spacings for gaps
    if D{n}.flag(1) == 1
        c.e_le_hub(n) = M.e_j(n); c.e_te_hub(n) = M.e_j(n); 
        c.o_le_hub(n) = c.d_hub(n) / (diff(xr_le{n}([1 end],2))); 
        c.o_te_hub(n) = c.d_hub(n) / (diff(xr_te{n}([1 end],2))); 
    end
    if D{n}.flag(2) == 1
        c.e_le_cas(n) = M.e_j(n); c.e_te_cas(n) = M.e_j(n); 
        c.o_le_cas(n) = c.d_cas(n) / (diff(xr_le{n}([1 end],2))); 
        c.o_te_cas(n) = c.d_cas(n) / (diff(xr_te{n}([1 end],2))); 
    end

end

% Spacings at inlet and outlet of the machine
c.o_in_hub(1) = 0.04; c.o_in_cas(1) = 0.04;
c.o_out_hub(nn) = 0.04; c.o_out_cas(nn) = 0.04;
c.e_in_hub(1) = M.e_f(1); c.e_in_cas(1) = M.e_f(1); 
c.e_out_hub(nn) = M.e_f(end); c.e_out_cas(nn) = M.e_f(end); 

% Calculate spacings on mixing planes and shroud lines
fac = 1;
for n = 1:nn
    
    % Mixing planes in intermediate rows
    if n < nn
        if length(xr_row{n}) == 3; i_le = 2; else i_le = 1; end;
        p = (xr_bc(n+1,1,1) - xr_row{n}(end,1,1)) / (xr_row{n+1}(i_le,1,1) - xr_row{n}(end,1,1));
        c.e_out_hub(n) = (p * c.e_le_hub(n+1)^fac + (1-p) * c.e_te_hub(n)^fac)^fac;
        c.o_out_hub(n) = p * (c.o_le_hub(n+1) - c.o_te_hub(n)) + c.o_te_hub(n);
        p = (xr_bc(n+1,end,1) - xr_row{n}(end,end,1)) / (xr_row{n+1}(i_le,end,1) - xr_row{n}(end,end,1));
        c.e_out_cas(n) = (p * c.e_le_cas(n+1)^fac + (1-p) * c.e_te_cas(n)^fac)^fac;
        c.o_out_cas(n) = p * (c.o_le_cas(n+1) - c.o_te_cas(n)) + c.o_te_cas(n);
        c.e_in_hub(n+1) = c.e_out_hub(n); c.o_in_hub(n+1) = c.o_out_hub(n);
        c.e_in_cas(n+1) = c.e_out_cas(n); c.o_in_cas(n+1) = c.o_out_cas(n);
    end
    
    % Z constant line by inlet if shrouded
    if size(xr_row{n},1) == 3
        p = (xr_row{n}(1,1,1) - xr_bc(n,1,1)) / (xr_row{n}(2,1,1) - xr_bc(n,1,1));
        c.e_sh_hub(n) = (p * c.e_le_hub(n)^fac + (1-p) * c.e_in_hub(n)^fac)^fac;
        c.o_sh_hub(n) = p * (c.o_le_hub(n) - c.o_in_hub(n)) + c.o_in_hub(n);
        p = (xr_row{n}(1,end,1) - xr_bc(n,end,1)) / (xr_row{n}(2,end,1) - xr_bc(n,end,1));
        c.e_sh_cas(n) = (p * c.e_le_cas(n)^fac + (1-p) * c.e_in_cas(n)^fac)^fac;
        c.o_sh_cas(n) = p * (c.o_le_cas(n) - c.o_in_cas(n)) + c.o_in_cas(n);
    end
end

% Determine whether to use butterfly fillet topologies
for n = 1:nn
    if F{n}.flag(1) == 1; c.f_hub_fly(n) = 1; else c.f_hub_fly(n) = 0; end;
    if F{n}.flag(2) == 1; c.f_cas_fly(n) = 1; else c.f_cas_fly(n) = 0; end;
end

% Determine whether to use high or low staggered topologies, default to compressor
if isfield(c,'stin') == 0 && isfield(c,'stout') == 0
    for n = 1:nn
        if c.rpm(n) == 0
            c.stin(n) = -1; c.stout(n) = 0;
        else
            c.stin(n) = 1; c.stout(n) = 0;
        end
    end
end

% Default to untwist inlet to first row and outlet to last - inlet number measured from blade LE
c.ui = zeros(1,nn); c.ui(1) = 0.5; c.uo = zeros(1,nn); c.uo(end) = 0.5;

% Calculate gap interpolation spanwise extent
c.o_hub = zeros(1,nn); c.o_cas = zeros(1,nn);
for n = 1:nn
    if D{n}.flag(1) == 1
        x = mean(B{n}.xrrt(:,1,1)); 
        c.o_hub(n) = 100 * (c.d_hub(n) / (interp1(xr_cas(:,1),xr_cas(:,2),x) - ...
            interp1(xr_hub(:,1),xr_hub(:,2),x)) + 0.01);
    elseif D{n}.flag(2) == 1
        x = mean(B{n}.xrrt(:,end,1)); 
        c.o_cas(n) = 100 * (1 - c.d_cas(n) / (interp1(xr_cas(:,1),xr_cas(:,2),x) - ...
            interp1(xr_hub(:,1),xr_hub(:,2),x)) - 0.01);        
    end
end     

% Default to calculate expansion ratio in AUTOGRID
if isfield(M,'r_b') == 0
    c.r_b = zeros(1,nn);
end

% Default 2% spanwise interpolation
if isfield(M,'i_j') == 0
    c.i_j = 2 * ones(1,nn);
end

% Default free outlet clustering
if isfield(M,'e_o') == 0
    c.e_o = zeros(1,nn);
end

% Assemble mesh parameters into configuration structure
varnames = fieldnames(M);
for v = 1:length(varnames)
    c.(varnames{v}) = M.(varnames{v});
end

% Write configuration text file
fid = fopen([dr.ag file_prefix 'config.dat'],'w');
varnames = fieldnames(c);
for v = 1:length(varnames)
    fprintf(fid,'%s',varnames{v}); 
    fprintf(fid,'\t%g',c.(varnames{v})); 
    fprintf(fid,'\n');
end
fclose(fid);


%% Run AUTOGRID to generate mesh

% Delete existing autogrid files
if exist([dr.ag file_prefix(1:end-1) '.g'],'file') ~= 0 && run_ag == 1
    delete([dr.ag file_prefix(1:end-1) '.g']);
    delete([dr.ag file_prefix(1:end-1) '.trb']);
end

% Loop over until plot3d file is created in case autogrid license is not available
n_ag = 1;
while exist([dr.ag file_prefix(1:end-1) '.g'],'file') == 0 && run_ag == 1

    % Export directory and file_prefix to bash variables
    setenv('directory',dr.ag);
    setenv('file_prefix',file_prefix);

    % Run AUTOGRID to generate mesh
    if run_ag == 1
        fprintf(['%% Running AUTOGRID on ' file_prefix(1:end-1) ' %%\n'])
        [~,~] = system('/opt/numeca/bin/igg -batch -autogrid5 -script ./AG/ag_create_machine.py');
%         [~,~] = system('/opt/numeca/bin/igg -autogrid5 -script ./AG/ag_create_machine.py');
    end
    
    % Check grid quality
    if exist([dr.ag file_prefix(1:end-1) '.qualityReport'],'file') ~= 0
        q = ag_quality_check([dr.ag file_prefix(1:end-1) '.qualityReport']);
        if q.skew < M.skew(1) || q.neg == 1
            disp(q);
%             error('%% Poor cells found %%');
%             job = [];
%             return;
            quality_pass = 0;
            n_ag = n_ag + 1;
        else
            quality_pass = 1;
        end
    end
    
    % Break the loop if a decent mesh can't be made
    if n_ag > n_ag_lim
        error('Unable to mesh below skewness limit');
    end

    % Run IGG to re-cluster and convert mesh
    if exist([dr.ag file_prefix(1:end-1) '.trb'],'file') ~= 0 && quality_pass == 1
        setenv('directory',dr.ag);
        setenv('file_prefix',file_prefix);
        [~,~] = system('/opt/numeca/bin/igg -batch -script ./AG/ag_igg_machine.py');
    end
    
    % Check if AUTOGRID has completed successfully, if not restretch blades
    if exist([dr.ag file_prefix(1:end-1) '.g'],'file') == 0
        
        % Pause to wait for license
        pause(60);
        
        % Re-write blade files with stretched coordinates
        for n = 1:nn
            
            % Calculate random stretch array
            xrrt = B{n}.xrrt; [ni,nj,~] = size(xrrt);
            tol = rand * 5 * F{n}.tol - F{n}.tol;
            xrrt(:,:,2) = xrrt(:,:,2) + ...
                repmat([linspace(-tol,0,round(nj/2)) linspace(0,tol,nj-round(nj/2))],[ni 1 1]);
            
            % Re-write AUTOGRID blade input files
            xyz = bl_pol2cart(xrrt);
            ag_write_dat(dr.ag,[file_prefix 'B' num2str(n) '.dat'],xyz(:,:,[2 3 1]));
        end
    end

end


%% Convert mesh to HDF5 TURBOSTREAM format on Darwin HPC

% Copy AUTOGRID files onto HPC
[~, ~] = system(['scp ' dr.ag file_prefix(1:end-1) '.bcs ' job.user '@' job.host ':' dr.hpc '.']);
[~, ~] = system(['scp ' dr.ag file_prefix(1:end-1) '.g ' job.user '@' job.host ':' dr.hpc '.']);

% Copy PYTHON conversion function to HPC
[~, ~] = system(['scp ./AG/ag_convert_ts.py ' job.user '@' job.host ':~/.']);

% Run PYTHON conversion script on HPC
if strcmp(job.arch,'pascal') == 1
    [~, ~] = system(['ssh ' job.user '@' job.host ' "source /rds/project/hpc/rds-hpc-pullan/ts3/ts' ...
        strrep(job.version,'.','') '/bashrc_module_ts' strrep(job.version,'.','') ...
        ' && python ~/ag_convert_ts.py ' dr.hpc file_prefix(1:end-1) '"']);
else
    [~, ~] = system(['ssh ' job.user '@' job.host ' "module load turbostream' job.version ...
        ' && python ~/ag_convert_ts.py ' dr.hpc file_prefix(1:end-1) '"']);
end

% Copy back and rename
[~, ~] = system(['scp ' job.user '@' job.host ':' dr.hpc file_prefix 'input.hdf5 ' dr.ag '.']);

% Delete autogrid and TURBOSTREAM files on HPC
[~, ~] = system(['ssh ' job.user '@' job.host ' rm ' dr.hpc file_prefix(1:end-1) '.bcs']);
[~, ~] = system(['ssh ' job.user '@' job.host ' rm ' dr.hpc file_prefix(1:end-1) '.g']);
[~, ~] = system(['ssh ' job.user '@' job.host ' rm ' dr.hpc file_prefix(1:end-1) '_input.hdf5']);


%% Pre-process mesh in MATLAB for running in TURBOSTREAM

% Check whether to run TURBOSTREM or not
if run_ts == 1

% Read in converted file
g = ts_read_hdf5([dr.ag file_prefix 'input.hdf5']);

% Delete external patches
for bid = ts_get_bids(g)
    for pid = ts_get_pids(g,bid)
        if g{bid+1}.patch{pid+1}.attribute.kind == 13
            g = ts_remove_patch(g,bid,pid);            
            break
        end
    end
end

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

% Find bids of stubs to add shroud leakage model to
bids = ts_get_row_bids(g,xr_hub,xr_cas); bid_stb = zeros(nn,1); row_stb = zeros(nn,1);
for n = 1:nn
    if isfield(D{n},'shroud') == 1 && sum(D{n}.shroud) > 0 && isfield(D{n},'gap') == 1
        for bid = bids{n}
            r = mean(g{bid+1}.r(:)); x = mean(g{bid+1}.x(:));
            if r < interp1(xr_hub(:,1),xr_hub(:,2),x)
                bid_stb(n) = bid;
            end
        end
        row_stb(n) = D{n}.shroud(1);
    end
end

% Add bleed patches
clear p;
for n = 1:nn
    if isfield(D{n},'bleed') == 1
        for bid = bids{n}
            if mean(g{bid+1}.r(:)) > interp1(xr_cas(:,1),xr_cas(:,2),mean(g{bid+1}.x(:)))
                p.attribute.bid = bid; p.attribute.pid = length(g{bid+1}.patch); 
                p.attribute.ist = 0; p.attribute.ien = g{bid+1}.attribute.ni;
                p.attribute.jst = 0; p.attribute.jen = 1;
                p.attribute.kst = 0; p.attribute.ken = g{bid+1}.attribute.nk;
                p.bleed_flow = D{n}.bleed; p.attribute.kind = 10;
                g{bid+1}.patch{end+1} = p; g{bid+1}.attribute.np = length(g{bid+1}.patch);    
            end
        end
    end
end

% Match up pairs of blocks to add shroud leakage model
bids = []; ns = [];
for n = 1:nn
    if bid_stb(n) ~= 0
        n1 = row_stb(n); n2 = row_stb(row_stb(n));
        bids = [bids ; bid_stb(n1) bid_stb(n2)];
        ns = [ns ; n1 n2];
    end
end
bids = unique(sort(bids,2),'rows');

% Add on shroud leakage patches
for n = 1:size(bids,1)
    
    % Initialise upstream patch
    bid = bids(n,1);
    p0.bid = bid; p0.pid = length(g{bid+1}.patch); p0.kind = 21;
    p0.ist = 0; p0.ien = g{bid+1}.attribute.ni;
    p0.jst = 0; p0.jen = 1;
    p0.kst = 0; p0.ken = g{bid+1}.attribute.nk;

    % Initialise downstream patch
    bid = bids(n,2);
    p1.bid = bid; p1.pid = length(g{bid+1}.patch); p1.kind = 21;
    p1.ist = 0; p1.ien = g{bid+1}.attribute.ni;
    p1.jst = 0; p1.jen = 1;
    p1.kst = 0; p1.ken = g{bid+1}.attribute.nk;

    % Connect patches
    p0.nxbid = p1.bid; p0.nxpid = p1.pid;
    p1.nxbid = p0.bid; p1.nxpid = p0.pid;

    % Add patches to mesh
    g{bids(n,1)+1}.patch{end+1}.attribute = p0;
    g{bids(n,2)+1}.patch{end+1}.attribute = p1;
    g{bids(n,1)+1}.attribute.np = length(g{bids(n,1)+1}.patch);
    g{bids(n,2)+1}.attribute.np = length(g{bids(n,2)+1}.patch);

    % Set patch variables
    for bid = bids(n,:)
        g{bid+1}.patch{end}.shroud_sealgap = D{ns(n,1)}.gap;    % Size of clearance
        g{bid+1}.patch{end}.shroud_nseal = 2;                   % Number of fins
        g{bid+1}.patch{end}.shroud_cfshroud = 0.0025;           % Skin friction underside of platform
        g{bid+1}.patch{end}.shroud_cfcasing = 0.0025;           % Skin friction on endwall 
        g{bid+1}.patch{end}.shroud_wshroud = 0;                 % Rotation of underside of platform
        g{bid+1}.patch{end}.shroud_wcase = max(I.rpm);          % Rotation of endwall
        g{bid+1}.patch{end}.shroud_pitchin = 80;                % Injection angle
    end
    g{bids(n,1)+1}.patch{end}.shroud_dir = 1;                   % Re-injection patch
    g{bids(n,2)+1}.patch{end}.shroud_dir = 0;                   % Ingestion patch
    
end

% Generate rotating and stationary meridional coordinates for rotational face setting
xr_rotating = []; xr_stationary = [];
for n = 1:nn
    if isfield(D{n},'shroud') == 1 && sum(D{n}.shroud) > 0
        if I.rpm(n) == 0
            xr_rotating = [xr_rotating ; squeeze(xr_stb{n}(:,1,:))];
            xr_stationary = [xr_stationary ; squeeze(xr_stb{n}(:,2,:))];
        else
            xr_rotating = [xr_rotating ; squeeze(xr_stb{n}(:,2,:))];
            xr_stationary = [xr_stationary ; squeeze(xr_stb{n}(:,1,:))];            
        end
    end
end

% Match periodic patches across the pitch in polar coordinates
% g = ts_match_periodics(g);
% g = ts_match_patches(g);

% Save patch data to file
ts_check_patches(g,2*plot_stuff,[dr.ts file_prefix 'patch.dat'])

% Apply initial settings
g = ts_create_initial(g,15000,I.poisson_nstep,I.guess);
    
% Loop over all rows and set rotor or stator rows
for n = 1:nn
    if I.rpm(n) ~= 0
        if D{n}.flag(2) == 0; rows(n) = 1.5; else; rows(n) = 1; end
    else
        if D{n}.flag(1) == 0; rows(n) = 0; else; rows(n) = 0.5; end
    end
end

% Set rpms on blade rows and then cavities
if exist('xr_rotating','var') ~= 0 && isempty(xr_rotating) == 0
    g = ts_set_rpms(g,rows,I.rpm,xr_hub,xr_cas,xr_stationary,xr_rotating);
else
    g = ts_set_rpms(g,rows,I.rpm,xr_hub,xr_cas);
end

% Set stationary faces at bottom of shroud wells
for bid = ts_get_bids(g)
    x = sum(sum(sum(g{bid+1}.x))) / numel(g{bid+1}.x);
    r = sum(sum(sum(g{bid+1}.r))) / numel(g{bid+1}.r);
    if r < interp1(xr_hub(:,1),xr_hub(:,2),x)
        g{bid+1}.bv.rpmj1 = 0;
    end
end

% Add periodics on hub and casing for 2D calculations and set rotational speed to zero
if isfield(job,'itwd') == 1 && job.itwd == 1
    for bid = ts_get_bids(g)
        Q1 = [0 g{bid+1}.attribute.ni ; g{bid+1}.attribute.nj-1 g{bid+1}.attribute.nj ; ...
            0 g{bid+1}.attribute.nk];
        Q2 = [0 g{bid+1}.attribute.ni ; 0 1 ; 0 g{bid+1}.attribute.nk];
        g = ts_add_periodic_patch(bid,bid,Q1,Q2,[0 1 2],[0 1 2],g,5);
        g{bid+1}.bv.rpmj1 = 0; g{bid+1}.bv.rpmj2 = 0;
    end
end
    
% Display rotating surfaces
if plot_stuff == 1
    ts_check_rpms(g);
end

% Apply inlet boundary conditions and rotational speed
g = ts_apply_inlet(g,I.inlet,0);

% Apply restart SA settings
if job.run_sa == 1
    g = ts_create_restart(g,I.nstep,2);
    g{1}.av.poisson_nstep = I.poisson_nstep;
end

% Set rfmix for high speed machines
if isfield(I,'rfmix') == 1
    g{1}.av.rfmix = I.rfmix;
end

% Set number of steps for 2D calculations
% if isfield(job,'itwd') == 1 && job.itwd == 1; g{1}.av.nstep = 0; end;

% Linear interpolation for coarse tip gap meshes
if min(min(M.n_bu),min(M.n_bl)) <= 61
    g{1}.av.nomatch_int = 0;
    disp('% Linear Non-Matching Patch Interpolation %');
end

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

% Set all application and block variables if a file is specified
if isfield(I,'vs') == 1; g = ts_set_vs(g,I.vs); g{1}.av.poisson_nstep = I.poisson_nstep; end;

% Add probe patches to all mixing planes
% g = ts_conv_probes(g,0,1,xr_hub,xr_cas);


%% Run file in TURBOSTREAM

% Decide whether to run SA first or Poisson, then ML, then SA
if job.run_sa == 0
    
    % Run poisson loop and coarse mixing length
    job.name = [file_prefix 'poisson'];
    job = ts_smonitor(ts_sbatch(job,g),0,1);

    % Run fine mixing length solution
    if isnan(job.finished) == 0
        g = ts_read_hdf5([job.directory job.name '.hdf5']);
        g = ts_smooth_poisson(g);
        g = ts_create_restart(g,15000,1);
        
        if isfield(I,'rfmix') == 1
            g{1}.av.rfmix = I.rfmix;
        end
        
        if isfield(I,'m') == 1
            [~,~,dpout] = ts_plot_conv([job.directory job.logname],0);
            g{bid_out+1}.patch{pid_out+1}.pout = g{bid_out+1}.patch{pid_out+1}.pout + dpout;
        end
        
        job.name = [file_prefix 'ML'];
        job = ts_smonitor(ts_sbatch(job,g)); 
    end

    % Run SA restart
    if isnan(job.finished) == 0
        g = ts_read_hdf5([job.directory job.name '.hdf5']);
        g = ts_create_restart(g,I.nstep,2);
        
        if isfield(I,'rfmix') == 1
            g{1}.av.rfmix = I.rfmix;
        end
        
        if isfield(I,'m') == 1
            [~,~,dpout] = ts_plot_conv([job.directory job.logname],0);
            g{bid_out+1}.patch{pid_out+1}.pout = g{bid_out+1}.patch{pid_out+1}.pout + dpout;
        end
        
        job.name = file_prefix(1:end-1);
        job = ts_smonitor(ts_sbatch(job,g));
    end
    
else
    
    % Run SA model straight away
    job.name = file_prefix(1:end-1);
    job = ts_sbatch(job,g); 
    
    % Monitor job or exit function
    if monitor_job == 1; job = ts_smonitor(job); end
    
    % Restart job with slipwalls instead of periodics
%     if isfield(job,'itwd') == 1 && job.itwd == 1
%         copyfile([job.directory job.logname],[job.directory strrep(job.logname,'.txt','_poisson.txt')]);
%         g = ts_read_hdf5([job.directory job.name '.hdf5']);
%         for bid = ts_get_bids(g)
%             g{bid+1}.patch{end}.attribute.kind = 7; 
%             g{bid+1}.patch{end-1}.attribute.kind = 7; 
%         end
%         g{1}.av.nstep = I.nstep; g{1}.av.poisson_nstep = 0; g{1}.av.poisson_restart = 1;
%         g{1}.av.ilos = 1;
%         job = ts_smonitor(ts_sbatch(job,g));
%     end
end

end

end
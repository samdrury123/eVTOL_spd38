%   Diffuser_Example
%                               
%   Generate an example stator blade in a diffuser to run in TURBOSTREAM
%
%   Written by James Taylor              
%   May 2020

% Clear workspace and load eVTOL environment
clear; close all; clc;
eVTOL

% Current job variables
job.priority = 'MILLER-SL2-GPU';
job.cards = 1;
job.walltime = 0.1;

% Options for using interactive node and command to activate on Wilkes
% sintr -p pascal -t 8:00:00 -A MILLER-SL2-GPU --gres=gpu:1 -N 1
% job.node = 'gpu-e-74';

% Mesh parameters
dx_duct = 0.2; tol = 2e-4;
ni_in = 13; ni_out = 21; ni_bl = 45; ni_sp = 13; rj = 1.5; nk = 41;
ek = 5e-5; ej = 7e-5; ei = 3e-4;

% Blade metal angles and count
chi_le = 25; chi_te = -5; N = 15;

% Boundary conditions
Po_in = 101325; To_in = 288.15; Alpha_in = 23; Vx_in = 80; d_hub = 0.1; d_cas = 0.1;

% Number of iterations to run
nstep = 30000; nstep_poisson = 2000; nstep_avg = 5000;

% Reduced number of iterations for faster convergence in TBLOCK
nstep = 5000;


%% Create the geometry

% Load stator example
load([dr.geom 'NRF_3b_splines.mat']); b = B{2};

% Increase 3D stacking
b.sweep.coefs = b.sweep.coefs * 1.5;
b.lean.coefs = b.lean.coefs * 2;

% Modify hub and casing lines
b.x_in = 0.1; b.x_out = 0.35; b.A_te = 1.1;
b.drdx_hub_le = -0.15; b.drdx_hub_te = -0.15;
b.drdx_cas_le = 0.15; b.drdx_cas_te = 0.15;

% Apply blade angles and count
b.chi_le = spap2(1,2,[0 1],[1 1] * chi_le);
b.chi_te = spap2(1,2,[0 1],[1 1] * chi_te);
b.N = N; dt = 2 * pi / N;

% Create and plot the blade
b = bl_construct_blade(b,0);

% Create the meridional lines and extend with constant area ducts
m = bl_construct_mer({b},tol); xr_hub = m.xr_hub; xr_cas = m.xr_cas;
x_in = linspace(b.x_in-dx_duct,b.x_in,dx_duct/tol)'; 
x_out = linspace(b.x_out,b.x_out+dx_duct,dx_duct/tol)';
xr_hub = [x_in xr_hub(1,2) * ones(size(x_in)) ; xr_hub(2:end-1,:) ; ...
    x_out xr_hub(end,2) * ones(size(x_out))];
xr_cas = [x_in xr_cas(1,2) * ones(size(x_in)) ; xr_cas(2:end-1,:) ; ...
    x_out xr_cas(end,2) * ones(size(x_out))];
b.xr_hub = xr_hub; b.xr_cas = xr_cas;

% Initialise splitter cell array
b.split = cell(0,1);

% Struct for splitter definition
b.split{1}.r_nondim = 0.3; b.split{1}.x_nondim = [0.03 0.96];

% Copy section parameters
varnames = {'dcam_le','dcam_te','qcam','thick_max','rad_le','s_thick_max',...
    'rad_thick_max','wedge_te','thick_te'};
for v = 1:length(varnames); b.split{1}.(varnames{v}) = b.(varnames{v}); end;

% Splitter metal angles - must have same angle at both ends
b.split{1}.chi_le = spap2(1,3,[0 0.5 1],[-10 -5 -10]);
b.split{1}.chi_te = spap2(1,3,[0 0.5 1],[-20 -20 -20]);

% Splitter 3D design - must have 0 lean or sweep and 1 chord_fac at ends
b.split{1}.sweep = spap2(1,3,[0 0.5 1],[0 0.02 0]);
b.split{1}.lean = spap2(1,3,[0 0.5 1],[0 -0.02 0]);
b.split{1}.fchord = spap2(1,3,[0 0.5 1],[1 0.7 1]);

% Half the thickness
b.thick_max.coefs = 0.5 * b.thick_max.coefs;

% % Copy the splitter for a second one
% b.split{2} = b.split{1}; b.split{2}.r_nondim = 0.7; 
% b.split{2}.chi_le.coefs = -b.split{2}.chi_le.coefs;
% b.split{2}.chi_te.coefs = -b.split{2}.chi_te.coefs;
% b.split{2}.lean.coefs = -b.split{2}.lean.coefs;
% 
% % Copy the splitter for a third
% b.split{3} = b.split{1}; b.split{3}.r_nondim = 0.5; 
% b.split{3}.chi_le.coefs(:) = 0; b.split{3}.chi_te.coefs(:) = 0;
% b.split{3}.lean.coefs(:) = 0;
% b.split = b.split([1 3 2]);

% Construct splitter geometry
b = bl_construct_splitter(b,0); 

% Count number of splitters
ns = length(b.split);

% Plot blade and splitter meshes
figure(); hold on; grid on; box on; axis equal; cols = lines(7);
mesh(b.xrrt(:,:,1),-b.xrrt(:,:,3),b.xrrt(:,:,2),...
    repmat(reshape(cols(1,:),[1 1 3]),[size(b.xrrt,1) size(b.xrrt,2) 1]))
for n = 1:ns
    mesh(b.split{n}.xrrt(:,:,1),-b.split{n}.xrrt(:,:,3),b.split{n}.xrrt(:,:,2),...
        repmat(reshape(cols(n+1,:),[1 1 3]),[size(b.split{n}.xrrt,1) size(b.split{n}.xrrt,2) 1]))
end


%% Prepare the blade surfaces for intersection calculation

% Interpolate the blade onto constant non-dimensional radius sections, helps high lean and sweep cases
sr = bl_passage_coords(b.xrrt,xr_hub,xr_cas); sr_new = linspace(0,1,size(sr,2));
for i = 1:size(b.xrrt,1)
    b.xrrt(i,:,:) = interp1(sr(i,:)',squeeze(b.xrrt(i,:,:)),sr_new','pchip');
end

% Calculate midspan non-dimensional axial spacing
x_mid = b.xrrt(:,round(size(b.xrrt,2)/2),1); 
[x_min,i_le] = min(x_mid); [x_max,i_te] = max(x_mid); 
x_mid = bl_split(x_mid,i_le,i_te);
sx_new = (x_mid - min(x_mid)) / (max(x_mid) - min(x_mid));

% Split into two sides at max and min axial coordinates and re-interpolate to same axial coordinates
xrrt_bl = cell(1,2); xr_le = zeros(length(sr_new),2); xr_te = zeros(length(sr_new),2); 
for j = 1:length(sr_new)
    
    % Split the section
    [~,i_le] = min(b.xrrt(:,j,1)); [~,i_te] = max(b.xrrt(:,j,1)); 
    [xrrt{1},xrrt{2}] = bl_split(squeeze(b.xrrt(:,j,:)),i_le,i_te);
    
    % Offset second side of the blade by a pitch
    xrrt{2}(:,3) = xrrt{2}(:,3) + xrrt{2}(:,2) * dt;
    
    % Record leading and trailing edge coordinates
    xr_le(j,:) = xrrt{1}(1,1:2); xr_te(j,:) = xrrt{1}(end,1:2);
    
    % Interpolate to the same non-dimensional axial spacing
    for n = 1:length(xrrt)
        sx = (xrrt{n}(:,1) - xrrt{n}(1,1)) / (xrrt{n}(end,1) - xrrt{n}(1,1));
        xrrt_bl{n}(:,j,:) = interp1(sx,xrrt{n},sx_new,'linear','extrap');
    end
end

% Trim and then extend the surfaces for intersection calculation
xrrt_surf = cell(1,2); sx_surf = [-1 ; sx_new(sx_new > 0.02 & sx_new < 0.98) ; 2];
for n = 1:length(xrrt)
    xrrt_surf{n} = reshape(interp1(sx_new,reshape(xrrt_bl{n},length(sx_new),[]),sx_surf,...
        'linear','extrap'),[length(sx_surf) length(sr_new) 3]);
    xrrt_surf{n}([1 end],:,2:3) = xrrt_surf{n}([2 end-1],:,2:3);
end

% Plot blade surface meshes
figure(); hold on; grid on; box on; axis equal; cols = lines(7); set(gca,'clippingstyle','rectangle'); 
for n = 1:numel(xrrt_bl)
    mesh(xrrt_bl{n}(:,:,1),-xrrt_bl{n}(:,:,3),xrrt_bl{n}(:,:,2),...
        repmat(reshape(cols(n,:),[1 1 3]),[size(xrrt_bl{n},1) size(xrrt_bl{n},2) 1]))
    mesh(xrrt_surf{n}(:,:,1),-xrrt_surf{n}(:,:,3),xrrt_surf{n}(:,:,2),zeros(size(xrrt_surf{n})))
end


%% Calculate mesh spacings

% Average distances in current design
j_mid = round(length(sr_new) / 2);
di = sum(sum(diff(xrrt_bl{1}(:,j_mid,:),1,1).^2,3).^0.5); 
dj = mean(xr_cas([1 end],2) - xr_hub([1 end],2));
dk = xrrt_bl{2}(1,j_mid,3) - xrrt_bl{1}(1,j_mid,3);

% Create spacings for current blade and splitter design 
[sk,rk] = geo_bunch(nk,0,1,ek/dk,ek/dk,0.31);
if ns ~= 0 
    dsp = b.split{1}.x_nondim(1) * di;
    [si_bl,ri] = geo_bunch(ni_bl,0,1,ei / (di - 2*dsp),ei / (di - 2 * dsp),0.5);
    [si_sp,ri_sp] = geo_bunch(ni_sp,0,1,0.5 * ei / dsp,0.5 * ei / dsp); ni = ni_bl + 2 * ni_sp - 2;
else
    [si_bl,ri] = geo_bunch(ni_bl,0,1,ei/di,ei/di,0.5); ni = ni_bl; ri_sp = []; dsp = 0.02 * di;
end

% Create radial spacings
r_nondim = 0; for n = 1:ns; r_nondim = [r_nondim b.split{n}.r_nondim]; end; r_nondim(1,end+1) = 1;
njs = (2:40) * 4 + 1;
for n = 1:ns+1
    rj_list = zeros(length(njs),1); sj_list = cell(length(njs),1); fj = r_nondim(n+1) - r_nondim(n);
    for m = 1:length(njs)
        [sj_list{m},rj_list(m)] = geo_bunch(njs(m),0,1,ej / (dj * fj),ej / (dj * fj),0.5);
    end
    m = find(rj_list < rj,1); nj(n) = njs(m); sj{n} = sj_list{m};
end 

% Spacings away from the blade
sk_in = linspace(0,1,nk); 
nj_in = sum(nj) - ns; sj_in = geo_bunch(nj_in,0,1,ej/dj,ej/dj,0.5);
si_out = geo_bunch(ni_out,0,1,2*ei/dx_duct);
si_in = flip(1 - geo_bunch(ni_in,0,1,2*ei/dx_duct),2);

% Endwall streamwise spacings with and without splitters
if ns == 0
    si_ew = si_bl;
else
    si_ew = geo_bunch(ni,0,1,ei/di,ei/di,0.5);
end

% Print the expansion ratios
disp(['Expansion ratios = ' num2str([ri rj rk ri_sp])]);


%% Calculate intersections and passage generate meshes

% Pitchwise coordinates
if ns == 0; st_new = [0 1]; else; st_new = linspace(0,1,size(b.split{1}.xrrt,2)-4); end;

% Trim splitter blades with intersection calculation and split into upper and lower sides
xrrt_sp = cell(ns,2); 
for n = 1:ns
    
    % Calculate passage coordinates for splitter
    b.split{n}.xrrt = b.split{n}.xrrt(:,3:end-2,:);
    [sr,sx,st] = bl_passage_coords(b.split{n}.xrrt,xr_hub,xr_cas,xr_le,xr_te,xrrt_surf{2},xrrt_surf{1});
    
    % Interpolate new pitchwise coordinates
    for i = 1:size(st,1)
        b.split{n}.xrrt(i,:,:) = interp1(st(i,:)',squeeze(b.split{n}.xrrt(i,:,:)),st_new','pchip');
    end
    
    % Loop over all sections
    for j = 1:length(st_new)
    
        % Split the section
        [~,i_le] = min(b.split{n}.xrrt(:,j,1)); [~,i_te] = max(b.split{n}.xrrt(:,j,1)); 
        [xrrt{1},xrrt{2}] = bl_split(squeeze(b.split{n}.xrrt(:,j,:)),i_le,i_te);

        % Interpolate to final non-dimensional axial spacing
        for m = 1:length(xrrt)
%             [~,sx] = bl_passage_coords(xrrt{m},xr_hub,xr_cas,xr_le,xr_te);
%             sx = (sx - sx(1)) / (sx(end) - sx(1));
% %             sx = (xrrt{m}(:,1) - xrrt{m}(1,1)) / (xrrt{m}(end,1) - xrrt{m}(1,1));
            sx = dist_2d(xrrt{m},1);
            xrrt_sp{n,m}(:,j,:) = interp1(sx,xrrt{m},si_bl,'pchip');      
        end
    end
    
    % Plot splitter surface meshes
    for m = 1:2
        mesh(xrrt_sp{n,m}(:,:,1),-xrrt_sp{n,m}(:,:,3),xrrt_sp{n,m}(:,:,2),...
            repmat(reshape(cols(m+2,:),[1 1 3]),[size(xrrt_sp{n,m},1) size(xrrt_sp{n,m},2) 1]))
    end
    
end

% Figure window for mesh segments
figure(); hold on; grid on; box on; axis equal; cols = lines(7); set(gca,'clippingstyle','rectangle'); 

% Break the blade into surfaces with the splitter
xrrt_bk = cell(ns+3,1); ks = [1 length(st_new)];
for n = 1:ns+1
    
    % Preallocate cell arrays
    xrrtj = cell(1,2); xrrtk = cell(1,2); sx = cell(2,2); sr = cell(2,2);
    
    % Determine endwall non-dimensional coordinates
    for m = 1:2

        % Get hub endwall curve or splitter for lower segment 
        if n == 1
            sr{m,1} = zeros(ni,1); sx{m,1} = si_ew';
        else
            [sr{m,1},sx{m,1}] = bl_passage_coords(xrrt_sp{n-1,2}(:,ks(m),:),xr_hub,xr_cas,xr_le,xr_te);
            sx{m,1} = [flip(1 - si_sp',1) * sx{m,1}(1) ; sx{m,1}(2:end-1) ; ...
                si_sp' * (1 - sx{m,1}(end)) + sx{m,1}(end)];
            sr{m,1} = [ones(ni_sp-1,1) * sr{m,1}(1) ; sr{m,1} ; ones(ni_sp-1,1) * sr{m,1}(end)];
        end
        
        % Get casing endwall curve or splitter for upper segment
        if n == ns+1
            sr{m,2} = ones(ni,1); sx{m,2} = si_ew';
        else
            [sr{m,2},sx{m,2}] = bl_passage_coords(xrrt_sp{n,1}(:,ks(m),:),xr_hub,xr_cas,xr_le,xr_te);
            sx{m,2} = [flip(1 - si_sp',1) * sx{m,2}(1) ; sx{m,2}(2:end-1) ; ...
                si_sp' * (1 - sx{m,2}(end)) + sx{m,2}(end)];
            sr{m,2} = [ones(ni_sp-1,1) * sr{m,2}(1) ; sr{m,2} ; ones(ni_sp-1,1) * sr{m,2}(end)];
        end
        
    end

    % Match extended coordinates at blade leading and trailing edges
    if ns ~= 0
    for o = 1:2
        sx_av = 0.5 * (sx{1,o} + sx{2,o}); sr_av = 0.5 * (sr{1,o} + sr{2,o});
        for m = 1:2
            sx{m,o}(1:ni_sp) = flip(si_sp',1) .* sx_av(1:ni_sp) + flip(1-si_sp',1) .* sx{m,o}(1:ni_sp);
            sr{m,o}(1:ni_sp) = flip(si_sp',1) .* sr_av(1:ni_sp) + flip(1-si_sp',1) .* sr{m,o}(1:ni_sp);
            sx{m,o}(ni-ni_sp+1:end) = si_sp' .* sx_av(ni-ni_sp+1:end) + ...
                (1-si_sp)' .* sx{m,o}(ni-ni_sp+1:end);
            sr{m,o}(ni-ni_sp+1:end) = si_sp' .* sr_av(ni-ni_sp+1:end) + ...
                (1-si_sp)' .* sr{m,o}(ni-ni_sp+1:end);
        end
    end
    end    
    
    % Interpolate blade surfaces
    for m = 1:2

        % Radial and axial spacing arrays
        sx_mat = sx{m,1} * (1 - sj{n}) + sx{m,2} * sj{n};
        sr_mat = sr{m,1} * (1 - sj{n}) + sr{m,2} * sj{n};
        
        % Correct splitter interface meshes
        if ns ~= 0
        
            % Bulge the axial spacing in vicinity of splitter leading and trailing edges
            dx = 0.2 * dsp / di;
            dx_2 = interp1([0 sx_mat(ni_sp,end) + [0 0.05] sx_mat(ni-ni_sp+1,end) + [-0.05 0] 1],...
                [0 -dx 0 0 dx 0],sx_mat,'pchip');
            dx_1 = interp1([0 sx_mat(ni_sp,1) + [0 0.05] sx_mat(ni-ni_sp+1,1) + [-0.05 0] 1],...
                [0 -dx 0 0 dx 0],sx_mat,'pchip');
            qr_2 = repmat(interp1([0 0.9 0.98 1],[0 0 1 0],sj{n},'pchip'),[ni 1]);
            qr_1 = repmat(interp1([0 0.02 0.1 1],[0 1 0 0],sj{n},'pchip'),[ni 1]);
            dx = dx_2 .* qr_2 + dx_1 .* qr_1;
            sx_mat = sx_mat + dx;

            % Smooth axial spacings
            qr = repmat(interp1([0 0.02 0.98 1],[0 1 1 0],sj{n},'pchip'),[ni 1]);
            sx_mat(3:end-2,:) = qr(3:end-2,:) .* (sx_mat(1:end-4,:) + sx_mat(2:end-3,:) + ...
                sx_mat(3:end-2,:) + sx_mat(4:end-1,:) + sx_mat(5:end,:)) / 5 + ...
                (1-qr(3:end-2,:)) .* sx_mat(3:end-2,:);
            sx_mat(2:end-1,:) = qr(2:end-1,:) .* (sx_mat(1:end-2,:) + sx_mat(2:end-1,:) + ...
                sx_mat(3:end,:)) / 3 + (1-qr(2:end-1,:)) .* sx_mat(2:end-1,:);
        
        end
        
        % Interpolate new blade coordinates
        for o = 1:3
            f = griddedInterpolant(repmat(sx_new,[1 length(sr_new)]),...
                repmat(sr_new,[length(sx_new) 1]),xrrt_bl{m}(:,:,o),'linear','linear');
            xrrtj{m}(:,:,o) = f(sx_mat,sr_mat);
        end
        
        % Plot blade mesh segment
        mesh(xrrtj{m}(:,:,1),-xrrtj{m}(:,:,3),xrrtj{m}(:,:,2),...
            repmat(reshape(cols(n,:),[1 1 3]),[ni nj(n) 1]))

    end
    
    % Get splitter or endwall surfaces
    i_le = 1:ni_sp; i_te = ni-ni_sp+1:ni;
    if n == 1 
        xrrt{1} = [xrrtj{1}(:,1,:) xrrtj{2}(:,1,:)];
    else
        xrrt{1} = xrrt_sp{n-1,2};
    end
    if n == ns+1
        xrrt{2} = [xrrtj{1}(:,end,:) xrrtj{2}(:,end,:)];
    else
        xrrt{2} = xrrt_sp{n,1};
    end
    
    % Interpolate splitter or endwall coordinates at desired pitchwise spacing
    for m = 1:2
        
        % Index for spanwise point on blade surface
        if m == 1; j = 1; else; j = nj(n); end;
        
        % Extend the splitter surfaces to the blade leading and trailing edges
        if size(xrrt{m},1) ~= ni
            
            % Calculate delta between splitter edges and up & downstream blade surface coordinates
            dxrrt = [[xrrtj{1}(i_le,j,:) xrrtj{2}(i_le,j,:)] - ...
                repmat(xrrt{m}(1,[1 end],:),[ni_sp 1 1]) ; zeros(ni_bl-2,2,3) ; ...
                [xrrtj{1}(i_te,j,:) xrrtj{2}(i_te,j,:)] - ...
                repmat(xrrt{m}(end,[1 end],:),[ni_sp 1 1])];
            
            % Linearly scale the delta across the passage
            dxrrt = repmat(dxrrt(:,1,:),[1 length(st_new) 1]) .* repmat(1 - st_new,[ni 1 3]) + ...
                repmat(dxrrt(:,end,:),[1 length(st_new) 1]) .* repmat(st_new,[ni 1 3]);
            
            % Apply the delta and calculate extended coordinates
            xrrt{m} = [repmat(xrrt{m}(1,:,:),[ni_sp-1 1 1]) ; xrrt{m} ; ...
                repmat(xrrt{m}(end,:,:),[ni_sp-1 1 1])] + dxrrt;
        end
        
        % Loop over all streamwise points and re-interpolate
        for i = 1:ni
            st = (xrrt{m}(i,:,3) - xrrt{m}(i,1,3)) ./ (xrrt{m}(i,end,3) - xrrt{m}(i,1,3));
            xrrtk{m}(i,:,:) = interp1(st',squeeze(xrrt{m}(i,:,:)),sk,'pchip');
        end
        
        % Correct leading edge axial coordinates for splitter sweep
        dx = 1.2 * (xrrtk{m}(1,:,1) - xrrtk{m}(1,1,1));
        dx(dx < 0) = 0; dx(dx > dsp) = dsp;
        if all(dx == 0)
            dx = interp1([0 0.02 0.98 1],[0 1 1 0] * 0.5 * dsp,sk,'pchip'); dx([1 end]) = 0;
        end
%         xrrtk{m}(i_le,:,1) = xrrtk{m}(i_le,:,1) - flip(si_sp,2)' * dx;
        xrrtk{m}(i_le,:,1) = xrrtk{m}(i_le,:,1) - flip(geo_bunch(ni_sp,0,1,0.1/ni_sp),2)' * dx;
        
        % Correct trailing edge coordinates
        dx = interp1([0 0.02 0.98 1],[0 1 1 0] * 0.5 * dsp,sk,'pchip'); dx([1 end]) = 0;
        xrrtk{m}(i_te,:,1) = xrrtk{m}(i_te,:,1) + geo_bunch(ni_sp,0,1,0.1/ni_sp)' * dx;
        
        % Correct endwall coordinates back
        if n == 1 && m == 1 
            xrrtk{m}(:,:,2) = interp1(xr_hub(:,1),xr_hub(:,2),xrrtk{m}(:,:,1),'pchip');
        elseif n == ns+1 && m == 2
            xrrtk{m}(:,:,2) = interp1(xr_cas(:,1),xr_cas(:,2),xrrtk{m}(:,:,1),'pchip');
        end
        
        % Bulge the mesh at both endwalls
%         dx = interp1([0 0.01 0.05 0.95 0.99 1],[0 0.5 1 1 0.5 0] * dsp,sk,'pchip');
%         dx = interp1([0 0.02 0.98 1],[0 1 1 0] * 0.5 * dsp,sk,'pchip'); dx([1 end]) = 0;
%         dx = linspace(0,dsp,round(nk/2)); dx = [dx flip(dx(1:end-1),2)]; 
%         xrrtk{m}(i_le,:,1) = xrrtk{m}(i_le,:,1) - flip(si_sp,2)' * dx;
%         xrrtk{m}(i_te,:,1) = xrrtk{m}(i_te,:,1) + si_sp' * dx;
        
        % Plot endwall mesh segment
        mesh(xrrtk{m}(:,:,1),-xrrtk{m}(:,:,3),xrrtk{m}(:,:,2),...
            repmat(reshape(cols(n,:),[1 1 3]),[ni nk 1]))
    end
    
    % Calculate block coordinates from blade surfaces
    xrrt_bk{n} = repmat(reshape(xrrtj{1},[ni nj(n) 1 3]),[1 1 nk 1]) ...
        .* repmat(reshape(1-sk,[1 1 nk 1]),[ni nj(n) 1 3]) + ...
        repmat(reshape(xrrtj{2},[ni nj(n) 1 3]),[1 1 nk 1]) .* ...
        repmat(reshape(sk,[1 1 nk 1]),[ni nj(n) 1 3]);
    
    % Calculate deltas from desired endwall and splitter surfaces
    dxrrt_1 = reshape(xrrtk{1},[ni 1 nk 3]) - xrrt_bk{n}(:,1,:,:);
    dxrrt_2 = reshape(xrrtk{2},[ni 1 nk 3]) - xrrt_bk{n}(:,end,:,:);
    dxrrt = repmat(1-sj{n},[ni 1 nk 3]) .* repmat(dxrrt_1,[1 nj(n) 1 1]) + ...
        repmat(sj{n},[ni 1 nk 3]) .* repmat(dxrrt_2,[1 nj(n) 1 1]);
    
    % Correct block coordinates to match endwalls and splitters
    xrrt_bk{n} = xrrt_bk{n} + dxrrt;
    
end


%% Extend passage meshes to create inlet and outlet blocks

% Collect passage inlet and outlet geometry
xrrt = cell(2); xrrt{1,2} = xrrt_bk{1}(1,:,:,:); xrrt{2,2} = xrrt_bk{1}(end,:,:,:);
for n = 1:ns
    xrrt{1,2} = [xrrt{1,2} xrrt_bk{n+1}(1,2:end,:,:)]; 
    xrrt{2,2} = [xrrt{2,2} xrrt_bk{n+1}(end,2:end,:,:)];
end

% Inlet and outlet surface meridional coordinates
xrrt{1,1} = zeros(size(xrrt{1,2})); xrrt{2,1} = zeros(size(xrrt{2,2})); 
xrrt{1,1}(1,:,:,1) = xr_hub(1,1); xrrt{2,1}(1,:,:,1) = xr_hub(end,1);
xrrt{1,1}(1,:,:,2) = repmat(sj_in * (xr_cas(1,2) - xr_hub(1,2)) + xr_hub(1,2),[1 1 nk 1]);
xrrt{2,1}(1,:,:,2) = repmat(sj_in * (xr_cas(end,2) - xr_hub(end,2)) + xr_hub(end,2),[1 1 nk 1]);

% Straight inlet and outlet coordinates
for n = 1:2
    xrrt{n,2}(:,:,:,3) = xrrt{n,2}(:,:,:,3) ./ xrrt{n,2}(:,:,:,2); t = squeeze(xrrt{n,2}(1,:,:,3));
    t_1 = sj_in * (t(end,1) - t(1,1)) + t(1,1);
    t_2 = sj_in * (t(end,end) - t(1,end)) + t(1,end);
    xrrt{n,1}(1,:,:,3) = repmat(reshape(1-sk_in,[1 1 nk]),[1 nj_in 1]) .* repmat(t_1,[1 1 nk]) + ...
        repmat(reshape(sk_in,[1 1 nk]),[1 nj_in 1]) .* repmat(t_2,[1 1 nk]);
end

% Generate inlet and outlet blocks
sis = {si_in,flip(1-si_out,2)}; nis = [ni_in ni_out];
for n = 1:2
    
    % Scale between domain edge and blade edge surfaces
    xrrt_bk{ns+n+1} = repmat(1-sis{n}',[1 nj_in nk 3]) .* repmat(xrrt{n,1},[nis(n) 1 1 1]) + ...
        repmat(sis{n}',[1 nj_in nk 3]) .* repmat(xrrt{n,2},[nis(n) 1 1 1]);
    
    % Add deltas to radius to correct hub and casing lines
    sr = (xrrt_bk{ns+n+1}(:,:,:,2) - repmat(xrrt_bk{ns+n+1}(:,1,:,2),[1 nj_in 1 1])) ./ ...
        repmat(xrrt_bk{ns+n+1}(:,end,:,2) - xrrt_bk{ns+n+1}(:,1,:,2),[1 nj_in 1 1]);
    r_hub = interp1(xr_hub(:,1),xr_hub(:,2),xrrt_bk{ns+n+1}(:,1,:,1),'pchip');
    r_cas = interp1(xr_cas(:,1),xr_cas(:,2),xrrt_bk{ns+n+1}(:,end,:,1),'pchip');
    dr_hub = r_hub - xrrt_bk{ns+n+1}(:,1,:,2);
    dr_cas = r_cas - xrrt_bk{ns+n+1}(:,end,:,2);
    xrrt_bk{ns+n+1}(:,:,:,2) = xrrt_bk{ns+n+1}(:,:,:,2) + sr .* repmat(dr_cas,[1 nj_in 1 1]) + ...
        (1-sr) .* repmat(dr_hub,[1 nj_in 1 1]);
    
    % Convert angular coordinate back to distance
    xrrt_bk{ns+n+1}(:,:,:,3) = xrrt_bk{ns+n+1}(:,:,:,3) .* xrrt_bk{ns+n+1}(:,:,:,2);
end

% Flip the outlet block
xrrt_bk{end} = flip(xrrt_bk{end},1);

% Plot block surfaces
figure(); hold on; grid on; box on; axis equal; cols = lines(7); set(gca,'clippingstyle','rectangle'); 
for n = 1:length(xrrt_bk)
    for k = [1 nk]
        mesh(xrrt_bk{n}(:,:,k,1),-xrrt_bk{n}(:,:,k,3),xrrt_bk{n}(:,:,k,2),...
            repmat(reshape(cols(n,:),[1 1 3]),[size(xrrt_bk{n},1) size(xrrt_bk{n},2) 1]))
    end
    for j = [1 size(xrrt_bk{n},2)]
        mesh(squeeze(xrrt_bk{n}(:,j,:,1)),-squeeze(xrrt_bk{n}(:,j,:,3)),squeeze(xrrt_bk{n}(:,j,:,2)),...
            repmat(reshape(cols(n,:),[1 1 3]),[size(xrrt_bk{n},1) size(xrrt_bk{n},3) 1]))
    end
end

% Check for negative axial spacings
for n = 1:length(xrrt_bk)
    q = diff(xrrt_bk{n}(:,:,:,1),1,1) < 0;
    if any(q(:)); error([num2str(sum(q(:))) ' Negative axial spacings']); end;
end


%% Patch together and preprocess for TURBOSTREAM

% Add mesh to the blocks in cell array
nn = length(xrrt_bk); g = cell(nn,1); varnames = {'x' 'r' 'rt'}; 
for n = 1:nn

    % Store coordinates
    for v = 1:length(varnames); g{n}.(varnames{v}) = xrrt_bk{n}(:,:,:,v); end; 

    % Block attributes
    a.bid = n-1; a.np = 2; a.ni = size(g{n}.x,1); a.nj = size(g{n}.x,2); a.nk = size(g{n}.x,3); 
    a.nc = 0; a.ncl = 0; a.procid = 0; a.threadid = 0;
    g{n}.attribute = a; g{n}.patch = cell(0,1); g{n}.bv.nblade = N; g{n}.bv.fblade = N;
end

% Approximate dynamic head
ro = (Po_in / (R * To_in));
dP = 0.5 * ro * (Vx_in / cosd(Alpha_in))^2;

% Calculate mass flow 
mass = pi * (xr_cas(1,2)^2 - xr_hub(1,2)^2) * ro * Vx_in;

% Default patch attributes
p.idir = 0; p.jdir = 1; p.kdir = 2; p.nface = 0; p.nt = 1; p.nxbid = 0; p.nxpid = 0;

% Apply exit boundary condition
p.ist = ni_out-1; p.ien = ni_out; p.jst = 0; p.jen = nj_in; p.kst = 0; p.ken = nk;
p.bid = 0; p.pid = 0; p.kind = 1;
g{end}.patch{1}.attribute = p;
g{end}.patch{1}.pout = Pref - dP; g{end}.patch{1}.throttle_type = 1; 
g{end}.patch{1}.throttle_target = mass; g{end}.patch{1}.throttle_k0 = 50; 
g{end}.patch{1}.throttle_k1 = 50; g{end}.patch{1}.throttle_k2 = 50; g{end}.patch{1}.ipout = 3;

% Turn off mass flow control
g{end}.patch{1}.throttle_type = 0;

% Inlet stagnation pressure profile
inlet.r_nondim = (g{end-1}.r(1,:,1) - g{end-1}.r(1,1,1)) / diff(g{end-1}.r(1,[1 end],1));
inlet.Po = interp1([0 d_hub 1-d_cas 1],[-dP 0 0 -dP] + Po_in,inlet.r_nondim,'pchip');

% Remaining inlet boundary condition application
inlet.Alpha = ones(1,nj_in) * Alpha_in; inlet.Beta = zeros(1,nj_in); inlet.To = ones(1,nj_in) * To_in;
p.ist = 0; p.ien = 1; p.pid = 0; p.kind = 0; g{end-1}.patch{1}.attribute = p;
g = ts_apply_inlet(g,inlet,1);

% Periodic patch across pitch upstream of the blade
g = ts_add_periodic_patch(nn-2,nn-2,[0 ni_in ; 0 nj_in ; 0 1],[0 ni_in ; 0 nj_in ; nk-1 nk],...
    [0 1 2],[0 1 2],g,5);

% Periodic patch across pitch downstream of the blade
g = ts_add_periodic_patch(nn-1,nn-1,[0 ni_out ; 0 nj_in ; 0 1],[0 ni_out ; 0 nj_in ; nk-1 nk],...
    [0 1 2],[0 1 2],g,5);

% Periodic patches from inlet and outlet to blade blocks
j1 = 0;
for n = 1:ns+1
    j2 = j1 + nj(n);
    g = ts_add_periodic_patch(n-1,nn-2,[0 1 ; 0 nj(n) ; 0 nk],[ni_in-1 ni_in ; j1 j2 ; 0 nk],...
        [0 1 2],[0 1 2],g,5);
    g = ts_add_periodic_patch(n-1,nn-1,[ni-1 ni ; 0 nj(n) ; 0 nk],[0 1 ; j1 j2 ; 0 nk],...
        [0 1 2],[0 1 2],g,5);
    j1 = j2 -1;
end

% Periodic patches between splitter blades
for n = 1:ns
    g = ts_add_periodic_patch(n-1,n,[0 ni_sp ; nj(n)-1 nj(n) ; 0 nk],[0 ni_sp ; 0 1 ; 0 nk],...
        [0 1 2],[0 1 2],g,5);
    g = ts_add_periodic_patch(n-1,n,[ni-ni_sp ni ; nj(n)-1 nj(n) ; 0 nk],[ni-ni_sp ni ; 0 1 ; 0 nk],...
        [0 1 2],[0 1 2],g,5);
end

% Create default TURBOSTREAM settings
g = ts_create_initial(g,nstep,nstep_poisson,[ro Vx_in 0 0 2e5 0.0002]);

% Reduce smoothing but run with mixng length turbulence model
g = ts_create_restart(g,nstep,1); g{1}.av.nstep_save_start = nstep - nstep_avg;
g{1}.av.poisson_nstep = nstep_poisson; g{1}.av.poisson_restart = 0;
g{1}.av.write_yplus = 1;

% Check patches are valid
ts_check_patches(g,2);

% Set the job name
name = ['test_blade_' num2str(ns) 'spl'];

% % Write mesh to PARAVIEW to check for negative volumes
% ts_export_paraview(g,[dr.ts name '_mesh.hdf5'],[],[],0);


%% Run the job in TBLOCK locally

% Change variables for running in TBLOCK
g{1}.av.nchange = 1000; g{1}.av.sfac = 0.01; g{1}.av.damp = 25; g{1}.av.facsec = 0.8; 
g{1}.av.superfac = 0; g{1}.av.ilos = 10; g{1}.av.rfvis = 0.5; g{1}.av.rfmix = 0.025;
g{1}.av.f_extrap = 0.8; g{1}.av.fangle = 0.8; g{1}.av.f_smth_mix = 0; g{1}.av.nstepup = 5;
g{1}.av.nlos = 5; g{1}.av.nsetvis = 3; g{1}.av.itwd = 0;
    
% Change variables to default in code rather than manual
%     g{1}.av.rfmix = 0.01; g{1}.av.f_extrap = 0.9; g{1}.av.fangle = 0.9; g{1}.av.f_smth_mix = 1;
    
% Set TBLOCK required inputs 
g{1}.av.timtot = 0;
for bid = ts_get_bids(g)
    g{bid+1}.bv.nsector = 1; g{bid+1}.bv.ile = 0; g{bid+1}.bv.ite = 0; 
end

% Write TBLOCK input files
directory = [dr.tb name filesep];
tb3d_write(g,directory)

% Save MATLAB file recording all input parameters
save([directory 'input.mat'],'g')

% Run the case
codedir = pwd; cd(directory);
% system('tblock.x < input.dat');                 % To output to screen
system('tblock.x < input.dat > log.txt');       % To output to file
cd(codedir);

% Plot the convergence history
tb_conv(directory);

% Read the output files
g = tb3d_export(directory,name);


%% Submit and monitor case to run in TURBOSTREAM on Wilkes

% % Reduce number of steps to debug
% % g{1}.av.poisson_nstep = 10; g{1}.av.nstep = 5;
% 
% % Submit the job to Wilkes
% job.name = name; job = ts_sbatch(job,g);
% 
% % Monitor the job and then transfer files back
% ts_smonitor(job);
% 
% % Plot the job convergence history
% ts_plot_conv(job);
% 
% % Read the output file
% g = ts_read_hdf5([dr.ts job.name '_avg.hdf5']);


%% Post-process the job

% Export for viewing in PARAVIEW
ts_export_paraview(g,[dr.pv name '.hdf5'],'Track',[],1);

% Take cuts at inlet and outlet
inlet = ts_structunstruct_cut(g,[b.x_in 0 ; b.x_in 10],[],2,1);
outlet = ts_structunstruct_cut(g,[b.x_out 0 ; b.x_out 10],[],2,1);

% Plot downstream wakes
figure(); ts_plot_wake(inlet,outlet,subplot(1,1,1),1,2,1,15);

% Plot radial whirl profile
figure(); hold on; grid on; box on; xlabel('Whirl Angle / ^\circ'); ylabel('Radius');
plot(ts_mass_average(outlet,'Alpha',2),ts_mass_average(outlet,'r_nondim',2));

% Calculate stagnation pressure loss coefficient
Yp = (ts_mass_average(inlet,'Po',3) - ts_mass_average(outlet,'Po',3)) / ...
    (ts_mass_average(inlet,'Po',3) - ts_mass_average(inlet,'P',3));
disp(['Loss Coefficient = ' num2str(Yp)]);




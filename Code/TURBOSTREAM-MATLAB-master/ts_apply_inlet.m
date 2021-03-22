function [g,bc] = ts_apply_inlet(g,inlet_file,plot_stuff)
% TS_APPLY_INLET  Apply an inlet boundary condition from file

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end


%% Get inlet plane data

% Determine file types
if isnumeric(inlet_file) == 1
    inlet.r_nondim = [0 1]; 
    inlet.Po = inlet_file(1) * ones(1,2); inlet.To = inlet_file(2) * ones(1,2); 
    inlet.Alpha = inlet_file(3) * ones(1,2); inlet.Beta = inlet_file(4) * ones(1,2); 
elseif isstruct(inlet_file) == 1
    inlet = inlet_file;
elseif isempty(regexp(inlet_file,'.mp.dat','match')) ~= 1
    inlet = read_mp(inlet_file);
elseif isempty(regexp(inlet_file,'Boundary_Conditions','match')) ~= 1
    inlet = read_q263(inlet_file);
elseif ischar(inlet_file) == 1
    inlet = read_bc(inlet_file);
end


%% Interpolate inlet boundary conditions onto to the mesh

% Find all inlet patches
for bid = ts_get_bids(g)
    for pid = ts_get_pids(bid,g)
        if g{bid+1}.patch{pid+1}.attribute.kind == 0
            
            % Size of patch
            nj = g{bid+1}.attribute.nj; nk = g{bid+1}.attribute.nk;
            
            % Calculate non-dimensional block coordinates
            r = squeeze(g{bid+1}.r(1,:,1));
            r_nondim = (r - min(r)) / (max(r) - min(r));
            
            % Record previous boundary condition
            if isfield(g{bid+1}.patch{pid+1},'yaw') ~= 0
                Alpha = g{bid+1}.patch{pid+1}.yaw(1:nj);
                Po = g{bid+1}.patch{pid+1}.pstag(1:nj);
                To = g{bid+1}.patch{pid+1}.tstag(1:nj);
                Beta = g{bid+1}.patch{pid+1}.pitch(1:nj);
            else
                Alpha = nan(size(r_nondim));
                Po = nan(size(r_nondim));
                To = nan(size(r_nondim));
                Beta = nan(size(r_nondim));
            end

            % Reinterpolate BCs onto the mesh
            bc.Po = interp1(inlet.r_nondim,inlet.Po,r_nondim,'linear','extrap');
            bc.To = interp1(inlet.r_nondim,inlet.To,r_nondim,'linear','extrap');
            bc.pitch = interp1(inlet.r_nondim,inlet.Beta,r_nondim,'linear','extrap');
            bc.alpha = interp1(inlet.r_nondim,inlet.Alpha,r_nondim,'linear','extrap');

            % Apply new inlet boundary conditions
            g{bid+1}.patch{pid+1}.pstag = repmat(bc.Po.',[nk 1]);
            g{bid+1}.patch{pid+1}.tstag = repmat(bc.To.',[nk 1]);
            g{bid+1}.patch{pid+1}.pitch = repmat(bc.pitch.',[nk 1]);
            g{bid+1}.patch{pid+1}.yaw = repmat(bc.alpha.',[nk 1]);
            
            % Set patch variables
            g{bid+1}.patch{pid+1}.sfinlet = 0.1;
            g{bid+1}.patch{pid+1}.rfin = 0.2;
            
            % Check patch attributes
            if g{bid+1}.patch{pid+1}.attribute.nt ~= 1; error('Check patch nt'); end;

        end
    end
end
            
% Plot the inlet boundary condition
if plot_stuff == 1
    
    % Stagnation temperature
    figure(); subplot(2,2,1); hold on; grid on; box on;
    ylabel('(r - r_{hub} / (r_{tip} - r_{hub})'); xlabel('To');
    plot(bc.To, r_nondim, 'r-')
    plot(To,  r_nondim, 'k-')

    % Stagnation pressure
    subplot(2,2,2); hold on; grid on; box on;
    ylabel('(r - r_{hub} / (r_{tip} - r_{hub})'); xlabel('Po');
    plot(bc.Po, r_nondim, 'r-')
    plot(Po,  r_nondim, 'k-')
    
    % Yaw angle
    subplot(2,2,3); hold on; grid on; box on;
    ylabel('(r - r_{hub} / (r_{tip} - r_{hub})'); xlabel('Yaw');
    plot(bc.alpha, r_nondim, 'r-')
    plot(Alpha,  r_nondim, 'k-')
    
    % Pitch angle
    subplot(2,2,4); hold on; grid on; box on;
    ylabel('(r - r_{hub} / (r_{tip} - r_{hub})'); xlabel('Pitch');
    plot(bc.pitch, r_nondim, 'r-')
    plot(Beta,  r_nondim, 'k-')

end


end


function inlet = read_q263(filename)
% Read from a Boundary_Conditions file

% Read in all text
f = fileread(filename);

% Find the inlet plane
i_in = regexp(f,'INLET\n\d+','end');
i_le = regexp(f,'LE\n','start');
f_inlet = f(i_in+2:i_le-1);

% Convert to a matrix
f_cell = regexp(f_inlet,'\s+','split');
f_cell = reshape(f_cell(2:end-1),7,[])';
f_num = str2double(f_cell);

% Get variables from the file
inlet.r = f_num(:,1); inlet.x = f_num(:,2); inlet.Vx = f_num(:,3); inlet.Vr = f_num(:,4); 
inlet.Vt = f_num(:,5); inlet.T = f_num(:,6); inlet.P = f_num(:,7); 

% Calculate primary flow variables
inlet.ro = inlet.P ./ (287 * inlet.T);
inlet.rovx = inlet.ro .* inlet.Vx;
inlet.rovr = inlet.ro .* inlet.Vr;
inlet.rorvt = inlet.ro .* inlet.r .* inlet.Vt;
inlet.V = (inlet.Vx.^2 + inlet.Vr.^2 + inlet.Vt.^2).^0.5;
inlet.roe = inlet.ro .* (inlet.T * 718 + 0.5 * inlet.V.^2);

% Secondary flow variables
inlet.rt = zeros(size(inlet.r));
inlet = ts_secondary(inlet);


end


function inlet = read_mp(filename)
% Read from a .mp.dat file

% Get MP Data from file
fid = fopen(filename,'r');
MP_variables = {'r', 'ro', 'Vx', 'Vt', 'Vr', 'P', 'spallart'};
line = fgetl(fid); A = regexp(line,'\s+','split');
nj = str2double(A{2}); rpm = str2double(A{4});
fclose(fid);

MP_Data = importdata(filename, ' ', 1);
MP_Data.data = MP_Data.data(~isnan(MP_Data.data));

for v = 1:length(MP_variables)
    inlet.(MP_variables{v}) = MP_Data.data((v-1)*nj + 1: v*nj);
end

% Defaults
inlet.omega = rpm * 2 * pi / 60;
inlet.ga = 1.400878;
inlet.R = 287.040283;

% Absolute quantities
inlet.V = (inlet.Vx.^2 + inlet.Vt.^2 + inlet.Vr.^2).^0.5;
inlet.M = inlet.V ./ (( inlet.ga*inlet.R* ( inlet.P ./ (inlet.ro*inlet.R) ) ).^0.5);
inlet.U = inlet.r * inlet.omega;

% Relative quantities
inlet.Vt = inlet.Vt;
inlet.Vt_rel = inlet.Vt - inlet.U;
inlet.V_rel = (inlet.Vx.^2 + inlet.Vt_rel.^2 + inlet.Vr.^2).^0.5;
inlet.M_rel = inlet.V_rel ./ (( inlet.ga*inlet.R* (inlet.P ./ (inlet.ro*inlet.R) ) ).^0.5);

% Pressures
inlet.Po = inlet.P .* ( 1 + (inlet.ga-1)/2 * (inlet.M.^2) ) .^ (inlet.ga/(inlet.ga-1));
inlet.Po_rel = inlet.P .* ( 1 + (inlet.ga-1)/2 * (inlet.M_rel.^2) ) .^ (inlet.ga/(inlet.ga-1));

% Temperatures
inlet.T = inlet.P ./ (inlet.R*inlet.ro);
inlet.To_rel = inlet.T .* (1 + (inlet.ga-1)/2 * inlet.M_rel.^2);
inlet.To = inlet.T .* (1 + (inlet.ga-1)/2 * inlet.M.^2);

% Blade Speed
inlet.M_blade = inlet.U ./ (inlet.ga*inlet.R*inlet.T).^0.5;

% Whirl angles
inlet.Alpha = asin(inlet.Vt ./ inlet.V) * 360 / (2*pi());
inlet.Alpha_rel = asin(inlet.Vt_rel ./ inlet.V_rel) * 360 / (2*pi());

% Pitch angle
inlet.Beta = atand(inlet.Vr./inlet.Vx);
inlet.Beta(1) = 0;
inlet.Beta(end) = 0;

% Non dimensional r-coordinate
inlet.r_nondim = (inlet.r - min(inlet.r)) / (max(inlet.r) - min(inlet.r));

% Meridional velocity
inlet.Vm = inlet.V .* cosd(inlet.Alpha);


end


function inlet = read_bc(filename)
% Read from a specific inlet matrix

% Read matrix
M = dlmread(filename,'',1,0);
inlet.r = M(:,1);
inlet.To = M(:,2);
inlet.Po = M(:,3);
inlet.Alpha = M(:,4); 
inlet.Beta = M(:,5); 
inlet.SA = M(:,6);

% Calculate non-dimensional radius
inlet.r_nondim = (inlet.r - min(inlet.r)) / (max(inlet.r) - min(inlet.r));


end



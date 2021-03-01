function [t_mean,loss_var,outlet] = ts_plot_wake(inlet,outlet,h_axes,loss_type,...
    n_passages,coord_type,l,t_interp)
% Function to plot loss in cut of the mesh

%% Prepare cutplanes and defaults

% Choose loss type
% 0 - Stagnation pressure loss coefficient
% 1 - Entropy loss coefficient
if exist('loss_type','var') == 0
    loss_type = 1; 
end

% Choose number of passage repeats
if exist('n_passages','var') == 0
    n_passages = 1;
end

% Choose coordinate type
% 0 - r vs rt
% 1 - z vs y
% 2 - r_nondim vs rt_nondim
if exist('coord_type','var') == 0
    coord_type = 0;
end

% Choose contour levels
if exist('l','var') == 0
    l = 20;
end

% Choose whether to interpolate onto constant theta edges
if exist('t_interp','var') == 0
    t_interp = 0;
end

% Calculate secondary variables
inlet = ts_secondary(inlet,1); outlet = ts_secondary(outlet,1);


%% Reinterpolate outlet cut onto constant theta edged mesh
    
nk = size(outlet.r,2); nj = size(outlet.r,1);
j_mid = round(nj/2);

% Reinterpolate data onto a mesh with constant theta edges
if t_interp == 1
    
    % Prepare mesh size
    t_min = outlet.t(j_mid,1); t_max = outlet.t(j_mid,end);
    dt = t_max - t_min;

    % Define new ultra-fine tangential mesh
    t_new = linspace(t_min,t_max,nk*5);

    % Shift current coordinates +/- a pitch
    t_interp = [outlet.t(:,1:end-1) - dt outlet.t outlet.t(:,2:end) + dt];
    
    % Loop over all primary variables
    varnames = {'x' 'r' 'ro' 'rovx' 'rovr' 'rorvt' 'roe'};
    for v = 1:length(varnames)
        
        % Copy current variable to +/- a pitch
        var_interp = [outlet.(varnames{v})(:,1:end-1) outlet.(varnames{v}) ...
            outlet.(varnames{v})(:,2:end)];
        
        % Loop over all j indices and interpolate flow and grid coordinates
        var_new = zeros(nj,length(t_new));
        for j = 1:size(outlet.x,1)
            var_new(j,:) = interp1(t_interp(j,:),var_interp(j,:),...
                t_new,'pchip');
        end
        
        % Assign new variable
        outlet.(varnames{v}) = var_new;
        
    end
    
    % Assign new tangential coordinates
    outlet.rt = repmat(t_new,[nj 1]) .* outlet.r;
    nk = size(outlet.rt,2);
    
    % Recalculate secondary properties
    outlet = ts_secondary(outlet,1);
end

% Calculate loss parameter
if loss_type == 1
    % Change in lost work between outlet and mass averaged inlet value
    s1 = ts_mass_average(inlet, 's', 3);
    T2 = ts_mass_average(outlet, 'T', 3);
    delta_Ts = T2*(outlet.s - s1);

    % Calculate entropy generation coefficent
    V1 = ts_mass_average(inlet, 'V', 3);
    loss_var = 2*delta_Ts / (V1^2);
elseif loss_type == 0
    % Find stagnation pressure loss coefficient in absolute or relative frame
    if inlet.bv.rpm == 0
        Po1 = ts_mass_average(inlet,'Po',3);
        P1 = ts_mass_average(inlet,'P',3);
        loss_var = (Po1 - outlet.Po) / (Po1 - P1);
    else
        Po1 = ts_mass_average(inlet,'Po_rel',3);
        P1 = ts_mass_average(inlet,'P',3);
        loss_var = (Po1 - outlet.Po_rel) / (Po1 - P1);
    end
elseif loss_type == 2
    % Stagnation pressure loss coefficient referenced to mass fraction
    % upstream
    [Po_in, mass_in] = ts_mass_average(inlet,'Po',2);
    [~, mass_out] = ts_mass_average(outlet,'Po',2);
    dP = ts_mass_average(inlet,'Po',3) - ts_mass_average(inlet,'P',3);
    
    % Calculate mass fractions
    mass_in = [0 cumsum(mass_in)]; mass_out = [0 cumsum(mass_out)];
    mass_in = mass_in / max(mass_in); mass_out = mass_out / max(mass_out);
    
    % Interpolate total pressure upstream at exit mass fractions
    Po_in = repmat(reshape(interp1(mass_in,Po_in,mass_out,'pchip'),[nj 1]),[1 nk]);
    loss_var = (Po_in - outlet.Po) / dP;
end

% Nondimensional plotting coordinates
r_nondim = (outlet.r - min(min(outlet.r))) /...
    (max(max(outlet.r)) - min(min(outlet.r)));
rt_nondim = (outlet.rt - repmat(min(outlet.rt,[],2),1,size(outlet.rt,2))) ./...
   (repmat(max(outlet.rt,[],2),1,size(outlet.rt,2)) - repmat(min(outlet.rt,[],2),1,size(outlet.rt,2)));
t = outlet.rt./outlet.r;

% Remove any NaNs created
rt_nondim(isnan(rt_nondim)) = 0;
loss_var(isnan(loss_var)) = 0;

% Switch to current axes
axes(h_axes)
hold on

% Check cut is structured
if isfield(outlet,'faces') == 0

    % Plot multiple passages
    if n_passages ~= 1
        nblade = double(outlet.bv.nblade);
        t_plot = [];
        for p = 0:n_passages-1
            t_plot = [t_plot (t(:,2:end) + p*2*pi/nblade)];
        end
        t = t_plot;
        r_nondim = repmat(r_nondim(:,2:end),1,n_passages);
        r = repmat(outlet.r(:,2:end),1,n_passages);
        loss_var = repmat(loss_var(:,2:end),1,n_passages);
    else
        r = outlet.r;
    end

    [~, k] = sort(t(1,:));
    t = t(:,k);
    r_nondim = r_nondim(:,k);
    r = r(:,k);
    loss_var = loss_var(:,k);
else
    r = outlet.r;
end

% Shift theta direction so the centre is at 0 degrees
[~, mid_out] = min(r_nondim(:,1) - 0.5);
t_mean = mean(t(mid_out,:));
t = t - t_mean;
rt = t .* r;

if coord_type == 1    
    [r,rt] = pol2cart(t,r);
end
if coord_type == 2
    rt = rt_nondim; r = r_nondim;
end

if isfield(outlet,'faces') == 0
    contourf(rt, r, loss_var,l);
else
    patch('Vertices',[rt r -ones(size(rt))],'Faces',outlet.faces,...
        'FaceVertexCData',loss_var,'FaceColor','interp','edgecolor','interp');
end

% Set axes limits
axis equal
if coord_type == 0 || coord_type == 1
    axis equal
    axis([min(min(rt)) max(max(rt)) min(min(r)) max(max(r))])
elseif coord_type == 2
    axis equal
    axis([0 1 0 1])
end
colormap(flipud(parula));

caxis([floor(min(min(loss_var*10)))/10 ceil(max(max(loss_var*10)))/10])
set(gca,'Visible','off')

end
function [] = exp_traverse_monitor(filename)
% Plot output from 5 hole probe traverse as it is running

% Load environment variables
load('current_traverse.mat');

if strcmp(trav_type,'FH') == 1
    load(calib_name)
end

% Set scrolling time delay
t_delay = 200;

% Choose current operating point
des = des_current;
% N_plot.Po = N.Po_Rint; N_plot.P = N.P_Rint; N_plot.T = N.To_Rint;


%% Open figure window
h.window = figure(); set(h.window,'Position',[1 29 1920 985]);

% Five hole probe monitoring
if strcmp(trav_type,'FH') == 1
    
    % Window layout
    if isfield(des,'phi') == 1
        pos.Loss = [0.08 0.33 0.40 0.64];
        pos.phi = [0.08 0.06 0.40 0.24];
    else
        pos.Loss = [0.08 0.45 0.40 0.53];
        pos.mrtp = [0.08 0.06 0.40 0.15];
        pos.nrt = [0.08 0.27 0.40 0.15];
    end
    pos.Alpha = [0.52 0.66 0.40 0.32];
    pos.Error = [0.52 0.53 0.40 0.11];
    pos.P = [0.52 0.06 0.40 0.20];
    pos.T = [0.52 0.32 0.40 0.15];
    
    % Loss traverse axes
    h.Loss = axes('Position',pos.Loss,'XTickLabel','','YTickLabel',''); hold on;
    title('Loss'); box on; colormap(h.Loss,flipud(bone));
    colorbar('WestOutside'); %caxis([-0.05 0.45]);

    % Flow angle traverse axes
    h.Alpha = axes('Position',pos.Alpha,'XTickLabel','','YTickLabel',''); hold on;
    title('Yaw / deg');  box on; colormap(h.Alpha,'jet');
    colorbar('EastOutside'); %caxis([min(min(Iota))-15 max(max(Iota))+15]);
    
    % Five hole error axes
    h.Error = axes('Position',pos.Error); hold on; grid on; box on; 
    xlabel('Time / s'); ylabel('FHP Error'); axis([-t_delay 0 -0.2 5]);
    plot([-t_delay 0],[0.5 0.5],'r-');
end

% Total temperature monitoring
if strcmp(trav_type,'To') == 1
    
    % Loss traverse axes
    h.Temp = axes('Position',[0.08 0.07 0.40 0.85],'XTickLabel','','YTickLabel',''); 
    hold on;
    title('Stagnation Temperature'); box on; colormap(h.Temp,'jet');
    colorbar('WestOutside');
    
end

% Hot wire monitoring
if strcmp(trav_type,'HW') == 1
    
    % Window layout
    if isfield(des,'phi') == 1
        pos.phi = [0.08 0.06 0.40 0.41];
    else
        pos.mrtp = [0.08 0.06 0.40 0.17];
        pos.nrt = [0.08 0.30 0.40 0.17];
    end
    pos.V1 = [0.08 0.49 0.27 0.49];
    pos.V2 = [0.36 0.49 0.27 0.49];
    pos.V3 = [0.64 0.49 0.27 0.49];
    pos.P = [0.52 0.06 0.40 0.20];
    pos.T = [0.52 0.32 0.40 0.15];
    
    % Voltage axes
    varnames = fieldnames(pos);
    for v = 1:length(varnames)
        if strcmp(varnames{v}(1),'V') == 1
            h.(varnames{v}) = axes('Position',pos.(varnames{v}),'XTickLabel','','YTickLabel',''); 
            hold on;
            title('HW Voltage / V'); box on; colormap(h.(varnames{v}),'jet');
            colorbar('EastOutside'); 
        end
    end
end

% Rig operating point windows
if isfield(des,'phi') == 1
    % Variation of flow coefficient with time
    h.phi = axes('Position',pos.phi); hold on; grid on; box on; 
    xlabel('Time / s'); ylabel('Flow Coef.');
    axis([-t_delay 0 des.phi*(1-err.phi*3) des.phi*(1+err.phi*3)]);
else
    h.mrtp = axes('Position',pos.mrtp); hold on; grid on; box on; 
    xlabel('Time / s'); ylabel('Flow Func.');  
    axis([-t_delay 0 des.mrtp*(1-err.mrtp*3) des.mrtp*(1+err.mrtp*3)]);
    h.nrt = axes('Position',pos.nrt); hold on; grid on; box on; 
    xlabel('Time / s'); ylabel('Non. Dim. Speed');  
    axis([-t_delay 0 des.nrt*(1-err.nrt*3) des.nrt*(1+err.nrt*3)]);
end

% Raw data windows
h.P = axes('Position',pos.P); hold on; grid on; box on; 
xlabel('Time / s'); ylabel('Raw Pressure'); axis([-t_delay 0 -7000 7000]);
h.T = axes('Position',pos.T); hold on; grid on; box on; 
xlabel('Time / s'); ylabel('Raw Temperature'); axis([-t_delay 0 290 315]);

% Legends for raw data
P_legend = cell(P_chan,1); T_legend = cell(tc08.nchan,1);
for n = 1:P_chan; P_legend{n} = num2str(n); end
for n = 1:T_chan; T_legend{n} = num2str(n); end
T_cols = distinguishable_colors(T_chan);
P_cols = distinguishable_colors(P_chan);


%% Plot traverse grid

% Wait until traverse file is first written
while exist([filename '.mat'],'file') == 0
    pause(1);
end

load([filename '.mat']);
N = e.N;

% Edit coordinates to allow plots if line traverses are being used
if size(e.r,1) == 1
    e.r(2,:,:) = e.r + 3e-3;
elseif size(e.r,2) == 1
    e.t(:,2,:) = e.t + 1;
end

% Plot traverse grid in cartesians
[x,y] = pol2cart(2*pi*(e.t(:,:,1)+90)/360,e.r(:,:,1));
[x_hub,y_hub] = pol2cart(2*pi*(linspace(min(e.t(1,:)),max(e.t(1,:)),100)+90)/360,ones(1,100)*r_hub);
[x_cas,y_cas] = pol2cart(2*pi*(linspace(min(e.t(end,:)),max(e.t(end,:)),100)+90)/360,ones(1,100)*r_cas);
if strcmp(trav_type,'FH') == 1
    axes(h.Loss);
    plot3(x,y,ones(nj,nk),'k-');
    plot3(x',y',ones(nj,nk)','k-');
    plot(x_hub,y_hub,'k--'); plot(x_cas,y_cas,'k--');
    v = axis; axis(v); daspect([max(daspect)*[1 1] 1]); axis equal; caxis('auto');
    axes(h.Alpha);
    plot3(x,y,ones(nj,nk),'k-');
    plot3(x',y',ones(nj,nk)','k-');
    plot(x_hub,y_hub,'k--'); plot(x_cas,y_cas,'k--');
    v = axis; axis(v); daspect([max(daspect)*[1 1] 1]); axis equal; caxis('auto');
end
if strcmp(trav_type,'To') == 1
    axes(h.Temp); 
    plot3(x,y,ones(nj,nk),'k-');
    plot3(x',y',ones(nj,nk)','k-');
    plot(x_hub,y_hub,'k--'); plot(x_cas,y_cas,'k--');
    v = axis; axis(v); daspect([max(daspect)*[1 1] 1]); axis equal; caxis('auto');
end
if strcmp(trav_type,'HW') == 1
    varnames = fieldnames(h);
    for n = 1:length(varnames)
        if strcmp(varnames{n}(1),'V') == 1
            axes(h.(varnames{n})); 
            plot3(x,y,ones(nj,nk),'k-');
            plot3(x',y',ones(nj,nk)','k-');
            plot(x_hub,y_hub,'k--'); plot(x_cas,y_cas,'k--');
            v = axis; axis(v); daspect([max(daspect)*[1 1] 1]); axis equal; caxis('auto');
        end
    end
end

% Calculate offset grid for plotting
x_off = 0.25 * (x(1:end-1,1:end-1) + x(1:end-1,2:end) + x(2:end,1:end-1) + x(2:end,2:end));
y_off = 0.25 * (y(1:end-1,1:end-1) + y(1:end-1,2:end) + y(2:end,1:end-1) + y(2:end,2:end));
x_off = [x(1,1) 0.5 * (x(1,1:end-1) + x(1,2:end)) x(1,end) ; ...
    0.5 * (x(1:end-1,1) + x(2:end,1)) x_off 0.5 * (x(1:end-1,end) + x(2:end,end)) ;...
    x(end,1) 0.5 * (x(end,1:end-1) + x(end,2:end)) x(end,end)];
y_off = [y(1,1) 0.5 * (y(1,1:end-1) + y(1,2:end)) y(1,end) ; ...
    0.5 * (y(1:end-1,1) + y(2:end,1)) y_off 0.5 * (y(1:end-1,end) + y(2:end,end)) ;...
    y(end,1) 0.5 * (y(end,1:end-1) + y(end,2:end)) y(end,end)];


%% Configure all axes for quicker plotting

% Set automatic plot properties to manual
% pn = {'ALimMode','CameraPositionMode','CameraTargetMode',...
%     'CameraUpVectorMode','CLimMode','TickDirMode','XLimMode',...
%     'YLimMode','ZLimMode','XTickMode','YTickMode',...
%     'ZTickMode','XTickLabelMode','YTickLabelMode','ZTickLabelMode'};
% pv = cell(size(pn)); for n = 1:length(pn); pv(n) = {'manual'}; end;
% 
% h_axes = findobj(h.window,'type','axes');
% for n = 1:length(h_axes)
%     set(h_axes(n),pn,pv);
% end


%% Initialise time series plots

% Plot desired operating points and errors
varnames = fieldnames(des);
for v = 1:length(varnames)
    axes(h.(varnames{v})); 
    plot([-t_delay 0],des.(varnames{v})*ones(1,2),'k-');
    plot([-t_delay 0],des.(varnames{v})*ones(1,2) * (1 + err.(varnames{v})),'r-');
    plot([-t_delay 0],des.(varnames{v})*ones(1,2) * (1 - err.(varnames{v})),'r-');
end

% Initialise rig operating point
nt = round(t_delay*5);
time_stream = -200*ones(nt,1);
if isfield(des,'phi') == 1
    axes(h.phi);
    phi = nan(nt,1);
    l.phi = line(time_stream,phi);
else
    axes(h.mrtp);
    mrtp = nan(nt,1);
    l.mrtp = line(time_stream,mrtp);
    axes(h.nrt);
    nrt = nan(nt,1);
    l.nrt = line(time_stream,nrt);
end

% Initialise streaming pressure data
time_trav = -200*ones(length(e.J),1);
if strcmp(trav_type,'FH') == 1
    A = nan(length(e.J),1); axes(h.Error); l.Error = line(time_trav,A);
    P = nan(length(e.J),P_chan);
    axes(h.P); for n = 1:P_chan; l.P(n) = line(time_trav,P(:,n),'color',P_cols(n,:)); end;
elseif strcmp(trav_type,'HW') == 1
    P = nan(nt,P_chan);
    axes(h.P); for n = 1:P_chan; l.P(n) = line(time_stream,P(:,n),'color',P_cols(n,:)); end;
end
legend(P_legend,'Location','West','FontSize',6);

% Streaming temperatures
T = nan(nt,T_chan);
axes(h.T); for n = 1:T_chan; l.T(n) = line(time_stream,T(:,n),'color',T_cols(n,:)); end;
legend(T_legend,'Location','West','FontSize',6);


%% Plot traverse results as they are recorded

% Initialise plotted array
u = zeros(size(e.J)); u_new = zeros(size(e.J));
Yp = zeros(size(e.J)); Alpha = zeros(size(e.J));
To = zeros(size(e.J)); 
V_plot = cell(3); for v = 1:3; V_plot{v} = zeros(size(e.J)); end; V_max = -inf; V_min = inf;

% Memory cap
mem = memory;
if mem.MemUsedMATLAB > 2e9
    exit;
end

while sum(u == 0) ~= 0
    % Load in latest traverse data
    file_loaded = 0;
    while file_loaded == 0
        try
            load([filename '.mat'])
            file_loaded = 1;
        catch
            file_loaded = 0;
        end
    end
    
    % Memory cap
    mem = memory;
    if mem.MemUsedMATLAB > 2e9
        exit;
    end
    
    % Current time
    time_current = clock;
    if isempty(e.rig.time) == 0
        time_current = e.rig.time(end,:);
        % Elapsed rig time
        nr = min(size(e.rig.time,1),nt);
        time_temp = e.rig.time(end-nr+1:end,:);
        for m = 1:size(time_temp,1)
            time_stream(m) = etime(time_temp(m,:),time_current);
        end
        % Get rig operating point results
        if isfield(des,'phi') == 1
            phi(1:nr,:) = e.rig.phi(end-nr+1:end,:);
        else
            mrtp(1:nr,:) = e.rig.mrtp(end-nr+1:end,:);
            nrt(1:nr,:) = e.rig.nrt(end-nr+1:end,:);
        end
        if strcmp(trav_type,'FH') == 1
            T(1:nr,:) = e.rig.T(end-nr+1:end,:);
        end
        if strcmp(trav_type,'HW') == 1
            T(1:nr,:) = e.rig.T(end-nr+1:end,:);
            P(1:nr,:) = e.rig.P(end-nr+1:end,:);
        end
    end
    
    % Memory cap
    mem = memory;
    if mem.MemUsedMATLAB > 2e9
        exit;
    end

    % Loop over every point in grid and check if already traversed and
    % if plotted
    for n = 1:length(u)
        j = J(n); k = K(n);
        if u(n) == 0 && isnan(e.time(j,k,1)) == 0
            
            % Elapsed time
            time_trav(n) = etime(squeeze(e.time(j,k,:)).',time_current);

            % Get voltages and update colour scales
            if strcmp(trav_type,'HW') == 1
                if sum(isnan(e.V_av(j,k,:))) == 0
                    V_plot{1}(n) = e.V_av(j,k,1);
                    i_mid = round(ni/2);
                    V_plot{2}(n) = e.V_av(j,k,i_mid);
                    V_plot{3}(n) = e.V_av(j,k,end);
                    V_max = max([V_max max(e.V_av(j,k,:))]);
                    V_min = min([V_min min(e.V_av(j,k,:))]) - 1e-6;
                else
                    continue
                end
            end
            
            % Get loss, angles and pressures
            if strcmp(trav_type,'FH') == 1
                
                % Get pressures for loss coefficient
                Po = e.P(j,k,N_plot.Po); Ps = e.P(j,k,N_plot.P);

                % Calculate five hole probe results
                d = exp_fhp_calculate(c,e.P(j,k,N.P_FHP),probe.(probe_name),e.Iota(j,k),0);
                
                % Get five hole probe error
                d_temp = exp_fhp_squares(c,e.P(j,k,N.P_FHP),probe.(probe_name),...
                    e.Iota(j,k),0,'hodson',0);
                A(n) = mean(mean(d_temp.A));
                
                % Calculate loss coefficient
                Yp(n) = (Po - d.Po) / (Po - Ps);

                % Yaw angle
                Alpha(n) = d.Alpha;

                % Record raw pressures            
                P(n,:) = e.P(j,k,:);
            end
            
            % Pull out total temperature
            if strcmp(trav_type,'To') == 1
                To(n) = e.T(j,k,N.ToP);
                T(n,:) = e.T(j,k,:);
            end

            % List to plot patches later
            u_new(n) = 1; u(n) = 1;
            
        elseif isnan(e.time(j,k,1)) == 0
            % Update elapsed time for already plotted results
            time_trav(n) = etime(squeeze(e.time(j,k,:)).',time_current);
        end
    end
    
    % Memory cap
    mem = memory;
    if mem.MemUsedMATLAB > 2e9
        exit;
    end

    % Plot new patches of loss and flow angle
    if sum(u_new == 1) ~= 0
        n_new = reshape(find(u_new == 1),1,[]);
        if strcmp(trav_type,'FH') == 1
            set(h.window,'CurrentAxes',h.Loss);
            for n = n_new
                patch([x_off(J(n),K(n)) ; x_off(J(n)+1,K(n)) ; x_off(J(n),K(n)+1)],...
                    [y_off(J(n),K(n)) ; y_off(J(n)+1,K(n)) ; y_off(J(n),K(n)+1)],...
                    zeros(3,1),ones(3,1)*Yp(n),'EdgeColor','none','FaceColor','flat')
                patch([x_off(J(n)+1,K(n)+1) ; x_off(J(n)+1,K(n)) ; x_off(J(n),K(n)+1)],...
                    [y_off(J(n)+1,K(n)+1) ; y_off(J(n)+1,K(n)) ; y_off(J(n),K(n)+1)],...
                    zeros(3,1),ones(3,1)*Yp(n),'EdgeColor','none','FaceColor','flat')        
            end

            set(h.window,'CurrentAxes',h.Alpha);
            for n = n_new
                patch([x_off(J(n),K(n)) ; x_off(J(n)+1,K(n)) ; x_off(J(n),K(n)+1)],...
                    [y_off(J(n),K(n)) ; y_off(J(n)+1,K(n)) ; y_off(J(n),K(n)+1)],...
                    zeros(3,1),ones(3,1)*Alpha(n),'EdgeColor','none','FaceColor','flat')
                patch([x_off(J(n)+1,K(n)+1) ; x_off(J(n)+1,K(n)) ; x_off(J(n),K(n)+1)],...
                    [y_off(J(n)+1,K(n)+1) ; y_off(J(n)+1,K(n)) ; y_off(J(n),K(n)+1)],...
                    zeros(3,1),ones(3,1)*Alpha(n),'EdgeColor','none','FaceColor','flat')      
            end
        end
        if strcmp(trav_type,'To') == 1
            set(h.window,'CurrentAxes',h.Temp);
            for n = n_new
                patch([x_off(J(n),K(n)) ; x_off(J(n)+1,K(n)) ; x_off(J(n),K(n)+1)],...
                    [y_off(J(n),K(n)) ; y_off(J(n)+1,K(n)) ; y_off(J(n),K(n)+1)],...
                    zeros(3,1),ones(3,1)*To(n),'EdgeColor','none','FaceColor','flat')
                patch([x_off(J(n)+1,K(n)+1) ; x_off(J(n)+1,K(n)) ; x_off(J(n),K(n)+1)],...
                    [y_off(J(n)+1,K(n)+1) ; y_off(J(n)+1,K(n)) ; y_off(J(n),K(n)+1)],...
                    zeros(3,1),ones(3,1)*To(n),'EdgeColor','none','FaceColor','flat')      
            end
        end
        if strcmp(trav_type,'HW') == 1
            for n = n_new
                for v = 1:3
                    set(h.window,'CurrentAxes',h.(['V' num2str(v)]));
                    patch([x_off(J(n),K(n)) ; x_off(J(n)+1,K(n)) ; x_off(J(n),K(n)+1)],...
                        [y_off(J(n),K(n)) ; y_off(J(n)+1,K(n)) ; y_off(J(n),K(n)+1)],...
                        zeros(3,1),ones(3,1)*V_plot{v}(n),'EdgeColor','none','FaceColor','flat')
                    patch([x_off(J(n)+1,K(n)+1) ; x_off(J(n)+1,K(n)) ; x_off(J(n),K(n)+1)],...
                        [y_off(J(n)+1,K(n)+1) ; y_off(J(n)+1,K(n)) ; y_off(J(n),K(n)+1)],...
                        zeros(3,1),ones(3,1)*V_plot{v}(n),'EdgeColor','none','FaceColor','flat')  
                    caxis([V_min V_max]);
                end
            end
        end

        % Refresh new to plot array
        u_new = zeros(size(u_new));

        % Update time series traces
        if isfield(des,'phi') == 1
            set(l.phi,'XData',time_stream,'YData',phi);
        else
            set(l.mrtp,'XData',time_stream,'YData',mrtp);
            set(l.nrt,'XData',time_stream,'YData',nrt);
        end
        if strcmp(trav_type,'FH') == 1
            for n = 1:length(l.P)
                set(l.P(n),'XData',time_trav,'YData',P(:,n));
            end
            set(l.Error,'XData',time_trav,'Ydata',A);            
        elseif strcmp(trav_type,'HW') == 1
            for n = 1:length(l.P)
                set(l.P(n),'XData',time_stream,'YData',P(:,n));
            end
        end
        for n = 1:length(l.T)
            set(l.T(n),'XData',time_stream,'YData',T(:,n));
        end
    end
    % Build in delay
    pause(1);
    
    % Memory cap
    mem = memory;
    if mem.MemUsedMATLAB > 2e9
        exit;
    end

end


%% Save figure and close matlab

% Get path for figure directory
[directory,filename,~] = fileparts(filename);
filename = fullfile(directory,'Figures',filename);

% Save figure window
saveas(h.window,[filename '.fig'])
exit


end
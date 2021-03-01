function [e,on_point] = exp_rig_set(s,N,des,err,r_mid)
% EXP_RIG_SET  Open a plotting window to monitor and set rig flow coefficient
%
%   [] = EXP_RIG_SET(s,N,des,r_mid)
%
%   s - data structure with open instrument connections
%   N - data structure of required channels
%   des - data structure of desired operating point
%   err - data structure of allowable error fractions
%   r_mid - midheight radius for flow coefficient

% Get air constants
exp_air

% Set point set flag as global
global g

% Size of plotting window
t_delay = 200;

% Open rig operating point setting window
h.window = figure(); set(h.window,'Position',[1 29 1920 985]);

% Initialise axes for operating points and raw pressures and temperatures
if isfield(des,'phi') == 1
    h.phi = subplot_tight(2,1,1,[0.06 0.06]); hold on; grid on; box on;
    xlabel('Time / s'); ylabel('Flow Coef.');  
    axis([-t_delay 0 des.phi*(1-err.phi*6) des.phi*(1+err.phi*6)]);
else
    h.mrtp = subplot_tight(2,2,1,[0.06 0.06]); hold on; grid on; box on;
    xlabel('Time / s'); ylabel('Flow Func.'); 
    dmrtp = err.mrtp*6; mrtp_min = des.mrtp*(1-dmrtp); mrtp_max = des.mrtp*(1+dmrtp);
    axis([-t_delay 0 mrtp_min mrtp_max]);
    h.nrt = subplot_tight(2,2,2,[0.06 0.06]); hold on; grid on; box on;
    xlabel('Time / s'); ylabel('Non. Dim. Speed');  
    dnrt = err.nrt*6; nrt_min = des.nrt*(1-dnrt); nrt_max = des.nrt*(1+dnrt);
    axis([-t_delay 0 nrt_min nrt_max]);
end
h.P_raw = subplot_tight(2,2,3,[0.06 0.06]); hold on; grid on; box on;
xlabel('Time / s'); ylabel('Raw Pressure');  axis([-t_delay 0 -7000 7000]);
h.T_raw = subplot_tight(2,2,4,[0.06 0.06]); hold on; grid on; box on;
xlabel('Time / s'); ylabel('Raw Temperature'); axis([-t_delay 0 280 310]);

% Plot desired operating points and errors
varnames = fieldnames(des);
for v = 1:length(varnames)
    axes(h.(varnames{v})); 
    plot([-t_delay 0],des.(varnames{v})*ones(1,2),'k-');
    plot([-t_delay 0],des.(varnames{v})*ones(1,2) * (1 + err.(varnames{v})),'r-');
    plot([-t_delay 0],des.(varnames{v})*ones(1,2) * (1 - err.(varnames{v})),'r-');
end

% Get current rig conditions
[e,on_point] = exp_rig_monitor(s,[],N,des,err,r_mid);
P_chan = length(e.P); T_chan = length(e.T);

% Initialise plotting variables
nt = round(t_delay*5);
time_diff = -200*ones(nt,1); time_raw = zeros(nt,6);
if isfield(des,'phi') == 1
    phi = e.phi*ones(nt,1);
else
    mrtp = e.mrtp*ones(nt,1);
    nrt = e.nrt*ones(nt,1);
end
P = repmat(e.P,[nt 1]);
T = repmat(e.T,[nt 1]);

% Pressure and temperature legends
P_legend = cell(P_chan,1); T_legend = cell(T_chan,1);
for n = 1:P_chan; P_legend{n} = num2str(n); end;
for n = 1:T_chan; T_legend{n} = num2str(n); end;
P_cols = distinguishable_colors(P_chan);
T_cols = distinguishable_colors(P_chan);

% Plot initialised values
if isfield(des,'phi') == 1
    axes(h.phi); l.phi = line(time_diff,phi); 
else
    axes(h.mrtp); l.mrtp = line(time_diff,mrtp); 
    axes(h.nrt); l.nrt = line(time_diff,nrt); 
end
axes(h.P_raw); for n = 1:P_chan; l.P_raw(n) = line(time_diff,P(:,n),'color',P_cols(n,:)); end;
legend(P_legend,'Location','West','FontSize',6);
axes(h.T_raw); for n = 1:T_chan; l.T_raw(n) = line(time_diff,T(:,n),'color',T_cols(n,:)); end;
legend(T_legend,'Location','West','FontSize',6);

% Button readings for current operating points
vw = get(h.window,'Position');
if isfield(des,'phi') == 1
    va = get(h.phi,'Position'); 
    h.phi_indicator = uicontrol('style','push','position',...
        round([va(1)*vw(3)+5 va(2)*vw(4)+5 70 40]),...
        'string',num2str(phi(1)));
else
    va = get(h.mrtp,'Position'); 
    h.mrtp_indicator = uicontrol('style','push','position',...
        round([va(1)*vw(3)+5 va(2)*vw(4)+5 70 40]),...
        'string',num2str(mrtp(1)));
    va = get(h.nrt,'Position'); 
    h.nrt_indicator = uicontrol('style','push','position',...
        round([va(1)*vw(3)+5 va(2)*vw(4)+5 70 40]),...
        'string',num2str(nrt(1)));
end

% Push button for set flow coefficient
g = 0; 
h.button = uicontrol('style','push','position',[20 20 60 20],...
    'string','Point Set','callback','global g; g=g+1;');

% Indicator led for on point
if on_point == 0
    set(h.button,'CData',cat(3,ones(100),zeros(100),zeros(100)))
else
    set(h.button,'CData',cat(3,zeros(100),ones(100),zeros(100)))
end

% Read values until throttle is set
while g == 0    
    
    % Read and calculate current rig conditions
    [e,on_point] = exp_rig_monitor(s,[],N,des,err,r_mid);
    
    % Update all histories
    P = [e.P ; P(1:end-1,:)];
    T = [e.T ; T(1:end-1,:)];
    if isfield(des,'phi') == 1
        phi = [e.phi ; phi(1:end-1)];
    else
        nrt = [e.nrt ; nrt(1:end-1)];
        mrtp = [e.mrtp ; mrtp(1:end-1)];
    end
    
    % Read single time
    time_read = clock;
    
    % Elapsed time
    time_raw = [time_read ; time_raw(1:end-1,:)];
    for n = 1:length(time_diff)
        time_diff(n,:) = etime(time_raw(n,:),time_read);
    end
    
    % Update all time series traces
    if isfield(des,'phi') == 1
        set(l.phi,'XData',time_diff,'YData',phi);
    else
        set(l.mrtp,'XData',time_diff,'YData',mrtp);
%         set(h.mrtp,'ylim',[min(mrtp_min,(1-dmrtp)*mrtp(1)) max(mrtp_max,(1+dmrtp)*mrtp(1))])
        set(l.nrt,'XData',time_diff,'YData',nrt);
%         set(h.nrt,'ylim',[min(nrt_min,(1-dnrt)*nrt(1)) max(nrt_max,(1+dnrt)*nrt(1))])
    end
    for n = 1:length(l.P_raw)
        set(l.P_raw(n),'XData',time_diff,'YData',P(:,n));
    end
    for n = 1:length(l.T_raw)
        set(l.T_raw(n),'XData',time_diff,'YData',T(:,n));
    end
  
    % Update indicator values
    if isfield(des,'phi') == 1
        set(h.phi_indicator,'string',num2str(phi(1)));
    else
        set(h.mrtp_indicator,'string',num2str(mrtp(1)));
        set(h.nrt_indicator,'string',num2str(nrt(1)));
    end
    
    % Update LED colour
    if on_point == 0
        set(h.button,'CData',cat(3,ones(100),zeros(100),zeros(100)))
    else
        set(h.button,'CData',cat(3,zeros(100),ones(100),zeros(100)))
    end
    
    % Check if flow coefficient set button has been pressed
    drawnow;

    % Build in delay
    pause(0.3);
end

% Close window
close(h.window);
pause(1);


end
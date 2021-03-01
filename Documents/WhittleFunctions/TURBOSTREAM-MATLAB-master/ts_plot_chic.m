function [inlets,outlets,c,h] = ts_plot_chic(directory,file_prefix,h,use_avg,coords,col,force_read,speed)
% TS_PLOT_CHIC  Plot a constant speed compressor pressure rise and efficiency characteristic
%
%   [inlets,outlets,c,h] = TS_PLOT_CHIC(directory,file_prefix,h,use_avg,coords,col)
%
%   directory - string of directory to search in
%   file_prefix - string of filenames to match
%   h - figure window handle for plotting
%   use_avg - 0 or 1 to select average TURBOSTREAM solutions
%   coords - matrix or cell array to choose cut planes
%   col - RGB colour vector
%   force_read - 0 or 1 to read whole file and take new cuts
%   inlets - cell array of inlet cuts
%   outlets - cell array of outlet cuts
%   c - struct with calculated 1D performance values

% Default to open new window
if exist('h','var') == 0 || isempty(h) == 1
    h = figure();
end

% Default to use average solutions
if exist('use_avg','var') == 0
    use_avg = 1;
end

% Default colour to black
if exist('col','var') == 0 || isempty(col) == 1
    col = [0 0 0];
end

% Default to read cut file
if exist('force_read','var') == 0
    force_read = 0;
end

% Default to high speed plotting
if exist('speed','var') == 0
    speed = 1;
end

% Set figure window size
figure(h); set(h,'Position',[1 82 700 895])
    

%% Process filenames, read solutions and take cuts

% Get all filenames
[F,I] = ts_get_files([directory file_prefix],use_avg);

% Preallocate cuts and back pressure vector
inlets = cell(length(F),1);
outlets = cell(length(F),1);
pouts = zeros(length(F),1);
c = [];

% Loop over all files of interest and pull out cuts
for n = 1:length(F)

    % Check if the calculation is stalling
    logname = ['log_' strrep(strrep(F{n},'.hdf5','.txt'),'_avg','')];
    if exist([directory logname],'file') ~= 0 && ts_check_stall([directory logname]) == 1
        disp([logname ' stalling'])
%         delete([directory F{n}]);
%         delete([directory strrep(F{n},'.hdf5','.mat')]);
%         continue
    end

    % Take inlet and outlet cuts
    [inlet,outlet,pout,rpm] = ts_cut_chic([directory F{n}],coords,force_read);
    
    % Store the cuts in the cell array
    outlets{n} = outlet; inlets{n} = inlet; pouts(n) = pout;
    
end

% Return if no files are found
if isempty(F) == 1; return; end;

% Remove empty parts of the cell array
q = pouts == 0;
pouts(q) = []; F(q) = []; inlets(q) = []; outlets(q) = []; 
varnames = fieldnames(I); for v = 1:length(varnames); I.(varnames{v})(q) = []; end;


%% Calculate 1D performance between cut planes and plot the characteristic
    
% Calculate performance metrics at all points
for n = 1:length(inlets)

    % Calculate performance metrics
    c_temp = ts_calc_chic(inlets{n},outlets{n},rpm);
   
    % Assemble into vectors
    varnames = fieldnames(c_temp);
    for v = 1:length(varnames)
        c.(varnames{v})(n,1) = c_temp.(varnames{v});
    end
    
end

% Sort data
[~,i] = sort(c.mrtp_in,1,'descend'); c.F = F; c.pout = pouts;  
varnames = fieldnames(c); for v = 1:length(varnames); c.(varnames{v}) = c.(varnames{v})(i); end;
varnames = fieldnames(I); for v = 1:length(varnames); c.I.(varnames{v}) = I.(varnames{v})(i); end;
inlets = inlets(i); outlets = outlets(i); 

% Determine whether to plot high speed or low speed measures
if c.Rpo(1) < 1.1 || speed == 0
    flowname = 'phi'; flowlabel = 'Flow Coefficient'; 
    pname = 'Cptt'; plabel = 'Stagnation Pressure Rise Coefficient';
else
    flowname = 'mrtp_in'; flowlabel = 'Inlet Flow Function'; 
    pname = 'Rpo'; plabel = 'Stagnation Pressure Ratio';
end

% Plot the pressure rise characteristic
figure(h); subplot_tight(2,1,2,[0.05 0.1]); hold on; grid on; box on;
plot(c.(flowname),c.(pname),'.-','Color',col)
xlabel(flowlabel); ylabel(plabel);

% Plot the efficiency characteristic
subplot_tight(2,1,1,[0.03 0.1]); hold on; grid on; box on; set(gca,'xticklabel',[]);
plot(c.(flowname),c.eta_poly,'.-','Color',col)
ylabel('Total-Total Polytropic Efficiency');


end






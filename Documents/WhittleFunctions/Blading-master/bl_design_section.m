function [c,h] = bl_design_section(directory,c,f,p_target,target_type,h,plot_stuff)
% BL_DESIGN_SECTION  Design a section iteratively in MISES for given flow parameters

% Set default to not plot
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Set default to absolute target
if exist('target_type','var') == 0
    target_type = 'abs';
end

% Make a directory for the MISES solutions if it doesn't already exist
if exist(directory,'dir') == 0
    mkdir(directory);
end

% Select gradient type
grad_type = 'double';


%% Define iteration parameters

% Iteration limit
n_lim = 5; cols = lines(n_lim);

% Set maximum movements per step
o_master.chi_le = 2; 
o_master.dcam_le = 0.5;
o_master.s_thick_max = 0.02;
o_master.dcam_te = 0.02;
o_master.chi_te = 2;
o_master.N = 7;

% Step sizes to define gradients to input parameters
d_master.chi_le = 1; 
d_master.dcam_le = 0.05;
d_master.s_thick_max = 0.05;
d_master.dcam_te = 0.02; 
d_master.chi_te = 0.2;
d_master.N = 2;

% Set error weightings
w_master.psi_stag = 1; 
w_master.s_Cp_max = 100; 
w_master.Alpha = 10; 
w_master.aHb_1 = 250;
w_master.aHb_2 = 250;
w_master.Hb_te = 100;

% Convergence tolerances
t_master.psi_stag = 0.1; 
t_master.s_Cp_max = 0.005; 
t_master.Alpha = 0.05; 
t_master.aHb_1 = 0.0005;
t_master.aHb_2 = 0.0005;
t_master.Hb_te = 0.001;

% Characterise design operation type
if length(fieldnames(p_target)) == 5

    % Design detailed section parameters
    fl_names = {'psi_stag' 's_Cp_max' 'aHb_1' 'aHb_2' 'Alpha'};
    fl_text = 'Psi = %3.1f s = %4.3f aHbar = %5.4f,%5.4f Alpha = %3.1f';
    bl_names = {'chi_le' 'dcam_le' 's_thick_max' 'dcam_te' 'chi_te'};
    bl_text = 'chi_le = %2.1f dcam_le = %2.2f s_thick_max = %2.2f dcam_te = %2.2f chi_te = %2.1f';
    
elseif length(fieldnames(p_target)) == 2
    
    % Design exit flow angle and trailing edge shape factor
    fl_names = {'Hb_te' 'Alpha'}; fl_text = 'Hb = %3.2f Alpha = %3.1f';
    bl_names = {'chi_te' 'N'}; bl_text = 'chi_te = %2.1f N = %i';
    
elseif length(fieldnames(p_target)) == 1
    
    % Design exit flow angle only
    fl_names = {'Alpha'}; fl_text = 'Alpha = %3.1f';
    bl_names = {'chi_te'}; bl_text = 'chi_te = %2.1f';
end

% Remove excess variables from data structures
for v = 1:length(fl_names)
    w.(fl_names{v}) = w_master.(fl_names{v}); t.(fl_names{v}) = t_master.(fl_names{v});
end
for v = 1:length(bl_names)
    o.(bl_names{v}) = o_master.(bl_names{v}); d.(bl_names{v}) = d_master.(bl_names{v});
end

%% Iterate on sections paramters until tolerances are met or limit is exceeded

% Print current directory
fprintf('Loc      :  %s\n',directory)

% Run first iteration and print targets
[p,h] = mis_run_section(directory,c,f,h,[0 0 0],plot_stuff);

% Set relative target if specified so
if strcmp(target_type,'rel') == 1
    
    % Calculate differences
    p_rel = p_target;
    for v = 1:length(fl_names)
        p_target.(fl_names{v}) = p.(fl_names{v}) + p_rel.(fl_names{v});
    end
    
    % Check current target is physical
    if isfield(p_target,'s_Cp_max') == 1 && p_target.s_Cp_max < 0
        p_target.s_Cp_max = p.s_Cp_max;
    end
end

% Print target and current flow values
for v = 1:length(fl_names); fl_temp(v) = p.(fl_names{v}); end;
fprintf(['Current  :  ' fl_text '\n'],fl_temp)
for v = 1:length(fl_names); fl_temp(v) = p_target.(fl_names{v}); end;
fprintf(['Target   :  ' fl_text '\n'],fl_temp)

% Print current design
for v = 1:length(bl_names); bl_temp(v) = c.(bl_names{v}); end;
fprintf(['Design   :  ' bl_text '\n'],bl_temp)

% Iterate and fix each problem in turn if necessary
n = 1; conv_flag = 0;
while n <= n_lim && conv_flag == 0

    % Remove non-design variable names
    varnames = fieldnames(p); 
    for v = 1:length(varnames); 
        if sum(strcmp(varnames{v},fl_names)) == 0
            p = rmfield(p,varnames{v});
        end
    end
    
    % Get gradients continuously for a more accurate solution
    g = get_gradients(directory,c,d,f,p,grad_type);
    
    % Match flow parameters with camberline
    [c,p] = gradient_design(c,p_target,g,w,f,o,p,directory,h,cols(n,:),plot_stuff);

    % Print current design
    for v = 1:length(bl_names); bl_temp(v) = c.(bl_names{v}); end;
    fprintf(['Design   :  ' bl_text '\n'],bl_temp)
    
    % Calculate errors for convergence check
    e = zeros(length(fl_names),1); t_conv = zeros(length(fl_names),1);
    for v = 1:length(fl_names)
        e(v) = p_target.(fl_names{v}) - p.(fl_names{v});
        t_conv(v) = t.(fl_names{v});
    end
    
    % Check if overall convergence is reached
    if sum(abs(e) < t_conv) == length(fl_names)
        conv_flag = 1;
    end
    
    % Increment iteration count
    n = n+1;

end

% Print final flow achievements
for v = 1:length(fl_names); fl_temp(v) = p.(fl_names{v}); end;
fprintf(['Final    :  ' fl_text '\n\n'],fl_temp)


end

function c = param_design_1d(c,bl,fl,o,g,w,p,p_target)
% Construct matrix of predicted values from gradients

% Resolution to generate grid of possible results
ni = 101; ei = (1:ni) - round(ni/2);

% Construct vector of blade parameters
x = linspace(c.(bl)-o.(bl),c.(bl)+o.(bl),ni);

% Values of flow from gradients
a = p.(fl) + g.(bl).(fl) * (x - c.(bl)) - p_target.(fl);

% Penalty function for choosing extreme results
e = ei.^2; e = 5 * e / max(e);

% Calculate weighted distance from the optimal
q = abs(a) .* w.(fl) + e;

% Find the minimum that is the optimal design
[~,i] = min(q); c.(bl) = x(i); 

end

function c = param_design_2d(c,bl_1,bl_2,fl_1,fl_2,o,g,w,p,p_target)
% Construct matrix of predicted values from gradients

% Resolution to generate grid of possible results
ni = 101; ei = (1:ni) - round(ni/2);

% Construct matrices of blade parameters
x = repmat(reshape(linspace(c.(bl_1)-o.(bl_1),c.(bl_1)+o.(bl_1),ni),[ni 1]),[1 ni]);
y = repmat(reshape(linspace(c.(bl_2)-o.(bl_2),c.(bl_2)+o.(bl_2),ni),[1 ni]),[ni 1]);

% Values of flow from gradients
a = p.(fl_1) + g.(bl_1).(fl_1) * (x - c.(bl_1)) + ...
    g.(bl_2).(fl_1) * (y - c.(bl_2)) - p_target.(fl_1);
b = p.(fl_2) + g.(bl_1).(fl_2) * (x - c.(bl_1)) + ...
    g.(bl_2).(fl_2) * (y - c.(bl_2)) - p_target.(fl_2);

% Penalty function for choosing extreme results
[e1,e2] = ndgrid(ei,ei); e = e1.^2 + e2.^2;
e = 5 * e / max(max(e));

% Calculate weighted distance from the optimal
q = abs(a).*w.(fl_1) + abs(b).*w.(fl_2) + e;

% Find the minimum that is the optimal design
[q1,i] = min(q,[],1); [~,j] = min(q1,[],2); i = i(j);
c.(bl_1) = x(i,j); c.(bl_2) = y(i,j);


end

function [c,p] = gradient_design(c,p_target,g,w,f,o,p,sec_directory,h,col,plot_stuff)
% Design a section in multiple stages

% Local incidence and position of peak suction
if isfield(w,'psi_stag') == 1 && isfield(w,'s_Cp_max') == 1
    c = param_design_2d(c,'chi_le','dcam_le','psi_stag','s_Cp_max',o,g,w,p,p_target);
    p = mis_run_section(sec_directory,c,f,[],[],0);
    %%% NOTE min(c.dcam_le+o.dcam_le,0)  has been removed
end

% Shape factor distribution
if isfield(w,'aHb_1') == 1 && isfield(w,'aHb_2') == 1
    c = param_design_2d(c,'s_thick_max','dcam_te','aHb_1','aHb_2',o,g,w,p,p_target);
    p = mis_run_section(sec_directory,c,f,[],[],0);
    
end

% Exit gas angle and trailing edge shape factor
if isfield(w,'Alpha') == 1 && isfield(w,'Hb_te') == 1
    c = param_design_2d(c,'chi_te','N','Hb_te','Alpha',o,g,w,p,p_target);
    c.N = round(c.N);
    p = mis_run_section(sec_directory,c,f,h,col,plot_stuff);
end

% Exit gas angle only
if isfield(w,'Alpha') == 1 && isfield(w,'Hb_te') == 0
    c = param_design_1d(c,'chi_te','Alpha',o,g,w,p,p_target);
    p = mis_run_section(sec_directory,c,f,h,col,plot_stuff);
end

% Pause to show plots
if plot_stuff == 1
    pause(1);
end


end

function g = get_gradients(directory,c,d,f,p,grad_type)
% Calculate two sided gradients for all parameters

% Get blade and flow variable names
bl_names = fieldnames(d); fl_names = fieldnames(p);

% Calculate all derivatives
for n = 1:length(bl_names)
    
    % Perturb section design by steps in both directions
    c_1 = c; c_1.(bl_names{n}) = c_1.(bl_names{n}) - d.(bl_names{n});    
    c_2 = c; c_2.(bl_names{n}) = c_2.(bl_names{n}) + d.(bl_names{n});
    
    % Run MISES on both of the new section designs if double sided gradients are selected
    p_1 = mis_run_section(directory,c_1,f,[],[],0);
    if strcmp(grad_type,'double') == 1
        p_2 = mis_run_section(directory,c_2,f,[],[],0); s = 2;
    else
        p_2 = p; s = 1;
    end
    
    % Calculate derivatives for each flow variables
    for m = 1:length(fl_names)
        g.(bl_names{n}).(fl_names{m}) = (p_2.(fl_names{m}) - p_1.(fl_names{m})) / (s * d.(bl_names{n}));
    end
end


end
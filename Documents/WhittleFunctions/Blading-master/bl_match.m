function [B,p] = bl_match(dr,datum_name,new_name,P,rows,job,run_max,plot_stuff)
% Run TURBOSTREAM iteratively to match a multi-stage high speed compressor 

% Set default number of iterations
if exist('run_max','var') == 0
    run_max = 4;
end

% Set default to plot
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

% Number of rows in calculations
nn = length(rows);

% Sections to design on
r_secs = [0 0.2 0.5 0.8 1]; 
nj = 50; r_nondim = linspace(0,1,nj);

% Convergence limits
tol.Alpha = 0.1; tol.mfunc = 0.01;

% Define colours for iterations
cols = lines(run_max+1);

% Always run SA model only
job.run_sa = 1;


%% Plot datum flow in 3D and 2D if specified

% Open figure window for spline definitions and TURBOSTREAM
h.geom = cell(nn,1);
for n = 1:nn
    h.geom{n} = figure(); h.blade{n} = figure();
end

% Plot the datum flow
[p_datum,~,C_datum] = bl_plot_run(dr,datum_name,1:nn,rows,r_secs,h,[0 0 0],plot_stuff);

% Plot the datum spline geometry
for n = 1:nn
    load([dr.geom datum_name '_B' num2str(n) '_spline.mat'])
    bl_spline_plot(b,h.geom{n},[0 0 0]);
    subtitle(rows{n});
end

% Plot the Mach triangles
h.mach = figure(); hold on; grid on; box on; axis equal;
plot_triangles(C_datum,h.mach,[0 0 0])


%% Prepare for TURBOSTREAM design loop

% Search for all files with the current name
A = dir([dr.ts new_name '_R*_avg*']);
filenames = cell(length(A),1);
for n = 1:length(A);
    filenames{n} = A(n).name;
end

% Get the latest run numbers from filenames
splits = regexp(filenames,'_','split');
run_nos = [];
for n = 1:length(splits)
    for m = 1:length(splits{n})
        if strcmp(splits{n}{m}(1),'R') == 1
            run_nos = [run_nos ; str2double(splits{n}{m}(2:end))];
        end
    end
end
run_nos(isnan(run_nos)) = [];
run_nos = unique(run_nos);

% Check whether design loop has begun
if isempty(run_nos) == 1
    
    % Start the design process
    run_no = 1;

    % Save the initial parameters to file
    save([dr.geom new_name '_R1_param.mat'],'P')
    
    % Plot the spline geometry
    for n = 1:nn
        bl_spline_plot(P.B{n},h.geom{n},cols(run_no,:));
    end

    % Run TURBOSTREAM on the initial parameters
    job = bl_create_machine(P.B,dr,[new_name '_R1_'],job,P.I,P.M,P.F,P.D,0,1);
    [p,~,C] = bl_plot_run(dr,[new_name '_R1'],1:nn,rows,r_secs,h,cols(1,:),plot_stuff);
    
    % Plot the Mach triangles
    plot_triangles(C,h.mach,cols(m,:))
    
    % Delete unused files
    delete([job.directory job.outname]); delete([job.directory job.outxdmf]);
    
else
    
    % Get the last run number
    run_no = run_nos(end);
    
    % Loop over all previous runs and plot
    for m = 1:run_no
        
        % Plot TURBOSTREAM results
        [p,~,C] = bl_plot_run(dr,[new_name '_R' num2str(m)],1:nn,rows,r_secs,h,cols(m,:),plot_stuff);
    
        % Plot the spline geometry
        for n = 1:nn
            load([dr.geom new_name '_R' num2str(m) '_B' num2str(n) '_spline.mat'])
            bl_spline_plot(b,h.geom{n},cols(m,:));
        end
        
        % Plot the Mach triangles
        plot_triangles(C,h.mach,cols(m,:))
    end
end


%% Run the TURBOSTREAM design loop

% Initialise convergence flag
conv_flag = 0;

% Begin design loop
while run_no <= run_max && conv_flag == 0
    
    % Get blade spline data from last run
    B = cell(nn,1);
    for n = 1:nn
        load([dr.geom new_name '_R' num2str(run_no) '_B' num2str(n) '_spline.mat'])
        B{n} = b;
    end
    
    % Find relative whirl angles
    dAlpha = zeros(nn,nj);
    for n = 1:nn
            
        % Calculate differences between current and target
%         dAlpha(n,:) = p_datum{n}.Alpha - p{n}.Alpha;
        dAlpha(n,:) = interp1(p_datum{n}.radial.r_nondim,p_datum{n}.radial.Alpha,r_nondim,'linear',...
            'extrap') - interp1(p{n}.radial.r_nondim,p{n}.radial.Alpha,r_nondim,'linear','extrap');
            
        % Flip relative exit angle if already negative in the case of rotor rows
        if mean(p_datum{n}.Alpha) < 0
            dAlpha(n,:) = - dAlpha(n,:);
        end
    end
    
    % Update blade exit angles
    for n = 1:nn
        
        % Exclude angle changes in the endwall boundary layers
        i = find(r_nondim > p{n}.r_bl(1),1,'first'); dAlpha(n,1:i) = dAlpha(n,i);
        i = find(r_nondim < p{n}.r_bl(2),1,'last'); dAlpha(n,i:end) = dAlpha(n,i);
        
        % Calculate new trailing edge angle
        chi_te = fnval(B{n}.chi_te,r_nondim) + dAlpha(n,:);
        
        % Fit trailing edge angle spline
        B{n}.chi_te = spap2(B{n}.chi_te.knots,B{n}.chi_te.order,r_nondim,chi_te);
    end
    
    % Increment run counter
    run_no = run_no + 1;
    
    % Plot the new spline geometry
    for n = 1:nn
        bl_spline_plot(B{n},h.geom{n},cols(run_no,:));
    end

    % Run TURBOSTREAM on the new design
    [~,~] = system(['rm ' dr.ag new_name '_R' num2str(run_no) '*']);
    job = bl_create_machine(B,dr,[new_name '_R' num2str(run_no) '_'],job,P.I,P.M,P.F,P.D,0,1);

    % Delete unused files
    delete([job.directory job.outname]); delete([job.directory job.outxdmf]);
    
    % Plot latest iteration
    [p,~,C] = bl_plot_run(dr,[new_name '_R' num2str(run_no)],1:nn,rows,r_secs,h,cols(run_no,:),...
        plot_stuff);
    plot_triangles(C,h.mach,cols(run_no,:))
    
    % Check convergence on later iterations
    if max(abs(dAlpha)) < tol.Alpha
        conv_flag = 1;
    end
end

% Move last iteration to final name
copyfile([dr.ts new_name '_R' num2str(run_no) '_avg.mat'],[dr.ts new_name '_avg.mat']);
copyfile([dr.ts 'log_' new_name '_R' num2str(run_no) '.txt'],[dr.ts 'log_' new_name '.txt']);

% Copy last geometry to final name
copyfile([dr.geom new_name '_R' num2str(run_no) '_xr_hub.mat'],[dr.geom new_name '_xr_hub.mat']);
copyfile([dr.geom new_name '_R' num2str(run_no) '_xr_cas.mat'],[dr.geom new_name '_xr_cas.mat']);
for n = 1:nn
    copyfile([dr.geom new_name '_R' num2str(run_no) '_B' num2str(n) '_spline.mat'],...
        [dr.geom new_name '_B' num2str(n) '_spline.mat']);
end

% Run characteristic on final design
job.name = new_name;
ts_run_chic(job,plot_stuff,[],0,0.001)


end

function [] = plot_triangles(C,h,col)
% Plot the midspan Mach triangles for comparison

% Change to the figure window
figure(h);

% Rotational speed
w = 2 * pi * C{1}.inlet.bv.rpm / 60;

% Gas constants
cp = C{1}.inlet.av.cp; ga = C{1}.inlet.av.ga; cv = cp / ga; R = cp - cv;

for n = 1:length(C)+1

    % Take correct cut
    if n == 1
        c = C{1}.inlet;
    else
        c = C{n-1}.outlet;
    end

    % Calculate midspan mass averages
    M = interp1(ts_mass_average(c,'r_nondim',2),ts_mass_average(c,'M',2),0.5);
    Alpha = interp1(ts_mass_average(c,'r_nondim',2),ts_mass_average(c,'Alpha',2),0.5);
    r = interp1(ts_mass_average(c,'r_nondim',2),ts_mass_average(c,'r',2),0.5);
    T = interp1(ts_mass_average(c,'r_nondim',2),ts_mass_average(c,'T',2),0.5);

    % Blade Mach number
    M_blade = r * w / (ga * R * T)^0.5;

    % Plot Mach triangles
    plot([0 M*cosd(Alpha)] + n,[0 M*sind(Alpha)],'-','color',col);
    plot(M*cosd(Alpha) * [1 1] + n,M*sind(Alpha) + [0 -M_blade],'-','color',col)
    plot([0 M*cosd(Alpha)] + n,[0 M*sind(Alpha) - M_blade],'-','color',col)

end

end
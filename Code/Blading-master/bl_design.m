function [B,p] = bl_design(dr,datum_name,new_name,P,n_des,rows,job,run_max,plot_stuff)
% Couple TURBOSTREAM and MISES to match a design of a multi-stage, high speed compressor

% Set default number of iterations
if exist('run_max','var') == 0
    run_max = 4;
end

% Set default to plot
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

% Sections to design on
r_secs = 0:0.1:1; nj = length(r_secs);

% Define limits in changes to camber parameters per iteration
l.chi_le = 5; l.chi_te = 3; l.dcam_le = 0.4; l.dcam_te = 0.1; l.s_thick_max = 0.05; l.N = 10;

% Define convergence tolerances in objective functions
tol.psi_stag = 1; tol.s_Cp_max = 0.002; tol.Alpha = 0.1; 
tol.aHb_1 = 0.0005; tol.aHb_2 = 0.0005; tol.Hb_te = 0.01;

% Rows to design
nn = length(n_des);

% Define colours for iterations
cols = lines(run_max);

% Always run SA model only
job.run_sa = 1;

% Names for data structures
misnames = {'M' 'Alpha' 'AVDR' 'Re'};
% tarnames = {'psi_stag' 's_Cp_max' 'aHb_1' 'aHb_2' 'Hb_te' 'Alpha'};
tarnames = {'psi_stag' 's_Cp_max' 'aHb_1' 'aHb_2' 'Alpha'};
blnames = {'chi_le' 'chi_te' 's_thick_max' 'dcam_le' 'dcam_te'};


%% Plot datum flow in 3D and 2D if specified

% Open figure window for spline definitions and TURBOSTREAM
h.geom = cell(nn,1);
for n = n_des
    h.geom{n} = figure(); h.blade{n} = figure();
end

% Plot the datum flow
p_datum = bl_plot_run(dr,datum_name,n_des,rows,r_secs,h,[0 0 0],plot_stuff);

% Plot the datum spline geometry
for n = n_des
    load([dr.geom datum_name '_B' num2str(n) '_spline.mat'])
    bl_spline_plot(b,h.geom{n},[0 0 0]);
    subtitle(rows{n});
end


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
    for n = n_des
        bl_spline_plot(P.B{n},h.geom{n},cols(run_no,:));
    end

    % Run TURBOSTREAM on the initial parameters
    job = bl_create_machine(P.B,dr,[new_name '_R1_'],job,P.I,P.M,P.F,P.D,0,1);
    [p,f] = bl_plot_run(dr,[new_name '_R1'],n_des,rows,r_secs,h,cols(1,:),plot_stuff);
    
    % Delete unused files
    delete([job.directory job.outname]); delete([job.directory job.outxdmf]);
    
else
    
    % Get the last run number
    run_no = run_nos(end);
    
    % Loop over all previous runs and plot
    for m = 1:run_no
        
        % Plot TURBOSTREAM results
        [p,f] = bl_plot_run(dr,[new_name '_R' num2str(m)],n_des,rows,r_secs,h,cols(m,:),plot_stuff);
    
        % Plot the spline geometry
        for n = n_des
            load([dr.geom new_name '_R' num2str(m) '_B' num2str(n) '_spline.mat'])
            bl_spline_plot(b,h.geom{n},cols(m,:));
        end
    end
end


%% Run the TURBOSTREAM design loop

% Initialise convergence flag
conv_flag = 0;

% Begin design loop
while run_no <= run_max && conv_flag == 0
    
    % Print progress to screen
    if plot_stuff == 1
        disp(['%% Designing Sections on R' num2str(run_no)])
    end
    
    % Get blade spline data from last run
    B = cell(nn,1);
    for n = n_des
        load([dr.geom new_name '_R' num2str(run_no) '_B' num2str(n) '_spline.mat'])
        B{n} = b;
    end
        
    % Get section spline definitions from blade spline definition
    C = cell(nn,nj);
    for n = 1:nn
        for j = 1:nj
            C{n,j} = bl_spline_eval(B{n},r_secs(j));
        end
    end
    
    % Find relative targets
    p_target = cell(nn,nj);
    for n = 1:nn
        for j = 1:nj
            
            % Calculate differences between current and target
            for v = 1:length(tarnames)
                p_target{n,j}.(tarnames{v}) = p_datum{n}.(tarnames{v})(j) - p{n}.(tarnames{v})(j);
            end
            
            % Flip relative exit angle if already negative in the case of rotor rows
            if p_datum{n}.Alpha(j) < 0
                p_target{n,j}.Alpha = - p_target{n,j}.Alpha;
            end
        end
    end

    % Design sections in MISES to achieve the delta's calculated from TURBOSTREAM
    C_new = cell(nn,nj);
    for n = 1:nn
        
        % Choose blade number based upon midspan trailing edge shape factor
%         dr.sec = [dr.mis new_name '/' rows{n} '/Jmid/R' num2str(run_no) '/'];
%         clear p_sec; p_sec.Alpha = p_target{n,j_mid}.Alpha; p_sec.Hb_te = p_target{n,j_mid}.Hb_te;
%         for v = 1:length(misnames); f_sec.(misnames{v}) = f{n}.(misnames{v})(j_mid); end;
%         c_new = bl_design_section(dr.sec,C{n,j_mid},f_sec,p_sec,'rel',[],0); N = c_new.N;
        
        % Design sections for other parameters
        for j = 1:nj
            
            % Construct unique directory for each section
            dr.sec = [dr.mis new_name '/' rows{n} '/J' num2str(r_secs(j)*100) '/R' num2str(run_no) '/'];
            
            % Extract MISES boundary conditions
            for v = 1:length(misnames); f_sec.(misnames{v}) = f{n}.(misnames{v})(j); end;
            
            % Extract target flow values
            clear p_sec; 
            if run_no == 1
                p_sec.Alpha = p_target{n,j}.Alpha;
            else
                for v = 1:length(tarnames)
                    p_sec.(tarnames{v}) = p_target{n,j}.(tarnames{v});
                end
            end
            
            % Update blade number
%             C{n,j}.N = N; p_sec = rmfield(p_sec,'Hb_te');
            
            % Run section designing function
            C_new{n,j} = bl_design_section(dr.sec,C{n,j},f_sec,p_sec,'rel',[],0);
        end
    end

    % Update blade spline defintions with new sections
    for n = 1:nn
        for v = 1:length(blnames)
            
            % Record vectors of old and new parameters
            c = zeros(1,nj); c_new = zeros(1,nj);
            for j = 1:nj
                c(j) = C{n,j}.(blnames{v});
                c_new(j) = C_new{n,j}.(blnames{v});
            end
            
            % Calculate differences and apply limits to movements
            e = c_new - c; c_new = c + min(abs(e),l.(blnames{v})) .* sign(e);
            
            % Update blade spline definition
            B{n}.(blnames{v}) = spap2(B{n}.(blnames{v}).knots,B{n}.(blnames{v}).order,r_secs,c_new);
        end
        
        % Update blade number
%         B{n}.N = C_new{n,1}.N;
    end
    
    % Increment run counter
    run_no = run_no + 1;
    
    % Plot the new spline geometry
    for n = n_des
        bl_spline_plot(B{n},h.geom{n},cols(run_no,:));
    end
    
    % Run TURBOSTREAM on the new design
    [~,~] = system(['rm ' dr.ag new_name '_R' num2str(run_no) '*']);
    job = bl_create_machine(B,dr,[new_name '_R' num2str(run_no) '_'],job,P.I,P.M,P.F,P.D,0,1);

    % Delete unused files
    delete([job.directory job.outname]); delete([job.directory job.outxdmf]);
    
    % Plot latest iteration
    [p,f] = bl_plot_run(dr,[new_name '_R' num2str(run_no)],n_des,rows,r_secs,h,cols(run_no,:),...
        plot_stuff);
    
    % Check convergence on later iterations
    if run_no > 2
    
        % Assemble convergence matrix
        conv_mat = zeros(nn,length(tarnames));
        for n = 1:nn
            for v = 1:length(tarnames)
                if max(abs(p{n}.(tarnames{v}) - p_target{n}.(tarnames{v}))) < tol.(tarnames{v})
                    conv_mat(n,v) = 1;
                end
            end
        end

        % Check if all parameters are converged
        if sum(sum(conv_mat)) == numel(conv_mat)
            conv_flag = 1;
        end
    else
        
        % Not converged on first pass
        conv_flag = 0;
    end
    
end

% Move last iteration to final name
copyfile([directory new_name '_R' num2str(run_no) '_avg.mat'],[directory new_name '_avg.mat']);
copyfile([directory 'log_' new_name '_R' num2str(run_no) '.txt'],[directory 'log_' new_name '.txt']);

% Run characteristic on final design
job.name = new_name;
ts_run_chic(job,plot_stuff,[],0,0.002)


end


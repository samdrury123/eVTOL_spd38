function job = ts_run_chic(job,coords,plot_stuff,h,external_monitor,phi_tol,col)
% TS_RUN_CHIC  Use TURBOSTREAM to complete a characteristic from a given file
%
%   job = TS_RUN_CHIC(filename,job,h,plot_stuff,external_monitor,tol_phi)

% Set default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Open new figure window if requires
if plot_stuff == 1 && (exist('h','var') == 0 || isempty(h) == 1)
    h = figure(); hold on; grid on; box on;
    xlabel('Flow Coefficient'); ylabel('Total-static Pressure Rise');
end

% Default colour
if plot_stuff == 1 && (exist('col','var') == 0 || isempty(col) == 1)
    col = lines(1);
end

% Monitor in function
if exist('external_monitor','var') == 0
    external_monitor = 0;
end

% Set flow coefficient resolution
if exist('phi_tol','var') == 0
    phi_tol = 0.002;
end

% Starting mass flow decrement and chic differences
dmass_start = 0.04; dmass = dmass_start; dphi = inf; m_start = 10000; i_current = 1;

% Set maximum iteration count
nstep_max = 200000;


%% Run function recursively to run multiple chics at the same time

% Check if array of jobs is input
if iscell(job) == 1 && external_monitor == 0
    
    % Define job array
    Jobs = job; Jobs_start = Jobs; finished = ones(size(Jobs)); chic_finished = zeros(size(Jobs));
    
    % Job colours
    cols = hsv(length(Jobs));
    
    % Continually call chic running function until all jobs are finished
    while all(chic_finished) == 0
        for n = 1:length(Jobs)

            % Continue characteristic if last job has finished
            if finished(n) ~= 0
                Jobs{n} = ts_run_chic(Jobs_start{n},coords,plot_stuff,h,1,phi_tol,cols(n,:));
            end

        end
        
        % Check for finished jobs
        [~,finished] = ts_smonitor(Jobs,1);
        
        % Check for finished chics
        chic_finished = zeros(size(Jobs));
        for n = 1:length(Jobs)
            if isfield(Jobs{n},'id') == 0
                chic_finished(n) = 1;
            end
        end
        
    end
    
    % Exit the function when all characteristics are complete
    return
    
end


%% Analyse the design point case

% Copy first file if name does not contain "M"
if isempty(strfind(job.name,['m' num2str(m_start)])) == 1
    
    % Copy average TURBOSTREAM file
    if exist([job.directory job.name '_avg.hdf5'],'file') ~= 0
        copyfile([job.directory job.name '_avg.hdf5'],[job.directory job.name '_M' ...
            num2str(m_start) '_avg.hdf5'])
    else 
        copyfile([job.directory job.name '_avg.mat'],[job.directory job.name '_M' ...
            num2str(m_start) '_avg.mat'])
    end
    
    % Copy logfile
    if exist([job.directory 'log_' job.name '.txt'],'file') ~= 0
        copyfile([job.directory 'log_' job.name '.txt'],[job.directory 'log_' job.name '_M' ...
            num2str(m_start) '.txt'])
    else
        copyfile([job.directory 'log_' job.name '.mat'],[job.directory 'log_' job.name '_M' ...
            num2str(m_start) '.mat'])
    end
    
    % Update name
    job.name = [job.name '_M' num2str(m_start)];
end

% Log file name
job.logname = ['log_' job.name '.txt'];

% Read previous run and analyse the chic
g = ts_read_hdf5([job.directory job.name '_avg.hdf5']);
[inlet,outlet,pout,rpm] = ts_cut_chic([job.directory job.name '_avg.hdf5'],coords);

% Record iteration count and walltime
nstep_start = double(g{1}.av.nstep); time_start = job.walltime;

% Calculate backpressure adjustment for first point
[~,nstep_conv(1),dpout] = ts_plot_conv([job.directory 'log_' job.name '.txt'],0);
pouts(1) = pout + dpout;

% Calculate characteristic for first point
c = ts_calc_chic(inlet,outlet,rpm);
phi(1) = c.phi; Cp(1) = c.Cpts;

% Calculate mass flow rate at exit
[~,mass_start] = ts_mass_average(outlet,'ro',3); mass_start = mass_start * double(outlet.bv.nblade);
mass_frac(1) = 1;


%% Run a characteristic of fixed mass flows

% Check if lines already exist on the plot and initialise
if plot_stuff == 1
    
    % Find handles of existing lines
    l_points = findobj(h,'color',col,'linestyle','-');
    l_interp = findobj(h,'color',col,'linestyle','--');
    
    % Plot the lines if they don't already exist
    if isempty(l_points) == 1; l_points = plot(phi,Cp,'.-','color',col); end;
    if isempty(l_interp) == 1; l_interp = plot(nan,nan,'--','color',col); end;
    
end

% Record file_prefix
file_prefix = strrep(job.name,['_M' num2str(m_start)],'');

% Run jobs until flow coefficient is within tolerance
while dphi > phi_tol && dmass > 1 / m_start
    
    % Interpolate to find current maximum pressure rise point
    if i_current ~= 1
        
        % Generate flow coefficient vector to interpolate
        phi_interp = phi(1); 
        for n = 1:length(phi)-1
            phi_interp = [phi_interp linspace(phi_interp(end),phi(n+1),100)];
        end
        phi_interp = unique(phi_interp);
        
        % Spline all points
        Cp_interp = interp1(phi,Cp,phi_interp,'spline');
        
        % Find the point previous to the maximum
        [~,i_max] = max(Cp_interp); phi_max = phi_interp(i_max);
        i_current = find(phi >= phi_max,1,'last');
        
        % Plot the interpolation
        if plot_stuff == 1
            set(l_interp,'xdata',phi_interp,'ydata',Cp_interp);
            drawnow; pause(0.1);
        end
    end
    
    % Run jobs until pressure rise falls
    dCp = 1;
    while dCp > 0
    
        % Reduce flow rate
        mass_frac(end+1,1) = mass_frac(i_current) - dmass;

        % Update job name
        current_name = [file_prefix '_M' num2str(round(mass_frac(i_current) * m_start))];
        job.name = [file_prefix '_M' num2str(round(mass_frac(end) * m_start))];
        job.logname = ['log_' job.name '.txt'];

        % Run the next point if it hasn't already been successful
        job.finished = 1;
        if exist([job.directory job.name '_avg.mat'],'file') == 0 && ...
            ((exist([job.directory job.name '_avg.hdf5'],'file') ~= 0 && ...
            H5F.is_hdf5([job.directory job.name '_avg.hdf5']) == 0) || ...
            exist([job.directory job.name '_avg.hdf5'],'file') == 0) || ...
            isempty(ts_plot_conv([job.directory 'log_' job.name '.txt'],0)) == 1

            % Read previous solution
            g = ts_read_hdf5([job.directory current_name '_avg.hdf5']);
            g{1}.av.poisson_nstep = 0; g{1}.av.poisson_restart = 1;

            % Update target mass flow and previous back pressure
            for bid = ts_get_bids(g)
                for pid = ts_get_pids(bid,g)
                    if g{bid+1}.patch{pid+1}.attribute.kind == 1
                        g{bid+1}.patch{pid+1}.pout = pouts(i_current); 
                        g{bid+1}.patch{pid+1}.throttle_target = mass_frac(end) * mass_start;
                        g{bid+1}.patch{pid+1}.throttle_k0 = 20;
                        g{bid+1}.patch{pid+1}.throttle_k1 = 50;
                        g{bid+1}.patch{pid+1}.throttle_k2 = 10;
                    end
                end
            end

            % Extrapolate number of iterations required for convergence
            if length(nstep_conv) > 2

                % Calculated estimated convergence point
%                 nstep_next = interp1(mass_frac(1:end-1),nstep_conv,mass_frac(end),'linear','extrap');
%                 nstep_next = ceil(nstep_next / 10000) * 10000 + 20000;
                nstep_next = ceil(nstep_conv(i_current) / 10000) * 10000 + 20000;

                % Update iteration number if required
                if nstep_next > double(g{1}.av.nstep)
                    g{1}.av.nstep = nstep_next; g{1}.av.nstep_save_start = nstep_next - 10000;
                elseif nstep_next > nstep_max
                    g{1}.av.nstep = nstep_max; g{1}.av.nstep_save_start = nstep_max - 10000;
                end

                % Adjust walltime
                job.walltime = ceil(time_start * double(g{1}.av.nstep) / nstep_start);

            end

            % Display job inputs
            disp(['Mass fraction = ' num2str(mass_frac(end)) '  nstep = ' num2str(g{1}.av.nstep)])

            % Submit job
            job = ts_sbatch(job,g);

            % Monitor job
            if external_monitor == 0
                job = ts_smonitor(job);
            else
                return
            end
        end

        % Check the calculation is not stalling
        disp(job.name);
        if ts_check_stall(job) == 1; job.finished = nan; end;

        % Check the calculation has converged
        [l_conv,nstep_conv_temp,dpout_temp] = ts_plot_conv([job.directory job.logname],0);
        disp(['Mass fraction = ' num2str(mass_frac(end)) ' Convergence at nstep = ' ...
            num2str(nstep_conv_temp) ' Current nstep = ' num2str(l_conv.nstep(end) + 50)])
    %     if nstep_conv_temp > g{1}.av.nstep; job.finished = nan; end;

        % Calculate the pressure rise and record if it is converged
        if isnan(job.finished) == 0

            % Attempt to read the hdf5 file
            try
                [inlet,outlet,pout,rpm] = ts_cut_chic([job.directory job.name '_avg.hdf5'],coords);
            catch
                delete([job.directory '*' job.name '*']);
                disp([job.name ' is corrupted']);
                continue;
            end
            
            % Calculate pressure rise
            c = ts_calc_chic(inlet,outlet,rpm);
            phi(end+1,1) = c.phi; Cp(end+1,1) = c.Cpts;
            
            % Record pressure rise value and update next point
            dCp = Cp(end) - Cp(end-1); i_current = length(Cp);
            
            % Record convergence information
            nstep_conv(end+1,1) = nstep_conv_temp; pouts(end+1,1) = pout + dpout_temp;
            
        else
            
            % If the job has failed exit the loop by setting pressure rise to fall
            mass_frac(end) = []; dCp = -1;
            
        end

        % Update characteristic plot
        if plot_stuff == 1
            set(l_points,'xdata',phi,'ydata',Cp)
            drawnow; pause(0.1);
        end

    end
    
    % Half the mass flow step
    dmass = 0.5 * dmass;
    
    % Sort the lists
    [phi,i] = unique(phi); phi = flip(phi,1); i = flip(i,1);
    Cp = Cp(i); pouts = pouts(i); nstep_conv = nstep_conv(i); mass_frac = mass_frac(i);

    % Update characteristic plot
    if plot_stuff == 1
        set(l_points,'xdata',phi,'ydata',Cp)
        drawnow; pause(1);
    end
    
    % Check flow coefficient tolerance
    dphi = (phi(1) - phi(2)) * dmass / dmass_start;
end


end


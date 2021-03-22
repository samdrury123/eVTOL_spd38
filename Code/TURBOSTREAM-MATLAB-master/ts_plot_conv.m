function [l,nstep_conv,dpout] = ts_plot_conv(filename,plot_stuff,plot_blocks,plot_poisson)
% TS_PLOT_CONV  Read in a TURBOSTREAM convergenece history and plot output
%
%   [l,nstep_conv,dpout] = ts_plot_conv(filename,plot_stuff,plot_blocks,plot_poisson)
%
%   filename - string of convergence file or job data structure
%   plot_stuff - 0 or 1 for showing general working
%   plot_blocks - 0 or 1 for showing block by block convergence
%   plot_poisson - 0 or 1 for showing poisson loop convergence
%   l - output convergence data structure
%   nstep_conv - estimated timestep of when convergence was or will be reached
%   dpout - pout correction if using throttle patch

% Set plotting defaults
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end
if exist('plot_blocks','var') == 0
    plot_blocks = 0;
end
if exist('plot_poisson','var') == 0
    plot_poisson = 0;
end

% Handle filename inputs for job structure
if isstruct(filename) == 1
    job = filename;
    if isfield(job,'logname') == 0; job.logname = ['log_' job.name '.txt']; end;
    filename = [job.directory job.logname];
end

% Make sure newest modified file is always read
A_txt = dir(filename); A_mat = dir(strrep(filename,'.txt','.mat'));
if numel(A_txt) == 1 && numel(A_mat) == 1 
    if A_txt.datenum > A_mat.datenum
        delete(strrep(filename,'.txt','.mat'))
    end
end

% Check to see if convergence mat file already exists and load data structure
if exist(strrep(filename,'.txt','.mat'),'file') ~= 0
    load(strrep(filename,'.txt','.mat'));
end

% Check to see if the file actually exists
if exist(filename,'file') == 0 && exist(strrep(filename,'.txt','.mat'),'file') == 0
    warning([filename ' does not exist']);
    l = []; nstep_conv = []; dpout = []; return;
end


%% Scan poisson and main loops from text file

% Scan through text file if not already saved to .mat file
if exist(strrep(filename,'.txt','.mat'),'file') == 0
    
    % Remove non characters from log file
    [~,~] = system(['tr -cd ''\11\12\15\40-\176'' < ' filename ' > ' strrep(filename,'.txt','.temp')]);
    movefile(strrep(filename,'.txt','.temp'),filename);

    % Determine number of blocks from first part of file
    [~,t] = system(['tail -n+1 ' filename ' | head -n200000']);
    bids = cellfun(@str2double,regexp(t,'BLOCK NO.\s+([\w\.\-\+]+)','tokens'))';
    nb = max(bids+1);

    % Determine number of steps and poisson steps from first part of file
    poisson_nstep = cellfun(@str2double,regexp(t,'poisson_nstep:\s+([\w\.\-\+]+)','tokens'))';
    dts = cellfun(@str2double,regexp(t,'dts:\s+([\w\.\-\+]+)','tokens'))';
    if isempty(dts)
        warning([filename ' is corrupted']);
        l = []; nstep_conv = []; dpout = []; return;
%         dts = 0;
    end
    if dts(1) == 0
        nstep = cellfun(@str2double,regexp(t,'nstep:\s+([\w\.\-\+]+)','tokens'))'; nstep = nstep(1);
    elseif dts(1) == 1
        nstep_cycle = cellfun(@str2double,regexp(t,'nstep_cycle:\s+([\w\.\-\+]+)','tokens'))';
        ncycle = cellfun(@str2double,regexp(t,'ncycle:\s+([\w\.\-\+]+)','tokens'))';
        nstep = nstep_cycle * ncycle;
    end
    
    % Determine turbulence model to size residual arrays
    ilos = cellfun(@str2double,regexp(t,'ilos:\s+([\w\.\-\+]+)','tokens'))';
%     davg_temp = regexp(log,'TOTAL DAVG\s+([\w\.\-\+]+)\s+([\w\.\-\+]+)\s+INLET','tokens');
%     if isempty(davg_temp)
%         ilos = 1;
%     end

    % Preallocate arrays
    if plot_poisson == 1
        l.poisson_nstep = nan(floor(poisson_nstep/50),1);
        l.poisson_davg_total = nan(floor(poisson_nstep/50),1);
        l.poisson_davg = nan(floor(poisson_nstep/50),nb);
    end
    if dts(1) == 0
        ni = floor(nstep/50);
    else
        ni = nstep;
    end
    l.nstep = nan(ni,1); l.nstep(1) = 0;
    if ilos(1) == 2
        l.davg_total = nan(ni,2);
    else
        l.davg_total = nan(ni,1);
    end
    l.m_in = nan(ni,1); l.m_out = nan(ni,1);
    l.Po_in = nan(ni,1); l.Po_out = nan(ni,1);
    l.To_in = nan(ni,1); l.To_out = nan(ni,1);
    l.Power = nan(ni,1); l.eta = nan(ni,1);
    if plot_blocks == 1
        if ilos(1) == 2
            l.davg = nan(ni,nb,2);
        else
            l.davg = nan(ni,nb,1);
        end
    end
    
    % Reduce file size of large unsteady logs by removing block residuals and bleed data
    if dts(1) == 1
        if plot_blocks == 0
            system(['sed ''/^BLOCK/ d'' < ' filename ' > log_temp.txt']);
            system('sed -i ''/^TOTAL MASS FLOW/ d'' log_temp.txt');
        else
            system(['sed ''/^TOTAL MASS FLOW/ d'' < ' filename ' >  log_temp.txt']);
        end
        current_name = 'log_temp.txt';
    else
        current_name = filename;
    end

    % Find lines that start each iteration
    [~,N] = system(['grep -n "STEP No" ' current_name ' | cut -d : -f 1']);
    N = textscan(N,'%d'); N = N{1};
    
    % Calculate number of lines in file
    [~,result] = system(['sed -n ''$='' ' current_name]);
    N_total = str2double(result); N_temp = 1e6;
    
    % Lines to start reading chunks from
    nn = ceil(N_total / N_temp); N_starts = 1:N_temp:nn*N_temp;
    for n = 2:length(N_starts); [~,i] = min(abs(N_starts(n) - N)); N_starts(n) = N(i); end;
    
    % Operate on N_temp lines at a time
    main_start = [0 0];
    for N_start = N_starts
        
        % Display line numbers and step number
        N_end = N_start + N_temp - 1;
        if plot_stuff == 1
            disp([N_start N_end l.nstep(find(isnan(l.nstep)==0,1,'last'))])
        end
        
        % Read lines into memory
        [~,t] = system(['tail -n+' num2str(N_start) ' ' current_name ' | head -n' num2str(N_temp) ...
            ' | grep -v "WRITING BLOCK"']);

        % Analyse and plot poisson results
        if plot_poisson == 1

            % Trim file to poisson loop results
            log_poisson = t;
            i = regexp(log_poisson,'INITIALISING MULTIGRID','end');
            log_poisson(1:i) = [];
            i = regexp(log_poisson,'INITIALISING TURBULENCE DAMPING FOR LAMINAR FLOW','start');
            log_poisson(i:end) = [];

            % Remove results for 49 out of 50 poisson iterations
            i = regexp(log_poisson,'POISSON STEP','start');
            q = zeros(size(log_poisson));
            for n = 1:50:length(i)-1
                q(i(n):(i(n+1)-1)) = 1;
            end
            log_poisson(q == 0) = [];

            % Read in step numbers
            l.poisson_nstep = cellfun(@str2double,regexp(log_poisson,...
                'POISSON STEP No\.\s+(\d+)','tokens'))';
            l.poisson_nstep = [0 ; l.poisson_nstep];

            % Read total residuals
            davg_temp = cellfun(@str2double,regexp(log_poisson,'TOTAL DAVG\s+([\w\.\-\+]+)','tokens'))';
            l.poisson_davg_total = [nan ; davg_temp];

            % Read block residuals
            davg_temp = cellfun(@str2double,regexp(log_poisson,'\d+ DAVG\s+([\w\.\-\+]+)','tokens'))';
            ns = length(l.poisson_nstep)-1;

            l.poisson_davg = reshape(davg_temp',[nb ns])';
            l.poisson_davg = [nan(1,nb) ; l.poisson_davg];    
        end

        % Find line number of end of last iteration
        N_end = find(cellfun(@isempty,regexp(regexp(t,'\n','split'),...
            'EFFICIENCY','match')) == 0,1,'last');

        % Calculate line numbers to operate on
        if isempty(N_end) == 1
            N_start = N_start + N_temp;
            N_end = N_start + N_temp;
            continue;
        else
            N_start = N_start + N_end;
            N_end = N_start + N_temp;
        end

        % Trim off top of file
        i = regexp(t,'INITIALISING SLIDING PLANES','end');
        t(1:i) = [];

        % Trim off bottom of file
        i = regexp(t,'EFFICIENCY = [\w\.\-]+','end');
        t(i(end)+1:end) = [];

        % Trim lines that indicating block writing
        i1 = regexp(t,'WRITING BLOCK NO. \d+\n','start');
        i2 = regexp(t,'WRITING BLOCK NO. \d+\n','end');
        q = zeros(size(t));
        for n = 1:length(i1)
            q(i1(n):i2(n)) = 1;
        end
        t(q == 1) = [];

        % Read in step numbers
        if dts(1) == 0
            temp = cellfun(@str2double,regexp(t,'STEP No\.\s+(\d+)','tokens'));
            i = 1 + temp / 50; i_davg = i;
        else
            temp = cellfun(@str2double,regexp(t,['EFFICIENCY =\s+[\w\.\-]+\s+' ...
            'OUTER STEP No\.\s+(\d+)'],'tokens'));
            i = temp;
        end

        if dts(1) == 1
            l.nstep(i+1) = temp;
            i = [i(1) i+1];
            i_davg = i;
        else
            l.nstep(i) = temp;
        end
        if main_start(1) == 0 && dts(1) == 0
            i = [1 i]; 
            main_start(1) = 1;
        end
        
        % Check if the file is corrupted by large jumps in step numbers
        if max(diff(i)) > 1
            warning([filename ' is corrupted']);
            l = []; nstep_conv = []; dpout = []; return;
        end

        % Read in inlet and outlet flow rates
        m_in = cellfun(@str2double,regexp(t,'INLET FLOW =\s+([\w\.\-]+)','tokens'))';
        i = i(1:length(m_in)); l.m_in(i) = m_in;
        l.m_out(i) = cellfun(@str2double,regexp(t,'OUTLET FLOW =\s+([\w\.\-]+)','tokens'))';

        % Read in stagnation pressures and temperatures
        l.Po_in(i) = cellfun(@str2double,regexp(t,'INLET STAG P =\s+([\w\.\-]+)','tokens'))';
        l.Po_out(i) = cellfun(@str2double,regexp(t,'OUTLET STAG P =\s+([\w\.\-]+)','tokens'))';
        l.To_in(i) = cellfun(@str2double,regexp(t,'INLET STAG T =\s+([\w\.\-]+)','tokens'))';
        l.To_out(i) = cellfun(@str2double,regexp(t,'OUTLET STAG T =\s+([\w\.\-]+)','tokens'))';

        % Power output and efficiency
        l.Power(i) = cellfun(@str2double,regexp(t,'POWER OUTPUT =\s+([\w\.\-]+)','tokens'))';
        l.eta(i) = cellfun(@str2double,regexp(t,'EFFICIENCY =\s+([\w\.\-]+)','tokens'))';

        % Total residuals
        if ilos(1) == 2
            if dts(1) == 0
                davg_temp = regexp(t,'TOTAL DAVG\s+([\w\.\-\+]+)\s+([\w\.\-\+]+)','tokens');
            else
                davg_temp = regexp(t,'TOTAL DAVG\s+([\w\.\-\+]+)\s+([\w\.\-\+]+)\s+INLET','tokens');
            end
            l.davg_total(i_davg(1:length(davg_temp)),:) = cellfun(@str2double,cat(1,davg_temp{:}));
        else
            l.davg_total(i_davg) = cellfun(@str2double,regexp(t,...
                'TOTAL DAVG\s+([\w\.\-\+]+)\s+INLET','tokens'))';
        end
        
        % Block residuals
        if plot_blocks == 1
            if ilos(1) == 2
                davg_temp = regexp(t,'\d+ DAVG\s+([\w\.\-\+]+)\s+([\w\.\-\+]+)','tokens');
                davg_temp = cellfun(@str2double,cat(1,davg_temp{:}));
            else
                davg_temp = cellfun(@str2double,regexp(t,'\d+ DAVG\s+([\w\.\-\+]+)','tokens'))';
            end
            ns = length(i_davg);
            nb = round(size(davg_temp,1)/ns);

            l.davg(i_davg,:,:) = permute(reshape(davg_temp',[ilos(1) nb ns]),[3 2 1]);
        end

        % Extract throttle boundary condition behaviour
        if isempty(regexp(t,'THROTTLE BOUNDARY CONDITION','match')) == 0;

            % Preallocate arrays
            if main_start(2) == 0
                l.m_target = nan(ni,1);
                l.m_current = nan(ni,1);
                l.dpout = nan(ni,1);
                main_start(2) = 1;
            end

            % Read in mass flows and pressure corrections
            l.m_target(i_davg) = cellfun(@str2double,regexp(t,'TARGET =\s+([\w\.\-]+)','tokens'))';
            l.m_current(i_davg) = cellfun(@str2double,regexp(t,'CURRENT =\s+([\w\.\-]+)','tokens'))';
            l.dpout(i_davg) = cellfun(@str2double,regexp(t,'CORRECTION =\s+([\w\.\-]+)','tokens'))';

        end
    end
end

% Store final pressure correction or return an empty array
if isfield(l,'dpout') == 1
    dpout = l.dpout(isnan(l.dpout) == 0); dpout = dpout(end);
else
    dpout = [];
end

%% Determine at which point convergence is reached or predict when it will occur

% Check if the solution has finished
dn = 10000; 
if isnan(l.m_in(end)) == 0 && l.nstep(end) > dn 

    % Convergence tolerance
    ddavg_tol = 0.1; nn = dn / 50; dmass_tol = 5e-5;

    % Calculate gradient of residuals trace
    n = l.nstep; 
%     davg = smooth(log10(l.davg_total(:,1)),nn);
    l.davg_total(1,:) = l.davg_total(2,:);
    davg = interp1(n(1:nn:end),log10(l.davg_total(1:nn:end,1)),n,'pchip','extrap');
    grad_davg = diff(davg) ./ diff(n); grad_davg = [grad_davg(1) ; grad_davg];
%     grad_davg = grad_mg(n,davg);

    % Find when gradient first falls below threshold
    nstep_conv = n(find(grad_davg > -ddavg_tol / dn & n > dn,1,'first'));

    % If unconverged extrapolate to predict it's occurance
    if isempty(nstep_conv) == 1

        % Extract last part of the run
        n_last = n(end-nn:end); grad_last = grad_davg(end-nn:end);

        % Fit line to the gradient
        [p,s,mu] = polyfit(n_last,grad_last,1); 

        % Extrapolate to tolerance value
        if abs(p(1)) > 1e-18 
            nstep_conv = round(interp1(polyval(p,n_last([1 end]),s,mu),n_last([1 end]),...
                -ddavg_tol / dn,'linear','extrap') / 50) * 50;
        end
    end

    % Check mass flow change over final portion of run
    dmass = abs(diff(l.m_out([end-nn end]))) / l.m_out(end);

    % If mass flow variation is too great predict convergence by polynomial fit
    if dmass > dmass_tol && 2*nn < length(n)

        % Fit quadratic
        [p,~,mu] = polyfit(n(end-2*nn:end),l.m_out(end-2*nn:end),2); 

        % Find step number where minimum will occur
        x_min = -p(2) / (2*p(1)); nstep_min = round((x_min * mu(2) + mu(1)) / 50) * 50;

        % Choose whichever convergence prediction is longer
        nstep_conv = max(nstep_conv,nstep_min);

    end
    
    % Limit convergence prediction to multiple of current number of steps
    nstep_conv = min(nstep_conv,2 * n(end));
else
    
    % Empty convergence prediction
    nstep_conv = length(l.nstep) * 50;
end


%% Plot the converenge history

% Plot poisson loop if required
if plot_poisson == 1
    figure(); hold on; grid on; box on;
    xlabel('Iteration'); ylabel('Poisson Residuals');
    for n = 1:size(l.poisson_davg,2)
        plot3(l.poisson_nstep,log10(l.poisson_davg(:,n)),ones(length(l.poisson_nstep),1)*(n-1),'k-')
    end
    axis([0 poisson_nstep (floor(log10(min(min(l.poisson_davg(:,:,1)))))) ...
        (ceil(log10(max(max(l.poisson_davg(:,:,1))))))])
end

% Set plot layout based upon throttle boundary condition presence
if isfield(l,'dpout') == 1
    n_plot = 3;
else
    n_plot = 2;
end

% Plot main loop convergence history
if plot_stuff == 1
    h = figure(); set(h,'Position',[1 82 1280 895]);
    
    % Total residuals
    subplot(n_plot,2,1); hold on; grid on; box on;
    semilogy(l.nstep,log10(l.davg_total(:,1)), 'k')
    if size(l.davg_total,2) == 2
        semilogy(l.nstep,log10(l.davg_total(:,2)), 'b')
    end
    ylabel('Total Residuals'); xlabel('Iteration');
%     axis([0 nstep (floor(log10(min(l.davg_total(:,1))))) ...
%         (ceil(log10(max(l.davg_total(:,1)))))])

    subplot(n_plot,2,2); hold on; grid on; box on;
    hold on
    plot(l.nstep,l.m_in, 'k')
    plot(l.nstep,l.m_out, 'b')
    ylabel('Mass Flow (kg/s)'); xlabel('Iteration'); legend('Inlet','Outlet','Location','South')
    
    subplot(n_plot,2,3); hold on; grid on; box on;
    hold on
    plot(l.nstep,l.Po_in, 'k')
    plot(l.nstep,l.Po_out, 'b')
    ylabel('Total Pressure (Pa)'); xlabel('Iteration'); legend('Inlet','Outlet','Location','South')
    
    subplot(n_plot,2,4); hold on; grid on; box on;
    hold on
    plot(l.nstep,l.To_in, 'k')
    plot(l.nstep,l.To_out, 'b')
%     plot(n,ft(n),'r')
    ylabel('Total Temperature (K)'); xlabel('Iteration'); legend('Inlet','Outlet','Location','South')    
    
    if n_plot == 3
        subplot(n_plot,2,5); hold on; grid on; box on;
        hold on
        plot(l.nstep,l.m_target,'r')
        plot(l.nstep,l.m_current,'k')
        ylabel('Outlet Mass Flow'); xlabel('Iteration');   
                
        subplot(n_plot,2,6); hold on; grid on; box on;
        hold on
        plot(l.nstep,l.dpout,'k')
        ylabel('Pressure Correction'); xlabel('Iteration');   
    end
end

% Plot individual block residuals if required
if plot_blocks == 1
    figure(); hold on; grid on; box on;
    xlabel('Iteration'); ylabel('Block Residuals');
    for n = 1:size(l.davg,2)
        plot3(l.nstep,log10(l.davg(:,n,1)),ones(length(l.nstep),1)*(n-1),'k-')
        if ilos(1) == 2
            plot3(l.nstep,log10(l.davg(:,n,2)),ones(length(l.nstep),1)*(n-1),'b-')
        end
    end
    axis([0 nstep (floor(log10(min(min(l.davg(:,:,1)))))) ...
        (ceil(log10(max(max(l.davg(:,:,1))))))])
end

% Save history file
save(strrep(filename,'.txt','.mat'),'l');

% Delete text file
if exist(filename,'file') ~= 0; delete(filename); end;


end


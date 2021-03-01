function is_stalling = ts_check_stall(job)
% TS_CHECK_STALL  Check if a calculation is stalling by negative mass flow, or falling flow in log file
%
%   is_stalling = TS_CHECK_STALL(job)
%
%   job - struct containing job.directory and job.logname, or string of logname
%   is_stalling - 0 or 1 depending if calculation is stalling

% Define mass flow tolerance as a fraction of the final value
% dm_tol = 6e-5;
dm_tol = 1e-3;

% Get convergence history
if isstruct(job) == 1
    l = ts_plot_conv([job.directory job.logname],0);
else
    l = ts_plot_conv(job,0);
end

% Check for NaNs at end of run
if isnan(l.m_in(end)) == 0

    % Portion of run to check
    n_end = find(isnan(l.m_in) == 0,1,'last');
    nn = 10000 / 50;

    % Change in mass flow rate over the final part of the calculation
    p = polyfit(linspace(0,1,nn)',l.m_in(n_end-nn+1:n_end),1);
    dm = diff(polyval(p,[0 1]));

    % Check mass flow variation
    if dm / l.m_in(n_end) < - dm_tol
        is_stalling = 1;
%         dm / l.m_in(n_end)
%         figure(); hold on; grid on; box on;
%         plot(linspace(0,1,nn),l.m_in(n_end-nn+1:n_end),'-')
%         plot(linspace(0,1,nn),polyval(p,linspace(0,1,nn)));
    elseif l.m_in(n_end) < 0
        is_stalling = 1;
    else
        is_stalling = 0;
    end
    
else
    is_stalling = 1;
end

end
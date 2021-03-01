function job = ts_run_shaft(job,dpout,Wx_tol,plot_stuff,use_avg,bpid_out)
% TS_RUN_SHAFT  Iterate on backpressure until shaft power consumed is zero
%
%   job = TS_RUN_SHAFT(job,dpout,shaft_tol,plot_stuff,use_avg)
%
%   job - data structure containing name of previous run
%   dpout - initial change in backpressure to probe the power gradient
%   shaft_tol - power in KW to converge calculation to
%   plot_stuff - optional 0 or 1 to show working
%   use_avg - optional 0 or 1 to use average TS solution
%   bpid_out - optional vector with bid and pid of outlet patch

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Default to not use the average
if exist('use_avg','var') == 0
    use_avg = 0;
end

% Open figure window if requried
if plot_stuff == 1
    figure(); hold on; grid on; box on;
    xlabel('Back Pressure / Pa'); ylabel('Shaft Power / KW');
end

% Number of cases to use in fit
n_fit = 4;

% Initialise back pressure and power vectors
pouts = []; Wxs = [];

% Average suffix
if use_avg == 1; avg_suff = '_avg'; else; avg_suff = ''; end

% Read in previous solution
g = ts_read_hdf5([job.directory job.name avg_suff '.hdf5']);

% Get outlet block and patch ids of the actual blocks to change pressure on if required
if exist('bpid_out','var') == 0
    for bid = ts_get_bids(g)
        if exist('bid_out','var') == 0
            for pid = ts_get_pids(g,bid)
                if g{bid+1}.patch{pid+1}.attribute.kind == 1
                    bid_out = bid; pid_out = pid;
                end
            end
        end
    end
else
    bid_out = bpid_out(1); pid_out = bpid_out(2); 
end

% Extract backpressure and power
Wx = calc_power(g);
pout = g{bid_out+1}.patch{pid_out+1}.pout;
pouts = [pouts ; pout]; Wxs = [Wxs ; Wx];

% Plot initial power requirements
l = plot(pouts,Wxs,'.-'); drawnow; pause(1);

% Perturb back pressure by specified value
pout_new = pouts + dpout;

% Continue to run TURBOSTREAM solutions until shaft power required is zero
while abs(Wxs(end)) > Wx_tol

    % Apply new back pressure
    g{bid_out+1}.patch{pid_out+1}.pout = pout_new;
    
    % Turn off poisson loop
    g{1}.av.poisson_nstep = 0; g{1}.av.poisson_restart = 1;
    
    % Resubmit with new backpressure
    job = ts_smonitor(ts_sbatch(job,g)); 
    
    % Read TURBOSTREAM solution
    g = ts_read_hdf5([job.directory job.name avg_suff '.hdf5']);
    
    % Add backpressure and power to vectors
    Wx = calc_power(g);
    pout = g{bid_out+1}.patch{pid_out+1}.pout;
    pouts = [pouts ; pout]; Wxs = [Wxs ; Wx];

    % Update plot
    if plot_stuff == 1
        set(l,'xdata',pouts,'ydata',Wxs); drawnow; pause(1);
    end
    
    % Interpolate new backpressure
%     if max(Wxs) > 0 && min(Wxs) < 0
%         pout_new = interp1(Wxs,pouts,0,'pchip');
%     else
%         pout_new = interp1(Wxs,pouts,0,'linear','extrap');
%     end
    
    % Find next backpressure by linear fit
    [~,i] = sort(abs(Wxs)); if length(i) > n_fit; i = i(1:n_fit); end
    pout_new = polyval(polyfit(Wxs(i),pouts(i),1),0);
    
end

end

function Wx = calc_power(g)

% Initialise shaft power
Wx = 0; cp = g{1}.av.cp;

% Add all outlet power fluxes
for bid = ts_get_bids(g)
    for pid = ts_get_pids(g,bid)
        if g{bid+1}.patch{pid+1}.attribute.kind == 1 
            outlet = ts_structured_cut(g,bid,'en','en',1,'en',1,'en');
            [To_out,m_out] = ts_mass_average(outlet,'To',3); m_out = m_out * double(outlet.bv.nblade);
            Wx = Wx + cp * To_out * m_out;
        elseif g{bid+1}.patch{pid+1}.attribute.kind == 10
            bleed = ts_structured_cut(g,bid,1,'en',1,1,1,'en');
            [To_bl,m_bl] = ts_mass_average(bleed,'To',3); m_bl = m_bl * double(bleed.bv.nblade);
            Wx = Wx + cp * To_bl * m_bl;
        end
    end
end

% Subtract all inlet power fluxes
for bid = ts_get_bids(g)
    for pid = ts_get_pids(g,bid)
        if g{bid+1}.patch{pid+1}.attribute.kind == 0
            inlet = ts_structured_cut(g,bid,1,1,1,'en',1,'en');
            [To_in,m_in] = ts_mass_average(inlet,'To',3); m_in = m_in * double(inlet.bv.nblade);
            Wx = Wx - cp * To_in * m_in;
        end
    end
end

% Return power in kW
Wx = Wx / 1000;


end
function c = exp_chic_cut(c,plot_stuff)
% EXP_CHIC_CUT  Truncate a set of characteristics to unstalled range

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Truncate characteristics to negative slope
phi_max = 0.63; dpsi = 0.06;

% Initialise point acceptance vector
i_delete = zeros(length(c.phi),1);

% Delete all points up to first crossing of max flow
i = find(c.phi(1:end-1) > phi_max - 0.005 & c.phi(1:end-1) < phi_max + 0.005 &...
    c.phi(2:end) < c.phi(1:end-1));
i_delete(1:min(i)) = 1;
i_start = [];
for n = 1:length(i)-1
    if i(n+1) - i(n) > 5
        i_start = [i_start ; i(n+1)];
    end
end

% Delete all points after stall point
i = find(c.psi(1:end-1) - c.psi(2:end) > dpsi & c.psi(1:end-1) > 0.2);
di = diff(i); di = [100 ; di];
i(di < 5) = [];
i_end = i;

% Create deletion vector
for n = 1:length(i_start)
    i_delete(i_end(n):i_start(n)) = 1;
end
i_delete(i_end(end):end) = 1;

% Plot start and end points of each characteristic to show working
if plot_stuff == 1
    figure(); hold on; grid on; box on;
    plot(c.phi,c.psi,'k-')
    plot(c.phi(i_delete == 1),c.psi(i_delete == 1),'g-')
    plot(c.phi(i_start),c.psi(i_start),'ro','markerfacecolor','r')
    plot(c.phi(i_end),c.psi(i_end),'bo','markerfacecolor','b')
end

% Delete out of range points
varnames = fieldnames(c);
for v = 1:length(varnames)
    c.(varnames{v})(i_delete == 1) = [];
end

% Put rows of nans in for clear plotting
i = find(diff(c.phi) > 0.1);
i = [0 ; i ; length(c.phi)];
for v = 1:length(varnames)
    c_new.(varnames{v}) = [];
    for n = 1:length(i)-1
        c_new.(varnames{v}) = [c_new.(varnames{v}) ; nan ; c.(varnames{v})(i(n)+1:i(n+1))];
    end
end
c = c_new;

end
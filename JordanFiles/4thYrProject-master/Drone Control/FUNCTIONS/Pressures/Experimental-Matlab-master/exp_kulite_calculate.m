function g = exp_kulite_calculate(e,c,x,r)
% Process a row of over tip kulite data

% Get number of kulites and operating points
nk = length(e.N.V_kulite); no = size(e.V,3);

% Calculate gradient of kulite calibrations
p = zeros(2,nk); P = squeeze(c.P);
for k = 1:nk
    
    % Average voltages in time
    V = squeeze(mean(c.V(:,k,:),1));

    % Calculate linear fit
    p(:,k) = polyfit(V,P,1);
end

% Calculate offsets for current run
p(2,:) = -mean(e.V_zero(:,e.N.V_kulite),1) .* p(1,:);

% Apply calibrations to kulite data
P = zeros(size(e.V,1),nk,no);
for k = 1:nk
    P(:,k,:) = e.V(:,e.N.V_kulite(k),:) * p(1,k) + p(2,k);
end

% Find falling edges of once per rev signals
V_trig = 0.5 * (max(e.V(:,e.N.V_shaft,1)) + min(e.V(:,e.N.V_shaft,1)));
a_log = e.V(2:end,e.N.V_shaft,:) < V_trig & e.V(1:end-1,e.N.V_shaft,:) > V_trig;

% Ensemble average kulite data
rate = 100e3; n_ring = 50; f = 3500/60;
nt = round(5 * rate / (n_ring * f));
P_av = zeros(nt,nk,no); w = zeros(1,nk,no);
for k = 1:nk
    for o = 1:no                  
        % Find indices of falling edges
        a_ind = 1:length(a_log(:,1,o)); a_ind = a_ind(a_log(:,1,o));

        % Remove spurious triggerings of once per rev signals
        dt = diff(a_ind); a_ind(dt < 0.5 * max(dt)) = [];

        % Average first nt points
        for m = 1:length(a_ind)-1
            P_av(:,k,o) = P_av(:,k,o) + P(a_ind(m):a_ind(m)+nt-1,k,o);
        end
        P_av(:,k,o) = P_av(:,k,o)/(length(a_ind)-1);
    end
end

% Calculate disk angular velocity
w = rate * 2 * pi / mean(diff(a_ind));

% Calculate rotor relative coordinates
t_rel = (1:nt)' * w ./ rate;
rt_rel = repmat(t_rel * r,[1 nk]);

% Repeat axial coordinates
x = repmat(x,[nt 1]);

% Record output variables
g.rt_rel = rt_rel;
g.x = x;
g.P_av = P_av;

end
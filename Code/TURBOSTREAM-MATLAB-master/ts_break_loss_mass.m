function [Yp,Zeta,c] = ts_break_loss_mass(inlet,outlet,mass_frac,is_exp)
% Break un the loss in a wake into different components

if exist('mass_frac','var') == 0
    mass_frac = 0.2;
end

% Determine flow properties to use
if exist('is_exp','var') == 0
    is_exp = 1;
end

% Get mass flow rates at inlet and exit
[r_1,m_1] = ts_mass_average(inlet,'r',2);
[r_2,m_2] = ts_mass_average(outlet,'r',2);

m_1_total = sum(m_1); m_2_total = sum(m_2);

% Calculate nearby indices
cum_m = [0 cumsum(m_1)];
[~,j_hub] = min(abs(cum_m - mass_frac*m_1_total));
[~,j_cas] = min(abs(cum_m - (1-mass_frac)*m_1_total));

% Calculate radiuses at desired mass fractions
r_1_hub = interp1(cum_m(j_hub-3:j_hub+1), r_1(j_hub-3:j_hub+1), mass_frac*m_1_total,'pchip');
r_1_cas = interp1(cum_m(j_cas-3:j_cas+1), r_1(j_cas-3:j_cas+1), (1-mass_frac)*m_1_total,'pchip');

r_2_hub = interp1([0 cumsum(m_2)], r_2, mass_frac*m_2_total,'pchip');
r_2_cas = interp1([0 cumsum(m_2)], r_2, (1-mass_frac)*m_2_total,'pchip');

% Interpolate new meshes at hub, mid and casing sections
varnames = {'x','rt','ro','rovx','rovr','rorvt','roe'};
if is_exp == 1
    varnames = [varnames ,'Cpo','Cp'];
end

for v = 1:length(varnames)
    if isequal(size(inlet.(varnames{v})),size(inlet.r)) == 1
        % Interpolate inlet
        var_new = zeros(size(inlet.r) + [2 0]);
        r_new = zeros(size(inlet.r) + [2 0]);
        for k = 1:size(inlet.r,2)
            r_new(:,k) = sort([inlet.r(:,k) ; r_1_hub ; r_1_cas]);
            var_new(:,k) = interp1(inlet.r(:,k),inlet.(varnames{v})(:,k),r_new(:,k),'pchip');
        end
        inlet.(varnames{v}) = var_new;
    end
end
inlet.r = r_new;
for v = 1:length(varnames)
    if isequal(size(outlet.(varnames{v})),size(outlet.r)) == 1
        % Interpolate outlet
        var_new = zeros(size(outlet.r) + [2 0]);
        r_new = zeros(size(outlet.r) + [2 0]);
        for k = 1:size(outlet.r,2)
            r_new(:,k) = sort([outlet.r(:,k) ; r_2_hub ; r_2_cas]);
            var_new(:,k) = interp1(outlet.r(:,k),outlet.(varnames{v})(:,k),r_new(:,k),'pchip');
        end
        outlet.(varnames{v}) = var_new;       
    end
end
outlet.r = r_new;

% Break up inlet and outlet into 3 sections
[~,j_1_hub] = min(abs(inlet.r(:,1) - r_1_hub));
[~,j_1_cas] = min(abs(inlet.r(:,1) - r_1_cas));
[~,j_2_hub] = min(abs(outlet.r(:,1) - r_2_hub));
[~,j_2_cas] = min(abs(outlet.r(:,1) - r_2_cas));

inlet_hub = inlet; inlet_mid = inlet; inlet_cas = inlet;
outlet_hub = outlet; outlet_mid = outlet; outlet_cas = outlet;

varnames = [varnames ,'r'];

for v = 1:length(varnames)
    inlet_hub.(varnames{v}) = inlet_hub.(varnames{v})(1:j_1_hub,:);
    inlet_mid.(varnames{v}) = inlet_mid.(varnames{v})(j_1_hub:j_1_cas,:);
    inlet_cas.(varnames{v}) = inlet_cas.(varnames{v})(j_1_cas:end,:);
    outlet_hub.(varnames{v}) = outlet_hub.(varnames{v})(1:j_2_hub,:);
    outlet_mid.(varnames{v}) = outlet_mid.(varnames{v})(j_2_hub:j_2_cas,:);
    outlet_cas.(varnames{v}) = outlet_cas.(varnames{v})(j_2_cas:end,:);    
end

inlet_hub = ts_secondary(inlet_hub); 
inlet_mid = ts_secondary(inlet_mid); 
inlet_cas = ts_secondary(inlet_cas);
outlet_hub = ts_secondary(outlet_hub); 
outlet_mid = ts_secondary(outlet_mid); 
outlet_cas = ts_secondary(outlet_cas);

% Determine flow parameters
if is_exp == 1
    Po_var = 'Cpo';
    P_var = 'Cp';
else
    if outlet.bv.rpm == 0
        Po_var = 'Po';
        P_var = 'P';
    else
        Po_var = 'Po_rel';
        P_var = 'P';
    end
end

% Calculate pressure loss coefficient for each region
Po_1 = ts_mass_average(inlet,Po_var,3);
P_1 = ts_mass_average(inlet,P_var,3);
Po_2 = ts_mass_average(outlet,Po_var,3);

Po_1_hub = ts_mass_average(inlet_hub,Po_var,3);
Po_2_hub = ts_mass_average(outlet_hub,Po_var,3);

Po_1_mid = ts_mass_average(inlet_mid,Po_var,3);
Po_2_mid = ts_mass_average(outlet_mid,Po_var,3);

Po_1_cas = ts_mass_average(inlet_cas,Po_var,3);
Po_2_cas = ts_mass_average(outlet_cas,Po_var,3);

Yp.hub = (Po_1_hub - Po_2_hub) / (Po_1 - P_1);
Yp.mid = (Po_1_mid - Po_2_mid) / (Po_1 - P_1);
Yp.cas = (Po_1_cas - Po_2_cas) / (Po_1 - P_1);

Yp.total = (Po_1 - Po_2) / (Po_1 - P_1);

% Calculate entropy loss coefficient for each region
V_1 = ts_mass_average(inlet,'V',3);
s_1 = ts_mass_average(inlet,'s',3);
s_2 = ts_mass_average(outlet,'s',3);
T_2 = ts_mass_average(outlet,'T',3);

s_1_hub = ts_mass_average(inlet_hub,'s',3);
s_2_hub = ts_mass_average(outlet_hub,'s',3);
T_2_hub = ts_mass_average(outlet_hub,'T',3);

s_1_mid = ts_mass_average(inlet_mid,'s',3);
s_2_mid = ts_mass_average(outlet_mid,'s',3);
T_2_mid = ts_mass_average(outlet_mid,'T',3);

s_1_cas = ts_mass_average(inlet_cas,'s',3);
s_2_cas = ts_mass_average(outlet_cas,'s',3);
T_2_cas = ts_mass_average(outlet_cas,'T',3);

Zeta.hub = T_2_hub * (s_2_hub - s_1_hub) / (0.5 * V_1^2);
Zeta.mid = T_2_mid * (s_2_mid - s_1_mid) / (0.5 * V_1^2);
Zeta.cas = T_2_cas * (s_2_cas - s_1_cas) / (0.5 * V_1^2);

Zeta.total = T_2 * (s_2 - s_1) / (0.5 * V_1^2);

% Record cuts
c.inlet.hub = inlet_hub; c.inlet.mid = inlet_mid; c.inlet.cas = inlet_cas;
c.outlet.hub = outlet_hub; c.outlet.mid = outlet_mid; c.outlet.cas = outlet_cas;

end
function c = ts_calc_chic(inlet,outlet,rpm)
% TS_CALC_CHIC  Calculate performance parameters between inlet and outlet cut planes

% Use relative blade speed for counter rotating machines
if inlet.bv.rpm ~= 0 && outlet.bv.rpm ~= 0
    rpm = abs(inlet.bv.rpm) + abs(outlet.bv.rpm);
end

% Midspan blade speed
r_mean = ts_area_average(inlet,'r',3);
U_mean = r_mean * rpm * 2 * pi / 60;

% Inlet density and axial velocity for low speed performance
ro_in = ts_area_average(inlet,'ro',3);
Vx_in = ts_area_average(inlet,'Vx',3);
Vx_mid = interp1(ts_mass_average(inlet,'r_nondim',2),ts_area_average(inlet,'Vx',2),0.5);

% Average pressures and temperatures
[Po_in, m_in, A_in] = ts_mass_average(inlet,'Po',3);
[Po_out, m_out, A_out] = ts_mass_average(outlet,'Po',3);
To_in = ts_mass_average(inlet,'To',3);
To_out = ts_mass_average(outlet,'To',3);
P_out = ts_area_average(outlet,'P',3);

% Gas constants
cp = inlet.av.cp; ga = inlet.av.ga;

% Flow rate measures
c.phi = Vx_in / U_mean;
c.mrtp_in = m_in * (cp * To_in)^0.5 / (A_in * Po_in);
c.mrtp_out = m_out * (cp * To_out)^0.5 / (A_out * Po_out);
c.mass = m_in * double(inlet.bv.nblade);

% Pressure rise measures
c.P_dyn = (0.5 * ro_in * U_mean^2);
c.Cpts = (P_out - Po_in) / c.P_dyn;
c.Cptt = (Po_out - Po_in) / c.P_dyn;
c.Rpo = Po_out / Po_in;
c.psi = cp * (To_out - To_in) / U_mean^2;

% Performance measures
if To_out > To_in
    c.eta_poly = ((ga-1) / ga) * (log(Po_out / Po_in) / log(To_out / To_in));
    To_out_is = To_in * (Po_out / Po_in) ^ ((ga-1)/ga);
    c.eta_is = (To_out_is - To_in) / (To_out - To_in);
else
    c.eta_poly = (ga / (ga-1)) * (log(To_out / To_in) / log(Po_out / Po_in));
end


end


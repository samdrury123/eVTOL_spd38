% Thrust_FoM_Calculate

% Calculate thrust components from hdf5 TS output for eVTOL rig
% Note: two methods either using the numerical pressure along cuts to get
% force or using the SFME equation 

% Megan McEveley

% % Clear the workspace
% close all; clear; clc;

% Load the eVTOL variables
eVTOL;
A =  0.0357;
rho = 1.2; 

% Make array to plot
diffuser_thrust = [];
inlet_thrust = [];
fan_thrust = [];
thrust = [];
efficiency = [];
power = [];
merit = [];
ideal_power = [];


%% Make control volume cuts
for i = 1:length(duct_names)
    
    g = duct_names{i};    

    % Check duct and block ids are correct
    ts_check_rpms(g);

%% Outlet CV

    % Outlet CV cuts for baseline outlet
    clear cv;
    
    cv.inlet = ts_structured_cut(g,22,1,1,1,'en',1,1);
    cv.diffuser = ts_structured_cut(g,22,1,'en','en','en',1,1);
    cv.cone = ts_structured_cut(g,22,1,'en',1,1,1,1);
    cv.outlet = ts_structured_cut(g,22,'en','en',1,'en',1,1);
    cv.tail = ts_structured_cut(g,24,1,1,1,20,1,1);

    tailnames = {'tail'};
    inletnames = {'inlet'};
    diffusernames = {'diffuser'};
    outletnames = {'outlet'};
    conenames = {'cone'};
    cutnames = fieldnames(cv);

    % Plot cuts on the meridional view
    for n = 1:length(cutnames)
        plot(cv.(cutnames{n}).x,cv.(cutnames{n}).r,'bo-');

    end
    
    % Cone pressure forces
    for n = 1:length(conenames)

        cut_sec = ts_secondary(cv.(conenames{n}));
        P = cut_sec.P;
        r = cut_sec.r;
        Pav = 0.5*(P(1:end-1)+P(2:end))-Pref;
        rav = 0.5*(r(1:end-1)+r(2:end));
        dr = r(2:end)-r(1:end-1);
        Ax = abs(2*pi*(rav.*dr));
        
        A_cone = sum(Ax(:)) ;

        forces.cone{n} = (Pav.*Ax);
        
    end
    
    
    %Tail pressure forces
    for n = 1:length(tailnames)

        cut_sec = ts_secondary(cv.(tailnames{n}));
        P = cut_sec.P;
        r = cut_sec.r;
        Pav = 0.5*(P(1:end-1)+P(2:end))-Pref;
        rav = 0.5*(r(1:end-1)+r(2:end));
        dr = r(2:end)-r(1:end-1);
        Ax = abs(2*pi*(rav.*dr));

        forces.tail{n} = (Pav.*Ax);

    end
    

    % Diffuser pressure forces
    for n = 1:length(diffusernames)

        cut_sec = ts_secondary(cv.(diffusernames{n}));
        P = cut_sec.P;
        r = cut_sec.r;
        Pav = 0.5*(P(1:end-1)+P(2:end))-Pref;
        rav = 0.5*(r(1:end-1)+r(2:end));
        dr = r(2:end)-r(1:end-1);
        Ax = abs(2*pi*(rav.*dr));
        
        A_diffuser = sum(Ax(:));
        
        forces.diffuser{n} = (Pav.*Ax);

    end


    % Inlet pressure & momentum forces
    for n = 1:length(inletnames)

        cut_sec = ts_secondary(cv.(inletnames{n}));
        P = cut_sec.P;
        r = cut_sec.r;
        rovx = cut_sec.rovx;
        Vx = cut_sec.Vx;

        Pav = 0.5*(P(1:end-1)+P(2:end));
        rav = 0.5*(r(1:end-1)+r(2:end));
        dr = r(2:end)-r(1:end-1);
        Ax = abs(2*pi*(rav.*dr));
        
        rovx_av = 0.5*(rovx(1:end-1)+rovx(2:end));
        vx_av = 0.5*(Vx(1:end-1)+Vx(2:end));
        mdot_av = rovx_av.*Ax;
        mass_flow = sum(mdot_av(:));
         
        P3 = sum(Pav(:))/length(Pav);
        A_in = sum(Ax(:));
        
        forces.inlet.pressure{n} = (Pav.*Ax);
        forces.inlet.momentum{n} = (mdot_av.*vx_av);

    end

    % Outlet pressure & momentum forces
    for n = 1:length(outletnames)

        cut_sec = ts_secondary(cv.(outletnames{n}));
        P = cut_sec.P;
        r = cut_sec.r;
        Vx = cut_sec.Vx;
        rovx = cut_sec.rovx;

        Pav = 0.5*(P(1:end-1)+P(2:end));
        rovx_av = 0.5*(rovx(1:end-1)+rovx(2:end));
        rav = 0.5*(r(1:end-1)+r(2:end));
        dr = r(2:end)-r(1:end-1);
        Ax = abs(2*pi*(rav.*dr));
        vx_av = 0.5*(Vx(1:end-1)+Vx(2:end));
        mdot_av = rovx_av.*Ax;
        A_out = sum(Ax(:));
      
        forces.outlet.pressure{n} = (Pref.*Ax);
        forces.outlet.momentum{n} = (mdot_av.*vx_av);
        Vx4 = vx_av;
        r4 =rav;
        alpha4 = cut_sec.Alpha;
    end

  
    cone_force = 0;
    for n = 1:length(conenames)
        cone_force = cone_force + sum(forces.cone{n}(:));
    end
    
       tail_force = 0;
    for n = 1:length(tailnames)
       tail_force = tail_force + sum(forces.tail{n}(:));
    end

    diffuser_force = 0;
    for n = 1:length(diffusernames)
        diffuser_force = diffuser_force + sum(forces.diffuser{n}(:));
    end

    outlet_force = 0;
    for n = 1:length(outletnames)
        outlet_force = outlet_force + sum(forces.outlet.momentum{n}(:)) + sum(forces.outlet.pressure{n}(:));
    end

    inlet_force = 0;
    for n = 1:length(inletnames)
        inlet_force = inlet_force + sum(forces.inlet.pressure{n}(:)) + sum(forces.inlet.momentum{n}(:));
    end
    
    diffuser_thrust_force(i) = diffuser_force + cone_force;
    diffuser_thrust_SFME(i) = outlet_force - inlet_force;
    
    if i == 1
        diffuser_thrust(i) = 0;  
    else 
        diffuser_thrust(i) = diffuser_thrust_force(i);    
    end 
    
 
    %% Inlet CV
     
    % Intake cuts
    clear cv;
    cv.intake_1 = ts_structured_cut(g,19,1,'en','en','en',1,1);
    cv.intake_2 = ts_structured_cut(g,23,1,'en','en','en',1,1);
    cv.outlet_1 = ts_structured_cut(g,23,'en','en',1,'en',1,1);
    cv.cone_1 = ts_structured_cut(g,23,1,'en',1,1,1,1);
    cv.cone_2 = ts_structured_cut(g,21,'en','en',1,'en',1,1);


    % Control volume names
    intakenames = {'intake_1' 'intake_2'}; 
    outletnames = {'outlet_1'};
    conenames = {'cone_1' 'cone_2'};
    cutnames = fieldnames(cv);

    % Plot cuts on the meridional view
    for n = 1:length(cutnames)
        plot(cv.(cutnames{n}).x,cv.(cutnames{n}).r,'bo-');
    end


    
    % Cone pressure forces
    for n = 1:length(conenames)

        cut_sec = ts_secondary(cv.(conenames{n}));
        P = cut_sec.P;
        r = cut_sec.r;
        Pav = 0.5*(P(1:end-1)+P(2:end))-Pref;
        rav = 0.5*(r(1:end-1)+r(2:end));
        dr = r(2:end)-r(1:end-1);
        Ax = abs(2*pi*(rav.*dr));
        forces.cone{n} = (Pav.*Ax);

    end

    % Intake pressure forces
    for n = 1:length(intakenames)

        cut_sec = ts_secondary(cv.(intakenames{n}));
        P = cut_sec.P;
        r = cut_sec.r;
        Pav = 0.5*(P(1:end-1)+P(2:end))-Pref;
        rav = 0.5*(r(1:end-1)+r(2:end));
        dr = r(2:end)-r(1:end-1);
        Ax = abs(2*pi*(rav.*dr));
        forces.intake{n} = (Pav.*Ax);

    end
    
    % Outlet pressure & momentum forces
    for n = 1:length(outletnames)

        cut_sec = ts_secondary(cv.(outletnames{n}));
        P = cut_sec.P;
        r = cut_sec.r;
        rovx = cut_sec.rovx;
        Vx = cut_sec.Vx;

        Pav = 0.5*(P(1:end-1)+P(2:end))-Pref;
        rav = 0.5*(r(1:end-1)+r(2:end));
        dr = r(2:end)-r(1:end-1);
        Ax = 2*pi*(rav.*dr);
        rovx_av = 0.5*(rovx(1:end-1)+rovx(2:end));
        vx_av = 0.5*(Vx(1:end-1)+Vx(2:end));
        mdot_av = rovx_av.*Ax;
        mass_flow = sum(mdot_av(:))
        
        forces.outlet.pressure{n} = (Pav.*Ax);
        forces.outlet.momentum{n} = abs(mdot_av.*vx_av);

    end

    cone_force = sum(forces.cone{1}(:)) + sum(forces.cone{2}(:)) ;
    intake_force = sum(forces.intake{1}(:)) + sum(forces.intake{2}(:));
    out_force = sum(forces.outlet.momentum{1}(:)) + sum(forces.outlet.pressure{1}(:));
    
    
    thrust_inlet_SFME(i) = out_force;
    thrust_inlet_force(i) = -(cone_force + intake_force);
    inlet_thrust(i) = thrust_inlet_force(i);
    %% Fan CV

    % Fan CV cuts
    clear cv;
    cv.inlet = ts_structured_cut(g,23,'en','en',1,'en',1,1);
    cv.outlet = ts_structured_cut(g,22,1,1,1,'en',1,1);

    % Plot cuts on the meridional view
    cutnames = fieldnames(cv);
    for n = 1:length(cutnames)
        plot(cv.(cutnames{n}).x,cv.(cutnames{n}).r,'go-');
    end

    % Inlet pressure & momentum forces
    cut_sec = ts_secondary(cv.inlet);
    P = cut_sec.P;
    r = cut_sec.r;
    rovx = cut_sec.rovx;
    Vx = cut_sec.Vx;
    To1 = cut_sec.To;
    Po1 = cut_sec.Po;

    To1_av = 0.5*(To1(1:end-1)+To1(2:end));
    Po1_av = 0.5*(Po1(1:end-1)+Po1(2:end));
    Pav = 0.5*(P(1:end-1)+P(2:end));
    rav = 0.5*(r(1:end-1)+r(2:end));
    dr = r(2:end)-r(1:end-1);
    Ax = abs(2*pi*(rav.*dr));
    rovx_av = 0.5*(rovx(1:end-1)+rovx(2:end));
    vx_av = 0.5*(Vx(1:end-1)+Vx(2:end));
    mdot_av = rovx_av.*Ax;

    forces.inlet.pressure{1} = abs(Pav.*Ax);
    forces.inlet.momentum{1} = abs(mdot_av.*vx_av);
    
    % Mass average flow properties at inlet
    To1 = sum(mdot_av .* To1_av) / sum(mdot_av);
    Po1 = sum(mdot_av .* Po1_av) / sum(mdot_av);

    % Outlet pressure & momentum forces
    cut_sec = ts_secondary(cv.outlet);
    P = cut_sec.P;
    r = cut_sec.r;
    rovx = cut_sec.rovx;
    Vx = cut_sec.Vx;
    To2 = cut_sec.To;
    Po2 = cut_sec.Po;

    To2_av = 0.5*(To2(1:end-1)+To2(2:end));
    Po2_av = 0.5*(Po2(1:end-1)+Po2(2:end));
    Pav = 0.5*(P(1:end-1)+P(2:end));
    rav = 0.5*(r(1:end-1)+r(2:end));
    dr = r(2:end)-r(1:end-1);
    Ax = 2*pi*(rav.*dr);
    rovx_av = 0.5*(rovx(1:end-1)+rovx(2:end));
    vx_av = 0.5*(Vx(1:end-1)+Vx(2:end));
    mdot_av = rovx_av.*Ax;
    
    forces.outlet.pressure{1} = abs(Pav.*Ax);
    forces.outlet.momentum{1} = abs(mdot_av.*vx_av);
    
   
    Vx3 = vx_av;
    r3 = rav;
    alpha3 = cut_sec.Alpha;
    
    % Mass average flow properties at outlet
    To2 = sum(mdot_av .* To2_av) / sum(mdot_av);
    Po2 = sum(mdot_av .* Po2_av) / sum(mdot_av);
    vx3 = sum(mdot_av .* vx_av) / sum(mdot_av);

    % Isentropic Temperature
    To2s = To1*(Po2./Po1).^(0.4/1.4);
    
    % Total forces
    outlet_force = sum(forces.outlet.momentum{1}(:)) + sum(forces.outlet.pressure{1}(:));
    inlet_force = sum(forces.inlet.pressure{1}(:)) + sum(forces.inlet.momentum{1}(:));
    
    % Fan thrust and efficiency
    fan_thrust(i) = outlet_force - inlet_force ;
    efficiency(i) = (To2s-To1)./(To2-To1);
    velocity(i) = vx3;
  
    
    %% Power calculation

    % Power from average temperature
    power(i) = sum(mdot_av(:))*cp*(To2-To1);
    efficiency_long = efficiency 
    ideal_power(i) = power(i)*efficiency(i);
    
    %% Figure of Merit calculation

    thrust(i) = inlet_thrust(i) + fan_thrust(i) + diffuser_thrust(i);
    merit(i) = (efficiency(1)/efficiency(i))*(thrust(i)/power(i))*((thrust(i)/(2*rho*A)).^0.5);
    
end

%% Thrust plots
% % Data from control volumes
'Diffuser Data'
efficiency
fan_thrust
inlet_thrust
diffuser_thrust
thrust
merit
power


% % Non-dimensionalised thrust
Vx =  30;
inlet_thrust = inlet_thrust/(A*rho*0.5*Vx*Vx);
diffuser_thrust = diffuser_thrust/(A*rho*0.5*Vx*Vx);
fan_thrust = fan_thrust/(A*rho*0.5*Vx*Vx);
thrustND = thrust/(A*rho*0.5*Vx*Vx);

   

   

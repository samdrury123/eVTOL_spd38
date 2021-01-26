% Instructions

% 1. Log current imbalance at part speed and appropriate filter settings
% 2. Add trial mass to specified position
% 3. MATLAB will calculate masses to be added to correct original imbalance from phase and magnitude.
% 4. Re-run rig at the same speed and check corrected imbalance

% Connect to PXIE
pxie.chan = [N.V_shaft N.V_accel];
s = exp_pxie_initialise(pxie);

% Log initial imbalance
uiwait(msgbox(['Start rig and set RPM = ' num2str(rpm)]));
[V_1,ph_1] = exp_imbalance_log(s,N,rpm);

% Determine whether to run trial mass test
if exist('r_calib','var') == 0 || exist('ph_trial','var') == 0

    % Log imbalance with trial mass
    uiwait(msgbox(['Add trial mass of ' num2str(mass_trial) ' g to ' num2str(th_trial) ...
        ' deg position and start rig']));
    [V_2,ph_2] = exp_imbalance_log(s,N,rpm); 

    % Calculate imbalance vectors
    xy_1 = [V_1*cosd(ph_1) V_1*sind(ph_1)];
    xy_2 = [V_2*cosd(ph_2) V_2*sind(ph_2)];
    xy_trial = xy_2 - xy_1;
    ph_trial = rad2deg(atan2(xy_trial(2),xy_trial(1)));
    V_trial = sum(xy_trial.^2).^0.5;

    % Find scale factor between logged voltages and mass using trial value
    r_calib = mass_trial / V_trial;

    % Display factors
    disp(['r_calib = ' num2str(r_calib) '  ph_trial = ' num2str(ph_trial)]);
    plot_trial = 1;
    
else
    
    % Don't plot trial mass results
    plot_trial = 0;
    
end

% Scale all voltages into mass values
V_1_new = V_1 * r_calib; 
if plot_trial == 1; V_2_new = V_2 * r_calib; V_trial_new = V_trial * r_calib; end

% Correct phase values using trial value
ph_1_new = ph_1 - ph_trial + th_trial; 
if plot_trial == 1; ph_2_new = ph_2 - ph_trial + th_trial; ph_trial_new = th_trial; end

% Plot imbalance vectors on phasor plot
if plot_trial == 1
    h = figure(); polar(0,max([V_2_new V_1_new]),'k-');
    hold on; grid on; box on;
    l(1) = polar(deg2rad([0 ph_1_new]),[0 V_1_new]);
    l(2) = polar(deg2rad([0 ph_2_new]),[0 V_2_new]);
    l(3) = polar(deg2rad([ph_1_new ph_2_new]),[V_1_new V_2_new],'k-');
    polar(deg2rad([0 ph_trial_new]),[0 V_trial_new],'k-');
else
    h = figure(); polar(0,V_1_new,'k-');
    hold on; grid on; box on;
    l(1) = polar(deg2rad([0 ph_1_new]),[0 V_1_new]);
end

% Flip original log to find imbalance value
xy_imbal = -[V_1_new*cosd(ph_1_new) V_1_new*sind(ph_1_new)];

% Find the angles between the imbalance and the allowable positions
n_bal = length(th_bal); th_bal = reshape(th_bal,[n_bal 1]);
xy_bal = [cosd(th_bal) sind(th_bal)];
a_bal = acosd(sum(xy_bal .* repmat(xy_imbal / (sum(xy_imbal.^2)^0.5) ,[n_bal 1]),2));

% Choose two positions in whose sector the imbalance lies within
[~,i] = sort(abs(a_bal)); i = i(1:2); a_bal = a_bal(i); th_bal = th_bal(i); xy_bal = xy_bal(i,:);

% Solve triangle for masses to add to two balance position
mass_bal = xy_bal' \ xy_imbal';

% Calculate coordinates to plot balance masses
xy_bal = [mass_bal .* cosd(th_bal) mass_bal .* sind(th_bal)];
xy_temp = sum(xy_bal,1);

% Print required mass instructions
for n = 1:length(mass_bal)
    disp(['Add mass of ' num2str(mass_bal(n)) ' g to ' num2str(th_bal(n)) ' deg position'])
end

% Add balance masses to the phasor plot
figure(h);
l(4) = polar(deg2rad([0 th_bal(1)]),[0 mass_bal(1)]);
l(5) = polar([deg2rad(th_bal(1)) atan2(xy_temp(2),xy_temp(1))],[mass_bal(1) sum(xy_temp.^2).^0.5]);

% Add legend to the phasor plot
if plot_trial == 1
    legend(l,'Imbalance','Imbalance + Trial','Trial','Mass 1','Mass 2','location','best');
else
    l(2:3) = []; legend(l,'Imbalance','Mass 1','Mass 2','location','best');
end

% Remove the trial mass
if plot_trial == 1; uiwait(msgbox('Remove trial mass')); end

% Continue to log imbalance using the same scale factors
while 1

    % Relog the imbalance value
    uiwait(msgbox('Add balance masses and start rig'));
    [V_3,ph_3] = exp_imbalance_log(s,N,rpm); 

    % Plot the final imbalance value
    V_3_new = V_3 * r_calib; ph_3_new = ph_3 - ph_trial + th_trial; 
    figure(h); l(6) = polar(deg2rad([0 ph_3_new]),[0 V_3_new]);
    
end

% Close pxie connection
% delete(s);

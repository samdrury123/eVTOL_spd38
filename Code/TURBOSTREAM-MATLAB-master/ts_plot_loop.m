function [F H] = ts_plot_loop(directory, filestring, coords,...
    plot_entropy, plot_wakes, plot_2D, plot_break, plot_curv, des, M_plot, H, cuts, normal)
% Plot the loss loop for a blade row defined by all files found by matching
% the filestring

% Set default behaviour
if exist('plot_entropy','var') == 0
    plot_entropy = 1;
end

if exist('plot_wakes','var') == 0
    plot_wakes = 0;
end

if exist('plot_2D','var') == 0 || plot_2D == 0
    plot_2D = -1;
end

if exist('plot_break','var') == 0
    plot_break = 0;
end

if exist('plot_curv','var') == 0
    plot_curv = 0;
end

if exist('H','var') == 0 || isempty(H) == 1
    H.loop = figure(); H.dev = figure(); H.wakes = figure(); H.mid = figure(); H.Cp_mid = figure();
    H.break = figure(); H.break_2 = figure; H.break_wakes = figure(); H.IM = figure();
else
%     H.loop = figures(1); H.wakes = figures(2); H.mid = figures(3); H.break = figures(4);...
%         H.break_2 = figures(5); H.break_wakes = figures(6); H.IM = figures(6); H.dev = figures(7);
end

if exist('cuts','var') == 0
    write_cuts = 0;
else
    write_cuts = 1;
end

%plot_entropy_mid = 0;
plot_entropy_2D = plot_entropy;

%% Read in files
% Get files to define the loops

[~, A] = system(['ls ' directory filestring]);
B = regexp(A,'\n','split'); 
A = cell2struct(B,'name'); A(end) = [];
B = regexp(B,'/','split');

for n = 1:length(B)-1
    A(n).name = B{n}{end};
end

if isempty(regexp(filestring,'?','match')) == 0
    filestring = regexp(filestring,'?','split');
    filestring = [filestring{1} '*'];
end

list.number = zeros(length(A),1);
list.name = cell(length(A),1);
for n = 1:length(A)
    list.name{n} = strrep(A(n).name,'.mat','.hdf5');
%     B = regexp(list.name{n},'\d+','match');
%     list.number(n) = str2double(B(end-1));
    B = regexp(list.name{n},'_','split');
    if strcmp(list.name{n}(end),'5') == 1
        list.number(n) = str2double(B{end-1});
    else
        list.number(n) = str2double(B{end});
    end
end

% Remove files not in the loops
% Failure cases
list.name(list.number < 1000) = [];
list.number(list.number < 1000) = [];

% Paraview files
par = cellfun(@isempty,regexp(list.name,'paraview','once'));
list.number(par == 0) = [];
list.name(par == 0) = [];

% Split into separate loops
list.M = round(list.number/100)/100;
[~, I J] = unique(list.M);
F = cell(length(I),1); G = cell(length(I),1); M_loop = zeros(length(I),1);
for m = 1:length(I)
    F{m}.M = list.M(I(m));
    M_loop(m) = list.M(I(m));
    F{m}.name = list.name(J == m);
    F{m}.number = list.number(J == m);
end

[~, m_des] = min(abs(M_loop - des.M));

for m = 1:length(F)
    % Initialise and read in cuts
    F{m}.Inlet = cell(length(F{m}.name),1);
    F{m}.Outlet = cell(length(F{m}.name),1);
    F{m}.Blade = cell(length(F{m}.name),1);
    F{m}.pout = zeros(length(F{m}.name),1);
    F{m}.m_exit = zeros(length(F{m}.name),1);
    F{m}.mf = zeros(length(F{m}.name),1);
    F{m}.mr = zeros(length(F{m}.name),1);
    G{m} = cell(length(F{m}.name),1);
    if write_cuts == 1
        F{m}.cuts = cell(length(F{m}.name),size(cuts,1));
    end
    
    if exist([directory strrep(filestring, '*', [num2str(F{m}.M*100) '00_cuts.mat'])],'file') == 0
        for n = 1:length(F{m}.name)
            filename = F{m}.name{n};
            if exist([directory strrep(filename, '.hdf5', '.mat')],'file') == 0
                g = ts_read_hdf5([directory filename]);
                save([directory strrep(filename, '.hdf5', '.mat')],'g','-v7.3')
            else
                load([directory strrep(filename, '.hdf5', '.mat')]);
            end
                
%             G{m}{n} = g;
            
            if exist('coords','var') == 0
                % If coords undefined use inlet and exit planes
                for bid = ts_get_bids(g)
                    for pid = ts_get_pids(g,bid)
                        if g{bid+1}.patch{pid+1}.attribute.kind == 0
                            bid1 = bid;
                            i1 = g{bid+1}.patch{pid+1}.attribute.ist;
                        end
                        if g{bid+1}.patch{pid+1}.attribute.kind == 1
                            bid2 = bid;
                            i2 = g{bid+1}.patch{pid+1}.attribute.ist;
                        end   
                    end
                end

                F{m}.Inlet{n} = ts_structured_cut(g,bid1,i1+1,i1+1,1,'en',1,'en');
                F{m}.Outlet{n} = ts_structured_cut(g,bid2,i2+1,i2+1,1,'en',1,'en');     

            elseif isstruct(coords) == 0
                % Use structured cuts if coords is a numerical array
                bid1 = coords(1,1); i1 = coords(1,2);
                bid2 = coords(2,1); i2 = coords(2,2);
                F{m}.Inlet{n} = ts_structured_cut(g,bid1,i1+1,i1+1,1,'en',1,'en');
                F{m}.Outlet{n} = ts_structured_cut(g,bid2,i2+1,i2+1,1,'en',1,'en');  
            else
                % Use unstructured cuts as defined in the data structure coords
                xr = coords.in.xr; bids = coords.in.bids;
                F{m}.Inlet{n} = ts_structunstruct_cut_2(g,xr,bids,2,1);
                xr = coords.out.xr; bids = coords.out.bids;
                F{m}.Outlet{n} = ts_structunstruct_cut_2(g,xr,bids,2,1);
            end

            % Determine back pressure and exit mass flow rate
            for bid = ts_get_bids(g)
                for pid = ts_get_pids(g,bid)
                    if g{bid+1}.patch{pid+1}.attribute.kind == 1
                        F{m}.pout(n) = g{bid+1}.patch{pid+1}.pout;
                        Exit = ts_secondary(ts_structured_cut(g,bid,'en','en',1,'en',1,'en'));
                        [~, mass] = ts_mass_average(Exit,'ro',3);
                        F{m}.m_exit(n) = mass * double(Exit.bv.nblade);
                    end
                end
            end
            
            % Determine downstream blade incidence if present
            [~, bid_os chi] = ts_find_chord(g,[min(min(F{m}.Inlet{n}.r)) max(max(F{m}.Inlet{n}.r))]);
            if length(bid_os) > 1
                if isfield(F{m},'rpm_2') == 0
                    F{m}.rpm_2 = zeros(length(F{m}.name),1);
                end
                F{m}.rpm_2(n) = g{bid_os(end)+1}.bv.rpm;
                F{m}.chi_2 = chi(1,end);
            end
            
            % Write out blade surface cut
            F{m}.Blade{n} = ts_structured_cut(g,bid_os(1),1,'en',1,'en',1,1);
           
            % Determine mass flow rates in inter platform gaps if present
            if size(coords,1) == 3
                Front = ts_secondary(ts_structured_cut(g,coords(3,1),1,1,1,'en',1,'en'));
                Rear = ts_secondary(ts_structured_cut(g,coords(3,2),1,1,1,'en',1,'en'));

                [~, F{m}.mf(n)] = ts_mass_average(Front,'ro',3);
                [~, F{m}.mr(n)] = ts_mass_average(Rear,'ro',3);
            end  
            
            % Write out other cuts as required
            for c = 1:size(cuts,1)
                F{m}.cuts{n,c} = ts_structured_cut(g,cuts(c,1),cuts(c,2),cuts(c,2),1,'en',1,'en');
            end
        end

        % Determine inlet and outlet angle
        [~, bid_os chi] = ts_find_chord(g,[min(min(F{m}.Inlet{n}.r)) max(max(F{m}.Inlet{n}.r))]);

        o = 1;
        for bid = bid_os
            if min(min(min(g{bid+1}.x))) > min(min(F{m}.Inlet{n}.x)) &&...
                    max(max(max(g{bid+1}.x))) < max(max(F{m}.Outlet{n}.x))
                F{m}.chi_1 = chi(1,o);
                F{m}.chi_1_out = chi(2,o);
            end
            o = o+1;
        end
        
        f = F{m};
        save([directory strrep(filestring, '*', [num2str(F{m}.M*100) '00_cuts.mat'])],'f','-v7.3')
    else
        load([directory strrep(filestring, '*', [num2str(F{m}.M*100) '00_cuts.mat'])])
        F{m} = f;
    end
end

% Set the colours for the loops
for m = 1:length(I)
    if length(I) > 1
        F{m}.colour = [(m-1)/(length(I)-1) 0 0];
    else
        F{m}.colour = [0 0 0];
    end
end

%% Write paraview hdf5 file if necessary
if plot_curv == 1 || plot_2D ~= -1 || plot_break ~= 0
    for m = 1:length(F)
        for n = 1:length(F{m}.number)
            pv_directory = strrep(strrep(directory,'TURBOSTREAM','PARAVIEW'),'TS','PV');
            delete([pv_directory 'x_intersect.*'])          
            pv_filename = strrep(F{m}.name{n},'.hdf5','_paraview.hdf5');
            if exist([pv_directory pv_filename],'file') == 0
                load([directory strrep(F{m}.name{n},'.hdf5','.mat')])
                if size(coords,1) == 0
                if exist('normal','var') ~= 0
                    ts_export_paraview(g,[pv_directory pv_filename],'Track',[0:6 coords(3,:)],normal);
                else
                    ts_export_paraview(g,[pv_directory pv_filename],'Track',[0:6 coords(3,:)]);
                end
                else
                if exist('normal','var') ~= 0
                    ts_export_paraview(g,[pv_directory pv_filename],'Track',[0:6],normal);
                else
                    ts_export_paraview(g,[pv_directory pv_filename],'Track',[0:6]);
                end
                end
            end
        end
    end
end
        

%% Plot total loss loops
% Calculate and plot the total loss loop
M = cell(size(F)); p = zeros(size(F)); q = zeros(size(F));
for m = 1:length(F)
    Inlet = ts_secondary(F{m}.Inlet);
    Outlet = ts_secondary(F{m}.Outlet);
    Alpha_mid = zeros(size(Inlet)); M_mid = zeros(size(Inlet)); Loss = zeros(size(Inlet));
    F{m}.Alpha_mid_out = zeros(size(Inlet));
    F{m}.V_mid_2 = zeros(length(F{m}.name),1); F{m}.Alpha_mid_2 = zeros(length(F{m}.name),1);
    
    for n = 1:length(Inlet)
        if plot_entropy == 0
            % Yp parameters
            if Inlet{n}.bv.rpm == 0
                Po1 = ts_mass_average(Inlet{n}, 'Po', 3);
                Po2 = ts_mass_average(Outlet{n}, 'Po', 3);
            else
                Po1 = ts_mass_average(Inlet{n}, 'Po_rel', 3);
                Po2 = ts_mass_average(Outlet{n}, 'Po_rel', 3); 
            end
            P1 = ts_mass_average(Inlet{n}, 'P', 3);
            Loss(n) = (Po1 - Po2) ./ (Po1 - P1);
        else
            % Zeta parameters
            V1 = ts_mass_average(Inlet{n}, 'V_rel', 3);
            s2 = ts_mass_average(Outlet{n}, 's', 3);
            s1 = ts_mass_average(Inlet{n}, 's', 3);
            T2 = ts_mass_average(Outlet{n}, 'T', 3);
            Loss(n) = 2 * (s2 - s1) .* T2 ./ (V1 .* V1);
        end

        % Mach number and incidence mass averaged over mid region
        mid = 0.1;
        r_in_nondim = (Inlet{n}.r(:,1) - min(min(Inlet{n}.r))) /...
            (max(max(Inlet{n}.r)) - min(min(Inlet{n}.r)));
        J_mid_in = find(r_in_nondim>(0.5-mid/2) & r_in_nondim<(0.5+mid/2));
        if Inlet{n}.bv.rpm == 0
            [Alpha1 mass_r_in] = ts_mass_average(Inlet{n}, 'Alpha', 2);
            M1 = ts_mass_average(Inlet{n}, 'M', 2);
        else
            [Alpha1 mass_r_in] = ts_mass_average(Inlet{n}, 'Alpha_rel', 2);
            M1 = ts_mass_average(Inlet{n}, 'M_rel', 2);
        end

        Alpha_mid(n) = sum(Alpha1(J_mid_in).*(mass_r_in(J_mid_in))) / sum(mass_r_in(J_mid_in));
        M_mid(n) = sum(M1(J_mid_in).*(mass_r_in(J_mid_in))) / sum(mass_r_in(J_mid_in));
        
        % Get exit flow angle
        r_out_nondim = (Outlet{n}.r(:,1) - min(min(Outlet{n}.r))) /...
            (max(max(Outlet{n}.r)) - min(min(Outlet{n}.r)));
        J_mid_out = find(r_out_nondim>(0.5-mid/2) & r_out_nondim<(0.5+mid/2));
        if Inlet{n}.bv.rpm == 0
            [Alpha2 mass_r_out] = ts_mass_average(Outlet{n}, 'Alpha', 2);
        else
            [Alpha2 mass_r_out] = ts_mass_average(Outlet{n}, 'Alpha_rel', 2);
        end        
        F{m}.Alpha_mid_2(n) = sum(Alpha2(J_mid_out).*(mass_r_out(J_mid_out))) ...
            / sum(mass_r_out(J_mid_out));
        
        % Get downstream absolute exit angle
%         [Alpha2 mass_r_out] = ts_mass_average(Outlet{n}, 'Alpha', 2);
        V2 = ts_mass_average(Outlet{n}, 'V', 2);
%         F{m}.Alpha_mid_2(n) = sum(Alpha2(J_mid_out).*(mass_r_out(J_mid_out))) ...
%             / sum(mass_r_out(J_mid_out));
        F{m}.V_mid_2(n) = sum(V2(J_mid_out).*(mass_r_out(J_mid_out))) / sum(mass_r_out(J_mid_out));
    end

    if exist('M_plot','var') == 0 || isempty(M_plot) == 1
        M_plot(m) = F{m}.M;
    end
    
    % Sort all variables in increasing yaw angle
    if Inlet{1}.bv.rpm ~= 0
        [F{m}.Alpha_mid n] = sort(Alpha_mid,'descend');
    else
        [F{m}.Alpha_mid n] = sort(Alpha_mid,'ascend');
    end
    
    % Record in data structures
    F{m}.M_mid = M_mid(n);
    F{m}.Loss = Loss(n);
    F{m}.Inlet = F{m}.Inlet(n);
    F{m}.Outlet = F{m}.Outlet(n);
    F{m}.pout = F{m}.pout(n);
    F{m}.number = F{m}.number(n);
    F{m}.name = F{m}.name(n);
    F{m}.Alpha_mid_2 = F{m}.Alpha_mid_2(n);

    if isempty(G{m}{1}) == 0
        G{m} = G{m}{n};
    end
    
    if isfield(F{m},'rpm_2') == 1
        F{m}.rpm_2 = F{m}.rpm_2(n);
        F{m}.V_mid_2 = F{m}.V_mid_2(n);
    end

    % Plot loss loops against incidence
    if Inlet{1}.bv.rpm ~= 0
        F{m}.Inc = - (F{m}.Alpha_mid - F{m}.chi_1);
    else
        F{m}.Inc = F{m}.Alpha_mid - F{m}.chi_1;
    end
    
    % Calculate deviation
    if Inlet{1}.bv.rpm ~= 0
        F{m}.Dev = - (F{m}.Alpha_mid_2 - F{m}.chi_1_out);
    else
        F{m}.Dev = F{m}.Alpha_mid_2 - F{m}.chi_1_out;
    end    
    
    % Calculate downstream incidence
    if isfield(F{m},'rpm_2') == 1
        Vm_2 = F{m}.V_mid_2 .* cosd(F{m}.Alpha_mid_2);
        r_mid_2 = 0.5*(max(max(Outlet{1}.r)) + min(min(Outlet{1}.r)));
        Vt_rel_2 = F{m}.V_mid_2 .* sind(F{m}.Alpha_mid_2) - (2*pi*F{m}.rpm_2 / 60) * r_mid_2;
        F{m}.Inc_2 = - (atand(Vt_rel_2 ./ Vm_2) - F{m}.chi_2);
    end

    figure(H.loop)
    hold on; view(2);
    p(m) = plot3(F{m}.Inc,F{m}.Loss,F{m}.number,'k.-','LineWidth',2,'MarkerSize',20);
    set(p(m),'Color',F{m}.colour)
    M{m} = ['M_{in} = ' num2str(F{m}.M)];
    
    figure(H.dev)
    hold on; view(2);
    q(m) = plot3(F{m}.Inc,F{m}.Dev,F{m}.number,'k.-','LineWidth',2,'MarkerSize',20);
    set(q(m),'Color',F{m}.colour)
end

% Add a legend
figure(H.loop)
legend(p,M);
title([filestring(1:2) ' Total Loss Loops'])
xlabel('\alpha_{1,mid} - \chi_{1,mid}','FontSize',16);

if plot_entropy == 1
    ylabel('\zeta_{s}','FontSize',16);
else
    ylabel('Y_{p}','FontSize',16);
end

[~, n_min] = min(abs(F{m_des}.Alpha_mid - des.Alpha));
Loss_min = F{m_des}.Loss(n_min); axis auto
V = axis; V(3) = 0; V(4) = 3*min(Loss_min); axis(V);
plot(V(1:2),[min(Loss_min) min(Loss_min)],'k-')
plot(V(1:2),[min(Loss_min)*1.5 min(Loss_min)*1.5],'k-')

% Deviation plot
figure(H.dev)
legend(q,M);
title([filestring(1:2) ' Deviation Loops'])
xlabel('\alpha_{1,mid} - \chi_{1,mid}','FontSize',16);
ylabel('\alpha_{2.mid} - \chi_{2,mid}','FontSize',16);


for m = 1:length(F)
    % Determine failure incidences
%     [Loss_min n_min] = min(F{m}.Loss);
    [~, n_min] = min(F{m}.Loss);

    Loss_1 = F{m}.Loss/Loss_min; Loss_2 = F{m}.Loss/Loss_min;
    Loss_1(n_min:end) = 0; Loss_2(1:n_min) = 0;
    
    % Remove values outside either ranges
    Fail_Alpha_1 = F{m}.Alpha_mid(Loss_1 ~= 0); Fail_pout_1 = F{m}.pout(Loss_1 ~= 0); 
    Fail_Loss_1 = Loss_1(Loss_1 ~= 0);
    Fail_Alpha_2 = F{m}.Alpha_mid(Loss_2 ~= 0); Fail_pout_2 = F{m}.pout(Loss_2 ~= 0); 
    Fail_Loss_2 = Loss_2(Loss_2 ~= 0);
    
    [Fail_Alpha_1 n] = sort(Fail_Alpha_1); Fail_pout_1 = Fail_pout_1(n); Fail_Loss_1 = Fail_Loss_1(n);
    [Fail_Alpha_2 n] = sort(Fail_Alpha_2); Fail_pout_2 = Fail_pout_2(n); Fail_Loss_2 = Fail_Loss_2(n);

    % Interpolate 1.5x limit
    if isempty(Fail_Loss_1) == 0 && max(Fail_Loss_1) > 1.25 && sum(isnan(Fail_Loss_1) == 1) == 0 &&...
            min(Fail_Loss_1) < 1.5
        [Loss n] = unique(Fail_Loss_1); Alpha = Fail_Alpha_1(n); pout = Fail_pout_1(n);
        F{m}.Alpha_fail(1) = interp1(Loss,Alpha,1.5,'pchip','extrap');
        F{m}.pout_fail(1) = interp1(Loss,pout,1.5,'pchip','extrap');
        Loss_fail(1) = min(Loss_min)*1.5;
    else
        [F{m}.Alpha_fail(1) n] = min(F{m}.Alpha_mid);
        F{m}.pout_fail(1) = F{m}.pout(n);
        Loss_fail(1) = F{m}.Loss(n);
    end
    if isempty(Fail_Loss_2) == 0 && max(Fail_Loss_2) > 1.25 && sum(isnan(Fail_Loss_2) == 1) == 0 &&...
            min(Fail_Loss_2) < 1.5
        [Loss n] = unique(Fail_Loss_2); Alpha = Fail_Alpha_2(n); pout = Fail_pout_2(n);
        F{m}.Alpha_fail(2) = interp1(Loss,Alpha,1.5,'pchip','extrap');
        F{m}.pout_fail(2) = interp1(Loss,pout,1.5,'pchip','extrap');
        Loss_fail(2) = min(Loss_min)*1.5;
    else
        [F{m}.Alpha_fail(2) n] = max(F{m}.Alpha_mid);
        F{m}.pout_fail(2) = F{m}.pout(n);
        Loss_fail(2) = F{m}.Loss(n);
    end

    if Inlet{1}.bv.rpm ~= 0
        F{m}.Inc_fail = - (F{m}.Alpha_fail + F{m}.chi_1);
    else
        F{m}.Inc_fail = F{m}.Alpha_fail - F{m}.chi_1;
    end

    % Plot on failure incidences
%     plot([F{m}.Inc(n_min) F{m}.Inc(n_min)], [0 min(Loss_min)],'k-')
    figure(H.loop)
    plot(F{m}.Inc_fail(1), Loss_fail(1), 'o','Color',F{m}.colour,'MarkerSize',10)
    plot(F{m}.Inc_fail(2), Loss_fail(2), 'o','Color',F{m}.colour,'MarkerSize',10)

    % Plot downstream wakes at 1, 1.5 x min loss
    if plot_wakes == 1
        figure(H.wakes)
        title([filestring(1:2) ' Wakes at High Incidence Failure'])
        % Determine cuts to plot
        [~, n_1] = min(abs(Loss_1-1.5));
        [~, n_2] = min(abs(Loss_2-1.5));

        % Plot wake at minimum loss
        hfig = subplot(length(F),3,m*3-1);
        ts_plot_wake(F{m}.Inlet{n_min},F{m}.Outlet{n_min},hfig,plot_entropy,1,1);
        xlabel(num2str(F{m}.Inc(n_min)))
        % Plot wake at -ve 1.5 x min loss
        hfig = subplot(length(F),3,m*3-2);
        ts_plot_wake(F{m}.Inlet{n_1},F{m}.Outlet{n_1},hfig,plot_entropy,1,1);
        xlabel(num2str(F{m}.Inc(n_1)))
        % Plot wake at +ve 1.5 x min loss
        hfig = subplot(length(F),3,m*3);
        ts_plot_wake(F{m}.Inlet{n_2},F{m}.Outlet{n_2},hfig,plot_entropy,1,1);    
        xlabel(num2str(F{m}.Inc(n_2)))

    elseif plot_wakes == 2
        % Wakes of all points
        figure(H.wakes)
        for n = 1:length(Inlet)
            hfig = subplot(4,5,n);
            ts_plot_wake(F{m}.Inlet{n},F{m}.Outlet{n},hfig,plot_entropy,1,1)
        end

    else
        if ishandle(H.wakes)
            close(H.wakes)
        end
    end
    
    % Calculate incidence at which maximum curvature is observed in the
    % total loss loop
%     dAlpha = 0.1;
%     % Fit both sides separately with a high order polynomial, differentiate
%     % twice and find the maximum curvature
%     P_1 = polyfit(Fail_Alpha_1,Fail_Loss_1,6);
%     Curv_Alpha_1 = min(Fail_Alpha_1):dAlpha:max(Fail_Alpha_1);
%     Curv_Loss_1 = polyval(P_1,Curv_Alpha_1)*Loss_min;
%     d2L_dA2_1 = (Curv_Loss_1(3:end) - 2*Curv_Loss_1(2:end-1) + Curv_Loss_1(1:end-2)) / dAlpha^2;
%     
%     P_2 = polyfit(Fail_Alpha_2,Fail_Loss_2,6);
%     Curv_Alpha_2 = min(Fail_Alpha_2):dAlpha:max(Fail_Alpha_2);
%     Curv_Loss_2 = polyval(P_2,Curv_Alpha_2)*Loss_min;
%     d2L_dA2_2 = (Curv_Loss_2(3:end) - 2*Curv_Loss_2(2:end-1) + Curv_Loss_2(1:end-2)) / dAlpha^2;
    
    % Calculate the gradient between each point, fit the result and
    % find the maximum change.
    dL_dA_1 = (Fail_Loss_1(2:end) - Fail_Loss_1(1:end-1))./(Fail_Alpha_1(2:end) - Fail_Alpha_1(1:end-1));
    dL_dA_2 = (Fail_Loss_2(2:end) - Fail_Loss_2(1:end-1))./(Fail_Alpha_2(2:end) - Fail_Alpha_2(1:end-1));

    Curv_Alpha_1 = (Fail_Alpha_1(2:end) + Fail_Alpha_1(1:end-1)) / 2;
    Curv_Alpha_2 = (Fail_Alpha_2(2:end) + Fail_Alpha_2(1:end-1)) / 2;
    
%     figure()
%     subplot(3,1,1)
%     hold on
%     plot(Fail_Alpha_1,Fail_Loss_1,'k.-');
%     plot(Fail_Alpha_2,Fail_Loss_2,'k.-');    
%     v1 = axis;
%     
%     subplot(3,1,2)
%     hold on
%     plot(Curv_Alpha_1,dL_dA_1,'k.-');
%     plot(Curv_Alpha_2,dL_dA_2,'k.-');
%     v = axis;
%     axis([v1(1:2) v(3:4)]);
    
    d2L_dA2_1 = (dL_dA_1(2:end) - dL_dA_1(1:end-1))./(Curv_Alpha_1(2:end) - Curv_Alpha_1(1:end-1));
    d2L_dA2_2 = (dL_dA_2(2:end) - dL_dA_2(1:end-1))./(Curv_Alpha_2(2:end) - Curv_Alpha_2(1:end-1));
    
    Curv_Alpha_1 = (Curv_Alpha_1(2:end) + Curv_Alpha_1(1:end-1)) / 2;
    Curv_Alpha_2 = (Curv_Alpha_2(2:end) + Curv_Alpha_2(1:end-1)) / 2;    

%     subplot(3,1,3)
%     hold on
%     plot(Curv_Alpha_1,d2L_dA2_1,'k.-');
%     plot(Curv_Alpha_2,d2L_dA2_2,'k.-');
%     v = axis;
%     axis([v1(1:2) v(3:4)]);

    
    % Assign into output structure
    if isempty(d2L_dA2_1) == 0 && max(Fail_Loss_1) > 1.25 && sum(isnan(Fail_Loss_1)) == 0
        [~, n_1] = max(d2L_dA2_1);
        F{m}.Alpha_curv(1) = Curv_Alpha_1(n_1);
        [Loss n] = unique(Fail_Loss_1); Alpha = Fail_Alpha_1(n); pout = Fail_pout_1(n);
        F{m}.pout_curv(1) = interp1(Alpha,pout,Curv_Alpha_1(n_1),'pchip');
        Loss_curv(1) = interp1(Alpha,Loss,Curv_Alpha_1(n_1),'pchip') * Loss_min;
    else
        [F{m}.Alpha_curv(1) n] = min(F{m}.Alpha_mid);
        F{m}.pout_curv(1) = F{m}.pout(n);        
        Loss_curv(1) = F{m}.Loss(n);
    end

    if isempty(d2L_dA2_2) == 0 && max(Fail_Loss_2) > 1.25 && sum(isnan(Fail_Loss_2)) == 0
        [~, n_2] = max(d2L_dA2_2);
        F{m}.Alpha_curv(2) = Curv_Alpha_2(n_2);
        [Loss n] = unique(Fail_Loss_2); Alpha = Fail_Alpha_2(n); pout = Fail_pout_2(n);
        F{m}.pout_curv(2) = interp1(Alpha,pout,Curv_Alpha_2(n_2),'pchip');
        Loss_curv(2) = interp1(Alpha,Loss,Curv_Alpha_2(n_2),'pchip') * Loss_min;
    else
        [F{m}.Alpha_curv(2) n] = max(F{m}.Alpha_mid);
        F{m}.pout_curv(2) = F{m}.pout(n);      
        Loss_curv(2) = F{m}.Loss(n);
    end
    
    % Calculate incidence
    if Inlet{1}.bv.rpm ~= 0
        F{m}.Inc_curv = - (F{m}.Alpha_curv + F{m}.chi_1);
    else
        F{m}.Inc_curv = F{m}.Alpha_curv - F{m}.chi_1;
    end
    
    % Plot on loss loops
    figure(H.loop)
    plot(F{m}.Inc_curv(1), Loss_curv(1),'o','Color',F{m}.colour,'MarkerSize',10)
    plot(F{m}.Inc_curv(2), Loss_curv(2),'o','Color',F{m}.colour,'MarkerSize',10)
end

%% Plot midheight loss loops
if plot_2D ~= -1
    M = cell(length(F)*2,1); p = zeros(length(F)*2,1);
        
    % Define MISES path and variables
    ext = 'mises';
    path = strrep(strrep(directory,'TURBOSTREAM','MISES'),'TS','Mises');
    path = [path 'Turbulent/N' num2str(F{1}.Inlet{1}.bv.nblade) '/R' num2str(plot_2D) '/'];
    ga = 1.4; R = 287; cp = 1005;
    
    % Plot a mid height slice of loss to compare with 2D flow solvers
    for m = 1:length(F)
        Inlet = ts_secondary(F{m}.Inlet); 
        Outlet = ts_secondary(F{m}.Outlet);
        Blade = ts_secondary(F{m}.Blade);
        Loss_mid = zeros(size(Inlet)); 
        
        r_nondim = ts_mass_average(Inlet{1},'r_nondim',2);
        [~, J] = min(abs(r_nondim - plot_2D));
        
        CLT = CLT_Blade2CLTCirc(Blade{1}.x(:,J),...
            Blade{1}.rt(:,J),0);
        iLE = CLT.surfaces.iLE; iTE = CLT.surfaces.iTE;
        F{m}.chi_1_2D = CLT.Geom.inlet_angle;
        F{m}.loc_2D = plot_2D;
        F{m}.s_p = CLT.surfaces.S1.s; F{m}.s_s = CLT.surfaces.S2.s;
        
        F{m}.Cp_s = zeros(length(Inlet),length(F{m}.s_s));
        F{m}.Cp_p = zeros(length(Inlet),length(F{m}.s_p));
        F{m}.Mis_s = zeros(length(Inlet),length(F{m}.s_s));
        F{m}.Mis_p = zeros(length(Inlet),length(F{m}.s_p));
        F{m}.Inc_2D = zeros(length(Inlet),1);
        F{m}.M_2D = zeros(length(Inlet),1);
        
        % Coordinates for streamsurfaces later
        F{m}.x_2D = Inlet{1}.x(J,1);
        F{m}.r_2D = Inlet{1}.r(J,1);
        F{m}.xrrt_2D_LE = [Blade{1}.x(iLE,J) - 1e-5 Blade{1}.r(iLE,J) Blade{1}.rt(iLE,J)];
            
        for n = 1:length(Inlet)
            % Pressures for Cp or Yp
            if Inlet{n}.bv.rpm == 0
                [Po1 mass_r_in] = ts_mass_average(Inlet{n}, 'Po', 2);
                P1 = ts_mass_average(Inlet{n}, 'P', 2);
                Po1_mid = sum(Po1(J).*(mass_r_in(J))) / sum(mass_r_in(J));
                P1_mid = sum(P1(J).*(mass_r_in(J))) / sum(mass_r_in(J));
            else
                [Po1_rel mass_r_in] = ts_mass_average(Inlet{n}, 'Po_rel', 2);
                P1 = ts_mass_average(Inlet{n}, 'P', 2);
                Po1_rel_mid = sum(Po1_rel(J).*(mass_r_in(J))) /...
                    sum(mass_r_in(J));
                P1_mid = sum(P1(J).*(mass_r_in(J))) / sum(mass_r_in(J));
            end
            
            if plot_entropy_2D == 0
                % Yp parameters
                if Inlet{n}.bv.rpm == 0
                    Po2 = ts_mass_average(Outlet{n}, 'Po', 2);
                    Po2_mid = sum(Po2(J).*(mass_r_in(J))) / sum(mass_r_in(J));
                    Loss_mid(n) = (Po1_mid - Po2_mid) ./ (Po1_mid - P1_mid);
                else
                    Po2_rel = ts_mass_average(Outlet{n}, 'Po_rel', 2);
                    Po2_rel_mid = sum(Po2_rel(J).*(mass_r_in(J))) /...
                        sum(mass_r_in(J));
                    Loss_mid(n) = (Po1_rel_mid - Po2_rel_mid) ./ (Po1_rel_mid - P1_mid);
                end
            else
                % Zeta parameters
                [V1 mass_r_in] = ts_mass_average(Inlet{n}, 'V_rel', 2);
                s2 = ts_mass_average(Outlet{n}, 's', 2);
                s1 = ts_mass_average(Inlet{n}, 's', 2);
                T2 = ts_mass_average(Outlet{n}, 'T', 2);
                V1_mid = sum(V1(J).*(mass_r_in(J))) / sum(mass_r_in(J));
                s2_mid = sum(s2(J).*(mass_r_in(J))) / sum(mass_r_in(J));
                s1_mid = sum(s1(J).*(mass_r_in(J))) / sum(mass_r_in(J));
                T2_mid = sum(T2(J).*(mass_r_in(J))) / sum(mass_r_in(J));
                
                if plot_break == 1
                    filename = F{m}.name{n};
                    if exist([directory 'Break_mid_' strrep(filename,'.hdf5','.mat')], 'file') == 0
                        if isempty(G{m}{n}) == 1
                            if exist([directory strrep(filename, '.hdf5', '.mat')],'file') == 0
                                g = ts_read_hdf5([directory filename]);
                            else
                                load([directory strrep(filename, '.hdf5', '.mat')]);
                            end
                        else
                            g = G{m}{n};
                        end
                        [Loss Surface h] = ts_break_loss_streamsurface(g, Inlet{n}, Outlet{n},...
                            normal, 1, 1, 'midheight');
                        close(h.domain)
                        close(h.mass)
                        close(h.bl)
                        save([directory 'Break_mid_' strrep(filename,'.hdf5','.mat')],...
                            'Loss','Surface','-v7.3');
                    else
                        load([directory 'Break_mid_' strrep(filename,'.hdf5','.mat')])
                    end
                    
                    Loss_mid(n) = Loss.Zeta.mid * Loss.inlet.mass / Loss.mid.mass_inlet;
                elseif plot_break == 2
                    pv_directory = strrep(strrep(directory,'TURBOSTREAM','PARAVIEW'),'TS','PV');
                    pv_filename = strrep(F{m}.name{n},'.hdf5','_paraview.hdf5');
                    if exist([pv_directory ...
                            strrep(pv_filename,'.hdf5','_hub_mid_cas_loss.mat')],'file') ~= 0
                        load([pv_directory ...
                            strrep(pv_filename,'.hdf5','_hub_mid_cas_loss.mat')])
                        Loss_mid(n) = 2 * (Mid.s2 - Mid.s1) .* Mid.T2 ./ (Mid.V1.^2);
                    else
                        Loss_mid = nan;
                    end
%                     if exist([pv_directory strrep(pv_filename,'.hdf5','_2D_cut_loss.mat')],'file') ~= 0
%                         load([pv_directory strrep(pv_filename,'.hdf5','_2D_cut_loss.mat')])
%                         Loss_mid(n) = Cut.Loss;
%                     else
%                         Loss_mid = nan;
%                     end
                else
                    Loss_mid(n) = 2 * (s2_mid - s1_mid) .* T2_mid ./ (V1_mid .* V1_mid);
                end
            end
            
            % Cp Plot
            P = Blade{n}.P(:,J);
            P_mid = sum(P.*repmat((mass_r_in(J)),size(P,1),1),2) / sum(mass_r_in(J));
            x_mid = Blade{n}.x(:,J);
            if iLE > iTE
                Pp = P_mid(iLE:-1:iTE);
                Ps = P_mid([iLE:end-1 1:iTE]);
                xp = x_mid(iLE:-1:iTE);
            else
                Ps = P_mid(iLE:iTE);
                Pp = P_mid([iLE:-1:1 end-1:-1:iTE]);
                xp = x_mid([iLE:-1:1 end-1:-1:iTE]);
            end
            
            Mf = 0.001:0.001:2; P_func = (1 + (ga-1)*Mf.^2/2).^(-ga/(ga-1));
            F{m}.Mis_s(n,:) = interp1(P_func,Mf,Ps/Po1_mid);
            F{m}.Mis_p(n,:) = interp1(P_func,Mf,Pp/Po1_mid);
            
            F{m}.Cp_s(n,:) = (Ps - P1_mid) / (Po1_mid - P1_mid);
            F{m}.Cp_p(n,:) = (Pp - P1_mid) / (Po1_mid - P1_mid);
            
            % Incidences
            Alpha_1 = ts_mass_average(Inlet{n}, 'Alpha', 2);
            M_1 = ts_mass_average(Inlet{n}, 'M', 2);
            F{m}.Inc_2D(n) = Alpha_1(J) - F{m}.chi_1_2D;
            F{m}.M_2D(n) = M_1(J);
        end

        F{m}.Loss_mid = Loss_mid;
        
        % Loss Loop plot
        figure(H.mid)
        hold on
        p(m) = plot3(F{m}.Inc_2D,F{m}.Loss_mid,F{m}.number,'k-','LineWidth',2,'MarkerSize',20);
        set(p(m),'Color',F{m}.colour);
        M{m} = ['M_{in} = ' num2str(F{m}.M) ' - TS'];
        
        % Cp plot
        figure(H.Cp_mid)
        hold on
        plot(F{m}.s_s,-F{m}.Cp_s,'Color',F{m}.colour);
        plot(F{m}.s_p,-F{m}.Cp_p,'Color',F{m}.colour);
        
        % Read in MISES data to plot on the same graph
        [MPolarx , ~, ~,]=ReadPolarx(ext,[path 'M' num2str(num2str(F{m}.M)) '/']);
        if isstruct(MPolarx) == 1 
            MI_P1_Po1 = MPolarx.p1pt; MI_P2_Po1 = MPolarx.p2pt;
            MI_M1 = MPolarx.minl; MI_M2 = MPolarx.mout;
    %         MI_Re = MIses.reyn; 
            F{m}.Alpha_mis = MPolarx.binl;
    %         MI_chi1 = atand(MBlade.sinl);
            MI_Yp = MPolarx.omega;


            % Variables to calculate entropy generation coefficient
            MI_P2_P1 = MI_P2_Po1 ./ MI_P1_Po1;
            MI_T2_To2 = 1./(1 + (ga-1)/2 * MI_M2.^2);
            MI_T1_To1 = 1./(1 + (ga-1)/2 * MI_M1.^2);
            MI_T2_T1 = MI_T2_To2 ./ MI_T1_To1;

            MI_Vfun = (ga-1)^0.5 * MI_M1 .* (1 + (ga-1)/2 * MI_M1.^2).^(-0.5);
            MI_dS = cp * log(MI_T2_T1) - R * log(MI_P2_P1);
            MI_Zeta = 2 * MI_dS .* MI_T2_To2 .* (MI_Vfun).^(-2) / cp;
            F{m}.Inc_mis = F{m}.Alpha_mis - F{m}.chi_1;

            if plot_entropy_2D == 1
                F{m}.Loss_mis = MI_Zeta;
            else
                F{m}.Loss_mis = MI_Yp;
            end

            % Plot mises loss loops
            figure(H.mid)
            p(m+length(F)) = plot(F{m}.Inc_mis,F{m}.Loss_mis,'--','LineWidth',2);
            set(p(m+length(F)),'Color',F{m}.colour)   
            M{m+length(F)} = ['M_{in} = ' num2str(F{m}.M) ' - MISES'];

            % Plot MISES Cp distributions
            figure(H.Cp_mid)
            F{m}.Cp_s_mis = zeros(length(F{m}.Inc_mis),length(MPolarx.ileb(1):MPolarx.iteb(1)));
            F{m}.Cp_p_mis = zeros(length(F{m}.Inc_mis),length(MPolarx.ileb(2):MPolarx.iteb(2)));
            Cp_s_mis = zeros(length(Inlet),length(MPolarx.ileb(1):MPolarx.iteb(1)));
            Cp_p_mis = zeros(length(Inlet),length(MPolarx.ileb(2):MPolarx.iteb(2)));
            th_p_mis = zeros(length(Inlet),length(MPolarx.ileb(2):MPolarx.iteb(2)));
            u_p_mis = zeros(length(Inlet),length(MPolarx.ileb(2):MPolarx.iteb(2)));
            F{m}.Inc_mis_close = zeros(length(Inlet),1);

            F{m}.s_s_mis = MPolarx.s{1}(MPolarx.ileb(1):MPolarx.iteb(1),1).';
            F{m}.s_p_mis = MPolarx.s{1}(MPolarx.ileb(2):MPolarx.iteb(2),2).';

            % Plot closest MISES Cp plot to 3D loop
            for n = 1:length(F{m}.Inc)
                [~, n_m] = min(abs(F{m}.Inc(n) - F{m}.Inc_mis));
                Cp_s_mis(n,:) = MPolarx.cp{n_m}(MPolarx.ileb(1):MPolarx.iteb(1),1).';
                Cp_p_mis(n,:) = MPolarx.cp{n_m}(MPolarx.ileb(2):MPolarx.iteb(2),2).';
                th_p_mis(n,:) = MPolarx.th{n_m}(MPolarx.ileb(2):MPolarx.iteb(2),2).';
                u_p_mis(n,:) = MPolarx.uedg{n_m}(MPolarx.ileb(2):MPolarx.iteb(2),2).';
                F{m}.Inc_mis_close(n) = F{m}.Inc_mis(n_m);
            end

            th_p_te_mis = th_p_mis(:,end);
            u_p_te_mis = u_p_mis(:,end);

            plot(F{m}.s_s_mis,-Cp_s_mis,'--','Color',F{m}.colour);
            plot(F{m}.s_p_mis,-Cp_p_mis,'--','Color',F{m}.colour);

            % Record all MISES Cp plots
            for n = 1:length(F{m}.Inc_mis)
                F{m}.Cp_s_mis(n,:) = MPolarx.cp{n}(MPolarx.ileb(1):MPolarx.iteb(1),1).';
                F{m}.Cp_p_mis(n,:) = MPolarx.cp{n}(MPolarx.ileb(2):MPolarx.iteb(2),2).';
            end

            % Boundary layer parameters for turbostream
            s_te_mis = MPolarx.s{1}(MPolarx.iteb(2),2);
            x_te_ts = interp1(CLT.surfaces.S1.s,xp,s_te_mis);
            th_p_te_ts = zeros(size(F{m}.number));

    %         for n = 1:length(F{m}.number)
    %             load([directory strrep(F{m}.name{n},'.hdf5','.mat')])
    %             cut = ts_secondary(ts_structunstruct_cut(g,[x_te_ts 0 ; x_te_ts 1],[0:6],2,1));
    %             V_mid = cut.V(J,round(size(cut.x,2)/2):end); 
    %             ro_mid = cut.ro(J,round(size(cut.x,2)/2):end); 
    %             rt_mid = cut.rt(J,round(size(cut.x,2)/2):end);
    %             V_mid(isnan(V_mid)) = interp1(find(isnan(V_mid) == 0),V_mid(isnan(V_mid) == 0),...
    %                 find(isnan(V_mid) == 1),'linear','extrap');
    %             ro_mid(isnan(ro_mid)) = interp1(find(isnan(ro_mid) == 0),ro_mid(isnan(ro_mid) == 0),...
    %                 find(isnan(ro_mid) == 1),'linear','extrap');
    %             rt_mid(isnan(rt_mid)) = interp1(find(isnan(rt_mid) == 0),rt_mid(isnan(rt_mid) == 0),...
    %                 find(isnan(rt_mid) == 1),'linear','extrap');
    % 
    %             [rt_mid i] = unique(double(rt_mid)); 
    %             ro_mid = double(ro_mid(i)); 
    %             V_mid = double(V_mid(i));
    % 
    %             rt_midf = interp1(linspace(1,length(rt_mid)*5,length(rt_mid)),...
    %                 rt_mid,1:length(rt_mid)*5,'pchip');
    %             ro_mid = interp1(rt_mid,ro_mid,rt_midf);
    %             V_mid = interp1(rt_mid,V_mid,rt_midf);
    %             rt_mid = rt_midf;
    %             [V_e k] = findpeaks(V_mid);
    %             ro_e = ro_mid(k);
    %             rt_e = rt_mid(k);
    %             
    % %             th_p_te_ts(n) = sum( diff(rt_mid(k-1:end)) .* (1 - V_mid(k:end) / V_e) .*...
    % %                 ro_mid(k:end) .* V_mid(k:end) / (ro_e * V_e) );
    %             
    %         end
        else
            M{m+length(F)} = ['M_{in} = ' num2str(F{m}.M) ' - MISES'];
            
        end
            
    end

    figure(H.Cp_mid)
    title([filestring(1:2) ' Mid Height Cp Distributions'])
    xlabel('s / c','FontSize',16);
    ylabel('C_p','FontSize',16);
    
    figure(H.mid)
    legend(p,M);
    title([filestring(1:2) ' Mid Height Loss Loops'])
    xlabel('\alpha_{1,mid} - \chi_{1,mid}','FontSize',16);
    if plot_entropy_2D == 1
        ylabel('\zeta_{s,mid}','FontSize',16);
    else
        ylabel('Y_{p,mid}','FontSize',16);
    end
    
    axis auto
    [~, n_min] = min(abs(F{m_des}.Alpha_mid - des.Alpha));
    Loss_min = F{m_des}.Loss_mid(n_min);
    V = axis; V(3) = 0; V(4) = 5*min(Loss_min); axis(V);
    
    % Interpolate mises failure incidences
    [~, n_min] = min(abs(F{m_des}.Alpha_mis - des.Alpha));
    Loss_min_mis = F{m_des}.Loss_mis(n_min);
    
    for m = 1:length(F)
    % Determine failure incidences
%     [Loss_min n_min] = min(F{m}.Loss);
    if isfield(F{m},'Loss_mis') == 1
        [~, n_min] = min(F{m}.Loss_mis);

        Loss_1 = F{m}.Loss_mis/Loss_min_mis; Loss_2 = F{m}.Loss_mis/Loss_min_mis;
        Loss_1(n_min:end) = 0; Loss_2(1:n_min) = 0;
        [Fail_Loss_1 n] = unique(Loss_1); Fail_Inc_1 = F{m}.Inc_mis(n);
        [Fail_Loss_2 n] = unique(Loss_2); Fail_Inc_2 = F{m}.Inc_mis(n);

        [Fail_Inc_1 n] = sort(Fail_Inc_1); Fail_Loss_1 = Fail_Loss_1(n);
        [Fail_Inc_2 n] = sort(Fail_Inc_2); Fail_Loss_2 = Fail_Loss_2(n);

        if length(Fail_Loss_1) > 2
            F{m}.Inc_mis_fail(1) = interp1(Fail_Loss_1,Fail_Inc_1,1.5,'pchip','extrap');
        else
            F{m}.Inc_mis_fail(1) = min(F{m}.Alpha_mid);
        end
        if length(Fail_Loss_2) > 2
            F{m}.Inc_mis_fail(2) = interp1(Fail_Loss_2,Fail_Inc_2,1.5,'pchip','extrap');
        else
            F{m}.Inc_mis_fail(2) = max(F{m}.Alpha_mid);
        end
        
        % Plot on failure incidences
        plot(F{m}.Inc_mis_fail(1), min(Loss_min_mis)*1.5,'o','Color',F{m}.colour)
        plot(F{m}.Inc_mis_fail(2), min(Loss_min_mis)*1.5,'o','Color',F{m}.colour)
    end
    end

else
    if ishandle(H.mid)
        close(H.mid)    
    end
    if ishandle(H.Cp_mid)
        close(H.Cp_mid)
    end
end

%% Plot incidence mach number diagram
if length(F) > 1

    % Get failure boundaries
    Inc_1 = zeros(size(F)); Inc_2 = zeros(size(F)); M_mid = zeros(size(F));
    for n = 1:length(F)
        Inc_1(n) = F{n}.Inc_fail(1); Inc_2(n) = F{n}.Inc_fail(2);
        M_mid(n) = F{n}.M;
    end
    
    % Sub NaN's for out of range results
    figure(H.loop)
    V = axis;
    
    Inc_1(Inc_1 < V(1) | Inc_1 > V(2)) = nan; Inc_2(Inc_2 < V(1) | Inc_2 > V(2)) = nan;
    
    % Plot on inc - mach plot
    figure(H.IM)
    hold on
    plot(M_mid,Inc_1,'k.-')
    plot(M_mid,Inc_2,'k.-')
    
    % Plot mises failure boundaries if available
    if plot_2D ~= -1
        Inc_1 = zeros(size(F)); Inc_2 = zeros(size(F)); M_mid = zeros(size(F));
        for n = 1:length(F)
            if isfield(F{n},'Inc_mis_fail') == 1
            Inc_1(n) = F{n}.Inc_mis_fail(1); Inc_2(n) = F{n}.Inc_mis_fail(2);
            M_mid(n) = F{n}.M;
            end
        end
        
        % Sub NaN's for out of range results
        figure(H.loop)
        V = axis;

        Inc_1(Inc_1 < V(1) | Inc_1 > V(2)) = nan; Inc_2(Inc_2 < V(1) | Inc_2 > V(2)) = nan;

        % Plot on inc - mach plot
        figure(H.IM)
        hold on
        plot(M_mid,Inc_1,'k.--')
        plot(M_mid,Inc_2,'k.--')
    end
    
    title([filestring(1:2) ' Incidence - Mach Number Envelope'])
    xlabel('M_{mid}','FontSize',16);
    ylabel('\alpha_{1,mid} - \chi_{1,mid}','FontSize',16);

else
    if ishandle(H.IM)
        close(H.IM)
    end
end

%% Calculate height of hub and casing boundary layers and theta coordinates 
%% of stagnation streamlines
% vorticity_limit = 700;
% frac_r = 0.05;
if plot_break == 2
for m = 1:length(F)
    if exist([directory strrep(filestring, '*', [num2str(F{m}.M*100) '00_bls.dat'])],'file') == 0
    r_hub = zeros(size(F{m}.number)); r_cas = zeros(size(F{m}.number)); 
    x_hub = zeros(size(F{m}.number)); x_cas = zeros(size(F{m}.number)); 
    rt_stag_hub = zeros(size(F{m}.number)); rt_stag_cas = zeros(size(F{m}.number)); 
    for n = 1:length(F{m}.number)
        Inlet = ts_secondary(F{m}.Inlet{n}); Blade = ts_secondary(F{m}.Blade{n});
        Vx1 = ts_mass_average(Inlet, 'Vx', 2);
        rovx1 = ts_mass_average(Inlet, 'rovx', 2);
        r1 = ts_mass_average(Inlet, 'r', 2);
        r1_nondim = (r1 - min(r1)) / (max(r1) - min(r1));
        J_mid_30_in = find(r1_nondim>(0.35) & r1_nondim<(0.65));
        
        % Vorticity based boundary layer search
%         dr1 = diff(r1);
%         r1_av = conv2(r1,[0.5 0.5]); r1_av = r1_av(2:end-1);
%         dVx1_dr1 = diff(Vx1) ./ dr1;
%         
%         Vort_1 = (dVx1_dr1(2:end) - vorticity_limit) .* (dVx1_dr1(1:end-1) - vorticity_limit);
%         J_Vort_1 = find(Vort_1 < 0); J_Vort_1(J_Vort_1>min(J_mid_30_in)) = [];
%         J_hub = max(J_Vort_1);
% 
%         r_hub_in = interp1(dVx1_dr1(J_hub-3:J_hub+3),r1_av(J_hub-3:J_hub+3),vorticity_limit);
%         r_hub(n) = r_hub_in + frac_r*(max(r1)-min(r1));
%         
%         Vort_2 = (dVx1_dr1(2:end) + vorticity_limit) .* (dVx1_dr1(1:end-1) + vorticity_limit);
%         J_Vort_2 = find(Vort_2 < 0); J_Vort_2(J_Vort_2 < max(J_mid_30_in)) = [];
%         J_cas = min(J_Vort_2);
% 
%         r_cas_in = interp1(dVx1_dr1(J_cas-3:J_cas+3),r1_av(J_cas-3:J_cas+3),-vorticity_limit);
%         r_cas(n) = r_cas_in - frac_r*(max(r1)-min(r1));
%         
%         [~, J_hub] = min(abs(r1 - r_hub(n)));
%         [~, J_cas] = min(abs(r1 - r_cas(n)));
        
        % Mass based boundary search
        rovx_av = 0.5*(rovx1(1:end-1) + rovx1(2:end));
        dr = diff(r1);
        mass = [0 cumsum(rovx_av .* dr)]; mass = mass / max(mass);
%         plot(mass,r1,'-','Color',[(n-10)/30 0 0])
        [~,J_hub] = min(abs(mass - 0.2)); [~,J_cas] = min(abs(mass - 0.8)); 
        r_hub(n) = interp1(mass,r1,0.2,'pchip');
        r_cas(n) = interp1(mass,r1,0.8,'pchip');
        x_hub(n) = Inlet.x(J_hub,1); x_cas(n) = Inlet.x(J_cas,1);
        
        % Angle based, stagnation streamline search
%         mean_x_blade = (max(Blade.x(:,J_hub)) - min(Blade.x(:,J_hub))) / 3 + min(Blade.x(:,J_hub));
%         V_blade = Blade.V(Blade.x(:,J_hub) < mean_x_blade,J_hub);
%         x_blade = Blade.x(Blade.x(:,J_hub) < mean_x_blade,J_hub);
%         rt_blade = Blade.rt(Blade.x(:,J_hub) < mean_x_blade,J_hub);
%         [~, i] = min(abs(V_blade));
%         x = x_blade(i); rt = rt_blade(i);
%         rt_stag_hub(n) = rt - (Inlet.Vt_rel(J_hub,1)/Inlet.Vx(J_hub,1)) * (x - Inlet.x(J_hub,1));
%         
%         mean_x_blade = (max(Blade.x(:,J_cas)) - min(Blade.x(:,J_cas))) / 3 + min(Blade.x(:,J_cas));
%         V_blade = Blade.V(Blade.x(:,J_cas) < mean_x_blade,J_cas);
%         x_blade = Blade.x(Blade.x(:,J_cas) < mean_x_blade,J_cas);
%         rt_blade = Blade.rt(Blade.x(:,J_cas) < mean_x_blade,J_cas);
%         [~, i] = min(abs(V_blade));
%         x = x_blade(i); rt = rt_blade(i);
%         rt_stag_cas(n) = rt - (Inlet.Vt_rel(J_cas,1)/Inlet.Vx(J_cas,1)) * (x - Inlet.x(J_cas,1));
        
        % Use paraview to track back the stagnation line to find rt stags
        % else 
        % at inlet
        CLT = CLT_Blade2CLTCirc(Blade.x(:,J_hub),Blade.rt(:,J_hub),0);
        iLE = CLT.surfaces.iLE; 
        xrrt = [Blade.x(iLE,J_hub) - 1e-5 Blade.r(iLE,J_hub) Blade.rt(iLE,J_hub)];
        [status, result] = system(['pvpython' ...
            ' /home/jvt24/Documents/PARAVIEW/ts_single_line.py "' pv_directory '" "'...
            pv_filename '" ' num2str(xrrt) ' ' num2str(x_hub(n))]);
        C = read_paraview([pv_directory 'x_intersect.csv'],1,0);
        delete([pv_directory 'x_intersect.*'])
        rt_stag_hub(n) = C.rt;
            
%         mean_x_blade = (max(Blade.x(:,J_hub)) - min(Blade.x(:,J_hub))) / 3 + min(Blade.x(:,J_hub));
%         V_blade = Blade.V(Blade.x(:,J_hub) < mean_x_blade,J_hub);
%         x_blade = Blade.x(Blade.x(:,J_hub) < mean_x_blade,J_hub);
%         r_blade = Blade.r(Blade.x(:,J_hub) < mean_x_blade,J_hub);
%         rt_blade = Blade.rt(Blade.x(:,J_hub) < mean_x_blade,J_hub);
%         [~, i] = min(abs(V_blade));
%         x = x_blade(i); r = r_blade(i); rt = rt_blade(i);
%         rt_stag_hub(n) = rt - (Inlet.Vt_rel(J_hub,1)/Inlet.Vx(J_hub,1)) * (x - Inlet.x(J_hub,1));
        
        CLT = CLT_Blade2CLTCirc(Blade.x(:,J_cas),Blade.rt(:,J_cas),0);
        iLE = CLT.surfaces.iLE; 
        xrrt = [Blade.x(iLE,J_cas) - 1e-5 Blade.r(iLE,J_cas) Blade.rt(iLE,J_cas)];
        [status, result] = system(['pvpython' ...
            ' /home/jvt24/Documents/PARAVIEW/ts_single_line.py "' pv_directory '" "'...
            pv_filename '" ' num2str(xrrt) ' ' num2str(x_cas(n))]);
        C = read_paraview([pv_directory 'x_intersect.csv'],1,0);
        delete([pv_directory 'x_intersect.*'])
        rt_stag_cas(n) = C.rt;        
%         mean_x_blade = (max(Blade.x(:,J_cas)) - min(Blade.x(:,J_cas))) / 3 + min(Blade.x(:,J_cas));
%         V_blade = Blade.V(Blade.x(:,J_cas) < mean_x_blade,J_cas);
%         x_blade = Blade.x(Blade.x(:,J_cas) < mean_x_blade,J_cas);
%         r_blade = Blade.r(Blade.x(:,J_cas) < mean_x_blade,J_cas);
%         rt_blade = Blade.rt(Blade.x(:,J_cas) < mean_x_blade,J_cas);
%         [~, i] = min(abs(V_blade));
%         x = x_blade(i); r = r_blade(i); rt = rt_blade(i);
%         rt_stag_cas(n) = rt - (Inlet.Vt_rel(J_cas,1)/Inlet.Vx(J_cas,1)) * (x - Inlet.x(J_cas,1));
        
    end
    
    % Robust regression to fit a quadratic
%     p = fit(F{m}.Inc, r_hub, 'poly2', 'Normalize', 'on', 'Robust', 'on');
%     r_hub = feval(p,F{m}.Inc);
%     p = fit(F{m}.Inc, r_cas, 'poly2', 'Normalize', 'on', 'Robust', 'on');
%     r_cas = feval(p,F{m}.Inc);
    
    F{m}.x_hub_bl = x_hub;
    F{m}.x_cas_bl = x_cas;
    F{m}.r_hub_bl = r_hub;
    F{m}.r_cas_bl = r_cas;
    
%     p = fit(F{m}.Inc, rt_stag_hub, 'poly2', 'Normalize', 'on', 'Robust', 'on');
%     rt_stag_hub = feval(p,F{m}.Inc);
%     p = fit(F{m}.Inc, rt_stag_cas, 'poly2', 'Normalize', 'on', 'Robust', 'on');
%     rt_stag_cas = feval(p,F{m}.Inc);
    F{m}.rt_stag_hub = rt_stag_hub;
    F{m}.rt_stag_cas = rt_stag_cas;

    % Write out coordinates
    dlmwrite([directory strrep(filestring, '*', [num2str(F{m}.M*100) '00_bls.dat'])],...
        [F{m}.number x_hub r_hub rt_stag_hub x_cas r_cas rt_stag_cas],'precision',10)
    else
        BLs = dlmread([directory strrep(filestring, '*', [num2str(F{m}.M*100) '00_bls.dat'])]);
        F{m}.x_hub_bl = BLs(:,2);
        F{m}.r_hub_bl = BLs(:,3);
        F{m}.rt_stag_hub = BLs(:,4);
        F{m}.x_cas_bl = BLs(:,5);
        F{m}.r_cas_bl = BLs(:,6);
        F{m}.rt_stag_cas = BLs(:,7);
    end
end

for m = 1:length(F)
    f = F{m};
    save([directory strrep(filestring, '*', [num2str(F{m}.M*100) '00_cuts.mat'])],'f','-v7.3')
end
end

%% Plot loss loops with streamline breakdown
if plot_break == 1
    % Plot full 3D loss break down loops using streamline tracking
    for m = 1:length(F)
        Inlet = ts_secondary(F{m}.Inlet);
        Outlet = ts_secondary(F{m}.Outlet);
        Losses = cell(size(Inlet)); Surfaces = cell(size(Inlet)); d = ceil(length(Inlet)/5);
        if sum(F{m}.M == M_plot) > 0
        F{m}.M
        if exist([directory 'Breaks_' filestring(1:end-1) num2str(F{m}.M*100) '00.mat']...
                ,'file') == 0
            for n = 1:length(Inlet)
                filename = F{m}.name{n};
                if exist([directory 'Break_single_' strrep(filename,'.hdf5','.mat')], 'file') == 0
                    if isempty(G{m}{n}) == 1
                        if exist([directory strrep(filename, '.hdf5', '.mat')],'file') == 0
                            g = ts_read_hdf5([directory filename]);
                        else
                            load([directory strrep(filename, '.hdf5', '.mat')]);
                        end
                    else
                        g = G{m}{n};
                    end

                    [Loss Surface h] = ts_break_loss_streamsurface(...
                        ts_secondary(g), Inlet{n}, Outlet{n}, normal, plot_entropy, 1);
                    
%                     save([directory 'Break_single_' strrep(filename,'.hdf5','.mat')],...
%                         'Loss','Surface','-v7.3');
%                     saveas(h.wake, [directory 'streamsurface_' strrep(filename,'.hdf5','.fig')], 'fig');
%                     saveas(h.domain, [directory 'streamlines_' strrep(filename,'.hdf5','.fig')], 'fig');
%                     close(h.domain)
%                     close(h.mass)
%                     close(h.bl)
                    close all
                else
                    load([directory 'Break_single_' strrep(filename,'.hdf5','.mat')])
%                     open([directory 'streamsurface_' strrep(filename,'.hdf5','.fig')])
                end
%                 Losses{n} = Loss; Surfaces{n} = Surface;
            end
            save([directory 'Breaks_' filestring(1:end-1) num2str(F{m}.M*100) '00.mat'],...
                'Losses','Surfaces','-v7.3');
        else
            load([directory 'Breaks_' filestring(1:end-1) num2str(F{m}.M*100) '00.mat']);
        end
        
        % Plot wakes with loss breakdown
        V = [0 0];
        for n = 1:length(Inlet)
            plot_break_wakes(Inlet{n},Outlet{n},Losses{n},Surfaces{n},plot_entropy,n,d,H)
            v = caxis;
            V(1) = min(V(1),v(1)); V(2) = max(V(2),v(2));
        end
        for n = 1:length(Inlet)
            subplot(d,5,n)
            caxis(V);
        end
        colorbar('EastOutside')
        
        if plot_entropy == 1
            F{m}.Zeta_break_total = zeros(size(Losses));
            F{m}.Zeta_hub = zeros(size(Losses));
            F{m}.Zeta_mid = zeros(size(Losses));
            F{m}.Zeta_cas = zeros(size(Losses));

            F{m}.Zeta_hub_specific = zeros(size(Losses));
            F{m}.Zeta_mid_specific = zeros(size(Losses));
            F{m}.Zeta_cas_specific = zeros(size(Losses));        

            for n = 1:size(Losses,1)
                F{m}.Zeta_break_total(n) = Losses{n}.Zeta.total;
                F{m}.Zeta_hub(n) = Losses{n}.Zeta.hub;
                F{m}.Zeta_mid(n) = Losses{n}.Zeta.mid;
                F{m}.Zeta_cas(n) = Losses{n}.Zeta.cas;

                F{m}.Zeta_hub_specific(n) = Losses{n}.Zeta.hub * ...
                    Losses{n}.inlet.mass / Losses{n}.hub.mass_inlet;
                F{m}.Zeta_mid_specific(n) = Losses{n}.Zeta.mid * ...
                    Losses{n}.inlet.mass / Losses{n}.mid.mass_inlet;
                F{m}.Zeta_cas_specific(n) = Losses{n}.Zeta.cas * ...
                    Losses{n}.inlet.mass / Losses{n}.cas.mass_inlet;            
            end

            figure(H.break)
            hold on

            plot3(F{m}.Inc, F{m}.Zeta_hub, F{m}.number, '--',...
                'LineWidth',2,'Color',F{m}.colour);
            plot3(F{m}.Inc, F{m}.Zeta_hub+F{m}.Zeta_mid, F{m}.number, '--',...
                'LineWidth',2,'Color',F{m}.colour);
            plot3(F{m}.Inc, F{m}.Zeta_hub+F{m}.Zeta_mid+F{m}.Zeta_cas, F{m}.number, '-',...
                'LineWidth',2,'Color',F{m}.colour);

            V = axis; V(3) = 0; V(4) = 3*min(F{m}.Zeta_break_total);
            axis(V);  
            ylabel('\zeta_{s}','FontSize',16);
            xlabel('\alpha_{1,mid} - \chi_{1,mid}','FontSize',16);

            V4 = max([3*min(F{m}.Zeta_hub) 3*min(F{m}.Zeta_mid) 3*min(F{m}.Zeta_cas)]);

            figure(H.break_2)
            subplot(1,3,1)
            hold on
            plot3(F{m}.Inc, F{m}.Zeta_hub, F{m}.number, '-',...
                'LineWidth',2,'Color',F{m}.colour);
            V = axis; V(3) = 0; V(4) = V4;
            axis(V);  
            ylabel('\zeta_{s}','FontSize',16);
            xlabel('\alpha_{1,mid} - \chi_{1,mid}','FontSize',16);
            title('Hub Loss - 3D','FontSize',16)

            subplot(1,3,2)
            hold on
            plot3(F{m}.Inc, F{m}.Zeta_mid, F{m}.number, '-',...
                'LineWidth',2,'Color',F{m}.colour);
            V = axis; V(3) = 0; V(4) = V4;
            axis(V);  
            ylabel('\zeta_{s}','FontSize',16);
            xlabel('\alpha_{1,mid} - \chi_{1,mid}','FontSize',16);
            title('Mid Loss - 2D','FontSize',16)

            subplot(1,3,3)
            hold on
            plot3(F{m}.Inc, F{m}.Zeta_cas, F{m}.number, '-',...
                'LineWidth',2,'Color',F{m}.colour);
            V = axis; V(3) = 0; V(4) = V4;
            axis(V);  
            ylabel('\zeta_{s}','FontSize',16);
            xlabel('\alpha_{1,mid} - \chi_{1,mid}','FontSize',16);
            title('Casing Loss - 3D','FontSize',16)
        else
            F{m}.Yp_break_total = zeros(size(Losses));
            F{m}.Yp_hub = zeros(size(Losses));
            F{m}.Yp_mid = zeros(size(Losses));
            F{m}.Yp_cas = zeros(size(Losses));

            F{m}.Yp_hub_specific = zeros(size(Losses));
            F{m}.Yp_mid_specific = zeros(size(Losses));
            F{m}.Yp_cas_specific = zeros(size(Losses));        

            for n = 1:size(Losses,1)
                F{m}.Yp_break_total(n) = Losses{n}.Yp.total;
                F{m}.Yp_hub(n) = Losses{n}.Yp.hub;
                F{m}.Yp_mid(n) = Losses{n}.Yp.mid;
                F{m}.Yp_cas(n) = Losses{n}.Yp.cas;

                F{m}.Yp_hub_specific(n) = Losses{n}.Yp.hub * ...
                    Losses{n}.inlet.mass / Losses{n}.hub.mass_inlet;
                F{m}.Yp_mid_specific(n) = Losses{n}.Yp.mid * ...
                    Losses{n}.inlet.mass / Losses{n}.mid.mass_inlet;
                F{m}.Yp_cas_specific(n) = Losses{n}.Yp.cas * ...
                    Losses{n}.inlet.mass / Losses{n}.cas.mass_inlet;            
            end

            figure(H.break)
            hold on

            plot3(F{m}.Inc, F{m}.Yp_hub, F{m}.number, '--',...
                'LineWidth',2,'Color',F{m}.colour);
            plot3(F{m}.Inc, F{m}.Yp_hub+F{m}.Yp_mid, F{m}.number, '--',...
                'LineWidth',2,'Color',F{m}.colour);
            plot3(F{m}.Inc, F{m}.Yp_hub+F{m}.Yp_mid+F{m}.Yp_cas, F{m}.number, '-',...
                'LineWidth',2,'Color',F{m}.colour);

            V = axis; V(3) = 0; V(4) = 3*min(F{m}.Yp_break_total);
            axis(V);  
            ylabel('Y_p','FontSize',16);
            xlabel('\alpha_{1,mid} - \chi_{1,mid}','FontSize',16);

            V4 = max([3*min(F{m}.Yp_hub) 3*min(F{m}.Yp_mid) 3*min(F{m}.Yp_cas)]);

            figure(H.break_2)
            subplot(1,3,1)
            hold on
            plot3(F{m}.Inc, F{m}.Yp_hub, F{m}.number, '-',...
                'LineWidth',2,'Color',F{m}.colour);
            V = axis; V(3) = 0; V(4) = V4;
            axis(V);  
            ylabel('Y_p','FontSize',16);
            xlabel('\alpha_{1,mid} - \chi_{1,mid}','FontSize',16);
            title('Hub Loss - 3D','FontSize',16)

            subplot(1,3,2)
            hold on
            plot3(F{m}.Inc, F{m}.Yp_mid, F{m}.number, '-',...
                'LineWidth',2,'Color',F{m}.colour);
            V = axis; V(3) = 0; V(4) = V4;
            axis(V);  
            ylabel('Y_p','FontSize',16);
            xlabel('\alpha_{1,mid} - \chi_{1,mid}','FontSize',16);
            title('Mid Loss - 2D','FontSize',16)

            subplot(1,3,3)
            hold on
            plot3(F{m}.Inc, F{m}.Yp_cas, F{m}.number, '-',...
                'LineWidth',2,'Color',F{m}.colour);
            V = axis; V(3) = 0; V(4) = V4;
            axis(V);  
            ylabel('Y_p','FontSize',16);
            xlabel('\alpha_{1,mid} - \chi_{1,mid}','FontSize',16);
            title('Casing Loss - 3D','FontSize',16)            
        end
        end
    end

elseif plot_break == 2
    % Run streamlines through paraview to pull out hub and casing 3D
    % losses, 2D losses in 1% of span and critical surface streamlines
    load([directory strrep(F{1}.name{1},'.hdf5','.mat')])
    pv_directory = strrep(strrep(directory,'TURBOSTREAM','PARAVIEW'),'TS','PV');
    H.dCp_limit = figure(); hold on;
    for m = 1:length(F)
    if sum(F{m}.M == M_plot) > 0
        F{m}.M
        delete([pv_directory 'x_intersect.*'])
        for n = 1:length(F{m}.number)
            pv_filename = strrep(F{m}.name{n},'.hdf5','_paraview.hdf5');
            
            % Calculate surface streamlines if necessary
            if exist([pv_directory strrep(pv_filename,'.hdf5','_hub_surflines0.csv')],'file') == 0 &&...
                exist([pv_directory strrep(pv_filename,'.hdf5','_hub_surflines0.mat')],'file') == 0
            [status, result] = system(['pvpython /home/jvt24/Documents/PARAVIEW/ts_surflines.py ' ...
                '"' pv_directory '" "' pv_filename '"']);
            end
            % Calculate hub and casing dividing streamlines
            if exist([pv_directory strrep(pv_filename,'.hdf5','_cas_bl_fine.csv')],'file') == 0 &&...
                exist([pv_directory strrep(pv_filename,'.hdf5','_cas_bl_fine.mat')],'file') == 0
            [status, result] = system(['pvpython /home/jvt24/Documents/PARAVIEW/ts_hub_mid_cas.py ' ...
                '"' pv_directory '" "' pv_filename '" ' num2str(F{m}.x_hub_bl(n)) ' ' ...
                num2str(F{m}.r_hub_bl(n)) ' ' num2str(F{m}.rt_stag_hub(n)) ' '...
                num2str(F{m}.x_cas_bl(n)) ' ' num2str(F{m}.r_cas_bl(n)) ' '...
                num2str(F{m}.rt_stag_cas(n))]);
            end
            % Calculate streamlines bounding a 2D section
            if plot_2D ~= -1 && exist([pv_directory strrep(pv_filename,'.hdf5','_2D_high_fine.csv')],'file') == 0 &&...
                exist([pv_directory strrep(pv_filename,'.hdf5','_2D_high_fine.mat')],'file') == 0
                [status, result] = system(['pvpython' ...
                    ' /home/jvt24/Documents/PARAVIEW/ts_single_line.py "' pv_directory '" "'...
                    pv_filename '" ' num2str(F{m}.xrrt_2D_LE) ' ' num2str(F{m}.x_2D)]);
                C = read_paraview([pv_directory 'x_intersect.csv'],1,0);
                delete([pv_directory 'x_intersect.*'])
                F{m}.rt_2D(n) = C.rt;
                [status, result] = system(['pvpython /home/jvt24/Documents/PARAVIEW/ts_2D_sec.py ' ...
                    '"' pv_directory '" "' pv_filename '" ' num2str(F{m}.x_2D) ' ' ...
                    num2str(F{m}.r_2D) ' ' num2str(F{m}.rt_2D(n))]); 
            end
            
            if exist([directory 'Lines_' strrep(F{m}.name{n},'.hdf5','.fig')],'file') == 0   
            hub_coarse = read_paraview([pv_directory ...
                strrep(pv_filename,'.hdf5','_hub_bl_coarse.csv')],[],1,[]);
            hub_fine = read_paraview([pv_directory ...
                strrep(pv_filename,'.hdf5','_hub_bl_fine.csv')],[],1,[]);
            cas_coarse = read_paraview([pv_directory ...
                strrep(pv_filename,'.hdf5','_cas_bl_coarse.csv')],[],1,[]);
            cas_fine = read_paraview([pv_directory ...
                strrep(pv_filename,'.hdf5','_cas_bl_fine.csv')],[],1,[]);
            
            H_lines = figure();
            ts_plot_surface(g,[0 0 0],H_lines)
            for i = 1:size(hub_coarse,3)
                plot3(hub_coarse(:,1,i),hub_coarse(:,3,i),hub_coarse(:,2,i),'r-')
            end
            for i = 1:size(hub_fine,3)
                plot3(hub_fine(:,1,i),hub_fine(:,3,i),hub_fine(:,2,i),'b-')
            end
            for i = 1:size(cas_coarse,3)
                plot3(cas_coarse(:,1,i),cas_coarse(:,3,i),cas_coarse(:,2,i),'r-')
            end
            for i = 1:size(cas_fine,3)
                plot3(cas_fine(:,1,i),cas_fine(:,3,i),cas_fine(:,2,i),'b-')
            end
            
            if plot_2D ~= -1
                low_coarse = read_paraview([pv_directory ...
                    strrep(pv_filename,'.hdf5','_2D_low_coarse.csv')],[],1,[]);
                high_coarse = read_paraview([pv_directory ...
                    strrep(pv_filename,'.hdf5','_2D_high_coarse.csv')],[],1,[]);
                low_fine = read_paraview([pv_directory ...
                    strrep(pv_filename,'.hdf5','_2D_low_fine.csv')],[],1,[]);
                high_fine = read_paraview([pv_directory ...
                    strrep(pv_filename,'.hdf5','_2D_high_fine.csv')],[],1,[]);
                
                for i = 1:size(low_coarse,3)
                    plot3(low_coarse(:,1,i),low_coarse(:,3,i),low_coarse(:,2,i),'g-')
                end
                for i = 1:size(high_coarse,3)
                    plot3(high_coarse(:,1,i),high_coarse(:,3,i),high_coarse(:,2,i),'g-')
                end
                for i = 1:size(low_fine,3)
                    plot3(low_fine(:,1,i),low_fine(:,3,i),low_fine(:,2,i),'g-')
                end
                for i = 1:size(high_fine,3)
                    plot3(high_fine(:,1,i),high_fine(:,3,i),high_fine(:,2,i),'g-')
                end
            end

            saveas(H_lines, [directory 'Lines_' strrep(F{m}.name{n},'.hdf5','.fig')], 'fig');
            close(H_lines)
            end
            
        end
        
        H.s_Cp = figure; hold on;
        H.limit_line = figure; hold on;
        ts_plot_surface(g,[0 0 0],H.limit_line);
        F{m} = surface_hub_limit(directory,pv_directory,F{m},H); 

%         F{m} = hub_mid_cas(pv_directory,F{m},H);
%         if plot_2D ~= -1
%             H.wakes_2D = figure(); hold on;
%             F{m} = slice_2D(pv_directory,F{m},H); 
%         end
    end
    end
else
    if ishandle(H.break)
        close(H.break)
    end
    if ishandle(H.break_2)
        close(H.break_2) 
    end
    if ishandle(H.break_wakes)
        close(H.break_wakes)
    end
end

for m = 1:length(F)
    f = F{m};
    save([directory strrep(filestring, '*', [num2str(F{m}.M*100) '00_cuts.mat'])],'f','-v7.3')
end

%% Calculate endwall curvature and secondary flow by mass flow rate
if plot_curv == 1
    pv_directory = strrep(strrep(directory,'TURBOSTREAM','PARAVIEW'),'TS','PV');
    
    figure(); hold on;
    for m = 1:length(F)
        [~, n_des] = min(abs(F{m}.number - F{m}.M * 10000));
        
        % Seed locations
        load([directory strrep(F{m}.name{n_des},'.hdf5','.mat')])
        S.x = g{7}.x(:,1,end); S.r = g{7}.r(:,1,end); S.rt = g{7}.rt(:,1,end); 
%         CLT = CLT_Blade2CLTCirc(g{5}.x(:,1,1),g{5}.rt(:,1,1),0); iLE = CLT.surfaces.iLE;
%         xrrt_hub_LE = [g{5}.x(iLE,1,1) g{5}.r(iLE,1,1) g{5}.rt(iLE,1,1)];

        F{m}.psi_hub = nan(size(F{m}.number)); F{m}.m_sec = nan(size(F{m}.number));
        for n = n_des:length(F{m}.number)
            pv_filename = strrep(F{m}.name{n},'.hdf5','_paraview.hdf5');
            curv_name = strrep(pv_filename,'.hdf5','_hub_curvature0.csv');
            if exist([pv_directory curv_name],'file') == 0 && ...
                    exist([pv_directory strrep(curv_name,'csv','mat')],'file') == 0
                [status, result] = system(['pvpython /home/jvt24/Documents/PARAVIEW/ts_hub_curv.py '...
                '"' pv_directory '" "' pv_filename '"']); 
            end

            C = read_paraview([pv_directory curv_name],0,1,S);

            p = struct('x',1);
            p.x = C.paths(:,1,:); p.r = C.paths(:,2,:); p.rt = C.paths(:,3,:);
            p.ro = C.primary(:,1,:); p.rovx = C.primary(:,2,:); p.rovr = C.primary(:,3,:); 
            p.rorvt = C.primary(:,4,:); p.roe = C.primary(:,5,:); p.omega = C.omega(:,1,:);

            % Squeeze up the arrays and flip directions
            varnames = fieldnames(p);
            for v = 1:length(varnames)
                varname = varnames{v};
                p.(varname) = squeeze(p.(varname));
                for j = 1:size(p.(varname),2) % Flip paths front to back
                    if isnan(p.(varname)(end,j)) == 1
                        j_nan = find(isnan(p.(varname)(:,j)),1,'first') - 1;
                    else
                        j_nan = size(p.(varname),1);
                    end
                    p.(varname)(1:j_nan,j) = p.(varname)(j_nan:-1:1,j);
                end
            end
            p = ts_secondary(p);

            % Plot lines coloured by angular velocity
%             ts_plot_surface(g,[0 0 0])
%             for i = 1:size(p.x,2)
%                 surface([p.x(:,i).';p.x(:,i).'],[p.rt(:,i).';p.rt(:,i).'],...
%                     [p.r(:,i).';p.r(:,i).'],[p.omega(:,i).';p.omega(:,i).'],...
%                     'facecol','no','edgecol','interp','linew',2);
%             end
            
            % Calculate angular velocities at periodic boundary
%             omega = zeros(size(S.x)); ro = zeros(size(S.x)); V = zeros(size(S.x));
%             rt = zeros(size(S.x)); r = zeros(size(S.x)); x = zeros(size(S.x));
%             for i = 1:length(S.x)
%                 d = ((p.x(:,i) - S.x(i)).^2 + (p.r(:,i) - S.r(i)).^2 + (p.rt(:,i) - S.rt(i)).^2).^0.5;
%                 [~, j] = min(abs(d)); 
%                 omega(i) = p.omega(j,i); ro(i) = p.ro(j,i); V(i) = p.V(j,i);
%                 x(i) = p.x(j,i); r(i) = p.r(j,i); rt(i) = p.rt(j,i);
%             end
%             Rad = V ./ omega;
%             plot3(x,rt,r,'k.')
%             plot(x,Rad,'-','Color',[n/length(F{m}.number) 0 0])
%             F{m}.radius(n) = mean(R(x > 0.182 & x < 0.195));

            % Calculate angular velocity along trailing edge streamline
%             i = 1;
%             while min(p.x(:,i)) > 0.18 && max(p.x(:,i)) < 0.195
%                 i = i+1;
%             end
%             figure()
%             hold on
%             plot3(p.x(:,i),p.rt(:,i),p.r(:,i),'r-')
%             
%             varnames = fieldnames(p);
%             for v = 1:length(varnames)
%                 varname = varnames{v};
%                 l.(varname) = smooth(p.(varname)(:,i),5,'rlowess');
%             end
%             l = ts_secondary(l);
%             plot3(l.x,l.rt,l.r,'k-')
%             
%             figure()
%             l.Rad = l.V ./ l.omega;
%             plot(l.x(:,i),l.Rad(:,i))
            
            % Calculate angle in mid-passage
            psi = zeros(size(S.x));
            rt = zeros(size(S.x)); r = zeros(size(S.x)); x = zeros(size(S.x));
            for i = 3:length(S.x)-3
                d = ((p.x(:,i) - S.x(i)).^2 + (p.r(:,i) - S.r(i)).^2 + (p.rt(:,i) - S.rt(i)).^2).^0.5;
                [~, j] = min(abs(d));
                if j > 3
                    psi(i) = 360 * atan2( diff(p.rt([j-3 j+3],i)), diff(p.x([j-3 j+3],i)) ) / (2*pi); 
                else
                    psi(i) = nan;
                end
                ro(i) = p.ro(j,i); V(i) = p.V(j,i);
                x(i) = p.x(j,i); r(i) = p.r(j,i); rt(i) = p.rt(j,i);
            end
%             plot(x,psi,'-','Color',[n/length(F{m}.number) 0 0])
%             F{m}.psi(n) = sum(psi(x > 0.182 & x < 0.195)) / sum(x > 0.182 & x < 0.195);
            F{m}.psi_hub(n) = min(psi(x > 0.176 & x < 0.198));

            % Mass flow through surface offset from suction surface
            Inlet = F{m}.Inlet{n}; [~,m_1] = ts_mass_average(Inlet,'V',3);
            load([directory strrep(F{m}.name{n},'.hdf5','.mat')])
%             SS = ts_secondary(ts_structured_cut(g,4,95,165,1,40,10,10));
            SS = ts_secondary(ts_structured_cut(g,4,95,165,10,10,1,15));
            
            rovx_av = 0.25 * (SS.rovx(1:end-1,1:end-1) + SS.rovx(2:end,1:end-1) + ...
                SS.rovx(1:end-1,2:end) + SS.rovx(2:end,2:end));
            rovr_av = 0.25 * (SS.rovr(1:end-1,1:end-1) + SS.rovr(2:end,1:end-1) + ...
                SS.rovr(1:end-1,2:end) + SS.rovr(2:end,2:end));
            rpm = SS.bv.rpm;
            omega = rpm * 2 * pi() / 60;
            SS.rovt = SS.rorvt ./ SS.r - (SS.r .* SS.ro * omega);
            rovt_av = 0.25 * (SS.rovt(1:end-1,1:end-1) + SS.rovt(2:end,1:end-1) + ...
                SS.rovt(1:end-1,2:end) + SS.rovt(2:end,2:end));

            % Calculate areas
            SS.t = SS.rt ./ SS.r;
            t_av = 0.25 * (SS.t(2:end,2:end) + SS.t(1:end-1,1:end-1) + SS.t(1:end-1,2:end) + SS.t(2:end,1:end-1));
            x1 = SS.x(2:end,2:end) - SS.x(1:end-1,1:end-1);
            x2 = SS.x(1:end-1,2:end) - SS.x(2:end,1:end-1);
            r1 = SS.r(2:end,2:end) - SS.r(1:end-1,1:end-1);
            r2 = SS.r(1:end-1,2:end) - SS.r(2:end,1:end-1);
            rt1 = SS.r(2:end,2:end) .* (SS.t(2:end,2:end) - t_av)...
                - SS.r(1:end-1,1:end-1) .* (SS.t(1:end-1,1:end-1) - t_av);
            rt2 = SS.r(1:end-1,2:end) .* (SS.t(1:end-1,2:end) - t_av)...
                - SS.r(2:end,1:end-1) .* (SS.t(2:end,1:end-1) - t_av);

            Ax = 0.5 * (r1.*rt2 - r2.*rt1);
            Ar = 0.5 * (x2.*rt1 - x1.*rt2);
            At = 0.5 * (x1.*r2 - x2.*r1);

            % Find mass flux
            dM = Ax.*rovx_av + Ar.*rovr_av + At.*rovt_av;
            dM(dM > 0) = 0;
            
            F{m}.m_sec(n) = -sum(sum(dM));
        end
        
        plot(F{m}.Inc,-F{m}.psi_hub,'.-','Color',F{m}.colour);
    end
    
    
end

for m = 1:length(F)
    f = F{m};
    save([directory strrep(filestring, '*', [num2str(F{m}.M*100) '00_cuts.mat'])],'f','-v7.3')
end

end

function f = slice_2D(pv_directory,f,H)
% Read in streamlines delimiting hub mid and casing regions, integrate up
% loss to get a broken down loss loop

n_step = 2;         plot_stuff = 0;
figure(H.wakes_2D)
n_tot = length(1:n_step:length(f.number));
wi = 4; he = ceil(n_tot / wi); 
h = zeros(n_tot,1);
for n = 1:n_tot
    h(n) = subplot(he,wi,n);
end


Zeta = zeros(size(f.number)); w = 0; L = cell(n_tot,1);
for n = 1:n_step:length(f.number)
    w = w+1;
    pv_filename = strrep(f.name{n},'.hdf5','_paraview.hdf5');
    Inlet = f.Inlet{n}; %Outlet = f.Outlet{n};
    Outlet = f.cuts{n,2};
    if exist([pv_directory strrep(pv_filename,'.hdf5','_2D_cut_loss.mat')],'file') == 0
        low_coarse = read_paraview([pv_directory strrep(pv_filename,'.hdf5','_2D_low_coarse.csv')],[],1,[]);
        low_fine = read_paraview([pv_directory strrep(pv_filename,'.hdf5','_2D_low_fine.csv')],[],1,[]);
        high_coarse = read_paraview([pv_directory strrep(pv_filename,'.hdf5','_2D_high_coarse.csv')],[],1,[]);
        high_fine = read_paraview([pv_directory strrep(pv_filename,'.hdf5','_2D_high_fine.csv')],[],1,[]);
        
%         ts_plot_surface(g,[0 0 0])
%         for i = 1:size(hub_coarse,3)
%             plot3(hub_coarse(:,1,i),hub_coarse(:,3,i),hub_coarse(:,2,i),'r-')
%         end
%         for i = 1:size(hub_fine,3)
%             plot3(hub_fine(:,1,i),hub_fine(:,3,i),hub_fine(:,2,i),'b-')
%         end
%         for i = 1:size(cas_coarse,3)
%             plot3(cas_coarse(:,1,i),cas_coarse(:,3,i),cas_coarse(:,2,i),'r-')
%         end
%         for i = 1:size(cas_fine,3)
%             plot3(cas_fine(:,1,i),cas_fine(:,3,i),cas_fine(:,2,i),'b-')
%         end
        
        % Assemble and re order
        low = nan(max(size(low_fine,1),size(low_coarse,1)),3,size(low_fine,3) + size(low_coarse,3));
        high = nan(max(size(high_fine,1),size(high_coarse,1)),3,size(high_fine,3) + size(high_coarse,3));
        
        low(1:size(low_fine,1),:,1:size(low_fine,3)) = low_fine;
        low(1:size(low_coarse,1),:,size(low_fine,3)+1:size(low_coarse,3)+size(low_fine,3)) = low_coarse;
        high(1:size(high_fine,1),:,1:size(high_fine,3)) = high_fine;
        high(1:size(high_coarse,1),:,size(high_fine,3)+1:size(high_coarse,3)+size(high_fine,3)) = high_coarse;
        
        [~, i] = sort(low(1,3,:)); low = low(:,:,i);
        [~, i] = sort(high(1,3,:)); high = high(:,:,i);
        
        % Intersect planes ensuring lines actually cross the inlet plane
        d_low = diff(low(1:2,:,:),1); d_high = diff(high(1:2,:,:),1);
        low = cat(1,low(1,:,:) - d_low,low);
        high = cat(1,high(1,:,:) - d_high,high);
%         squeeze(diff(hub(1:10,1,1:2),1))
        
        low_in = paths_intersect_plane(permute(low,[3 2 1]),Inlet);
        high_in = paths_intersect_plane(permute(high,[3 2 1]),Inlet);
        low_out = paths_intersect_plane(permute(low,[3 2 1]),Outlet);
        high_out = paths_intersect_plane(permute(high,[3 2 1]),Outlet);
        
        % Take a single pitch
        % Find middle of inlet or outlet cut and add half a pitch on each side
        [~,J_low_in] = min(abs(mean(low_in(:,2)) - Inlet.r(:,1)));
        [~,J_high_in] = min(abs(mean(high_in(:,2)) - Inlet.r(:,1)));
        [~,J_low_out] = min(abs(mean(low_out(:,2)) - Outlet.r(:,1)));
        [~,J_high_out] = min(abs(mean(high_out(:,2)) - Outlet.r(:,1)));
        t_low_in = ( (Inlet.rt(J_low_in,1) / Inlet.r(J_low_in,1)) +...
            (Inlet.rt(J_low_in,end) / Inlet.r(J_low_in,end)) ) / 2;
        t_high_in = ( (Inlet.rt(J_high_in,1) / Inlet.r(J_high_in,1)) +...
            (Inlet.rt(J_high_in,end) / Inlet.r(J_high_in,end)) ) / 2;
        t_low_out = ( (Outlet.rt(J_low_out,1) / Outlet.r(J_low_out,1)) +...
            (Outlet.rt(J_low_out,end) / Outlet.r(J_low_out,end)) ) / 2;
        t_high_out = ( (Inlet.rt(J_high_out,1) / Inlet.r(J_high_out,1)) +...
            (Outlet.rt(J_high_out,end) / Outlet.r(J_high_out,end)) ) / 2;
        
        nblade = double(Inlet.bv.nblade);
%         d = (diff(hub_in(:,1)).^2 + diff(hub_in(:,2)).^2 + diff(hub_in(:,3)).^2).^0.5;
%         [~, i] = min(d); t = hub_in(i,3) / hub_in(i,2);
        t_max = t_low_in + pi / nblade; t_min = t_low_in - pi / nblade; 
        t = low_in(:,3) ./ low_in(:,2); low_in(t > t_max | t < t_min,:) = [];
%         
%         d = (diff(cas_in(:,1)).^2 + diff(cas_in(:,2)).^2 + diff(cas_in(:,3)).^2).^0.5;
%         [~, i] = min(d); t = cas_in(i,3) / cas_in(i,2);
        t_max = t_high_in + pi / nblade; t_min = t_high_in - pi / nblade; 
        t = high_in(:,3) ./ high_in(:,2); high_in(t > t_max | t < t_min,:) = [];
%         
%         d = (diff(hub_out(:,1)).^2 + diff(hub_out(:,2)).^2 + diff(hub_out(:,3)).^2).^0.5;
%         [~, i] = min(d); t = hub_out(i,3) / hub_out(i,2);
        t_max = t_low_out + pi / nblade; t_min = t_low_out - pi / nblade; 
        t = low_out(:,3) ./ low_out(:,2); low_out(t > t_max | t < t_min,:) = [];
%         
%         d = (diff(cas_out(:,1)).^2 + diff(cas_out(:,2)).^2 + diff(cas_out(:,3)).^2).^0.5;
%         [~, i] = min(d); t = cas_out(i,3) / cas_out(i,2);
        t_max = t_high_out + pi / nblade; t_min = t_high_out - pi / nblade; 
        t = high_out(:,3) ./ high_out(:,2); high_out(t > t_max | t < t_min,:) = [];        
        
        % Repeat +/- a pitch
        low_in = [low_in(:,1:2) low_in(:,3) - low_in(:,2) * pi * 2 / nblade ; low_in ; ...
            low_in(:,1:2) low_in(:,3) + low_in(:,2) * pi * 2 / nblade];
        high_in = [high_in(:,1:2) high_in(:,3) - high_in(:,2) * pi * 2 / nblade ; high_in ; ...
            high_in(:,1:2) high_in(:,3) + high_in(:,2) * pi * 2 / nblade];
        low_out = [low_out(:,1:2) low_out(:,3) - low_out(:,2) * pi * 2 / nblade ; low_out ; ...
            low_out(:,1:2) low_out(:,3) + low_out(:,2) * pi * 2 / nblade];
        high_out = [high_out(:,1:2) high_out(:,3) - high_out(:,2) * pi * 2 / nblade ; high_out ; ...
            high_out(:,1:2) high_out(:,3) + high_out(:,2) * pi * 2 / nblade];
        
        % Reduce number of points
        d = 0; tol = 1e-5;
        while min(d) < tol
            d = (diff(low_in(:,1)).^2 + diff(low_in(:,2)).^2 + diff(low_in(:,3)).^2).^0.5;
            j = find(d < tol);
            low_in(j(1:2:end),:) = [];
        end

        d = 0;
        while min(d) < tol
            d = (diff(high_in(:,1)).^2 + diff(high_in(:,2)).^2 + diff(high_in(:,3)).^2).^0.5;
            j = find(d < tol);
            high_in(j(1:2:end),:) = [];
        end

        d = 0;
        while min(d) < tol
            d = (diff(low_out(:,1)).^2 + diff(low_out(:,2)).^2 + diff(low_out(:,3)).^2).^0.5;
            j = find(d < tol);
            low_out(j(1:2:end),:) = [];
        end

        d = 0;
        while min(d) < tol
            d = (diff(high_out(:,1)).^2 + diff(high_out(:,2)).^2 + diff(high_out(:,3)).^2).^0.5;
            j = find(d < tol);
            high_out(j(1:2:end),:) = [];
        end

        
        % Re order
%         [~, i] = sort(hub_in(:,3)); hub_in = hub_in(i,:);
%         [~, i] = sort(cas_in(:,3)); cas_in = cas_in(i,:);
%         [~, i] = sort(hub_out(:,3)); hub_out = hub_out(i,:);
%         [~, i] = sort(cas_out(:,3)); cas_out = cas_out(i,:);
        
        mid_in = [high_in ; flipud(low_in)];
        mid_out = [high_out ; flipud(low_out)];
        
        % Mass Averages
        [Cut.s2 Cut.mass_out Cut.area_out] = ts_mass_average_nonuniform(Outlet, 's', mid_out, 2,...
            plot_stuff);
        Cut.T2 = ts_mass_average_nonuniform(Outlet, 'T', mid_out, 2, 0);    

%         [Cut.s1 Cut.mass_in Cut.area_in] = ts_mass_average_nonuniform(Inlet, 's', mid_in, 2,...
%             0);
%         Cut.V1 = ts_mass_average_nonuniform(Inlet, 'V', mid_in, 2, 0);    
        Inlet = ts_secondary(Inlet);
        r_in = mean(mid_in(:,2)); r_1 = Inlet.r_nondim(:,1);
        Cut.s1 = interp1(r_1,ts_mass_average(Inlet,'s',2),r_in);
        Cut.V1 = interp1(r_1,ts_mass_average(Inlet,'V',2),r_in);
        
        [~, m_1, A_1] = ts_mass_average(Inlet,'ro',3);
        
        % Loss Coefficients
        Cut.Loss = 2 * (Cut.s2 - Cut.s1) .* Cut.T2 ./ (Cut.V1.^2);
        l.low_in = low_in; l.low_out = low_out; l.high_in = high_in; l.high_out = high_out;
        save([pv_directory strrep(pv_filename,'.hdf5','_2D_cut_loss.mat')],'Cut','m_1','l','-v7.3');
    else
        load([pv_directory strrep(pv_filename,'.hdf5','_2D_cut_loss.mat')])
    end
    
    % Plot wakes with streamlines
    l.hub_in = l.low_in; l.cas_in = l.high_in; l.hub_out = l.low_out; l.cas_out = l.high_out; 
    plot_break_wakes_2(Inlet,Outlet,l,1,h(w))
    L{w} = l;
    Zeta(w) = Cut.Loss;
end

% Plot new slice of loss on the original mises, 2D slice plot
Inc = f.Inc(Zeta ~= 0);
Zeta = Zeta(Zeta ~= 0);

figure(H.mid)
plot(Inc,Zeta,'-.','Color',f.colour)

f.Loss_2D_cut = Zeta;
f.Inc_2D_cut = Inc;
f.Lines = L;


end



function f = hub_mid_cas(pv_directory,f,H)
% Read in streamlines delimiting hub mid and casing regions, integrate up
% loss to get a broken down loss loop

n_step = 1;
figure(H.break_wakes)
n_tot = length(1:n_step:length(f.number));
wi = 4; he = ceil(n_tot / wi); 
h = zeros(n_tot,1);
for n = 1:n_tot
    h(n) = subplot(he,wi,n);
end

Zeta = cell(size(f.number)); w = 0; L = cell(n_tot,1);
for n = 1:n_step:length(f.number)
    w = w+1;
    pv_filename = strrep(f.name{n},'.hdf5','_paraview.hdf5');
    Inlet = f.Inlet{n}; %Outlet = f.Outlet{n};
    Outlet = f.cuts{n,2};
    if exist([pv_directory strrep(pv_filename,'.hdf5','_hub_mid_cas_loss.mat')],'file') == 0
        hub_coarse = read_paraview([pv_directory strrep(pv_filename,'.hdf5','_hub_bl_coarse.csv')],[],1,[]);
        hub_fine = read_paraview([pv_directory strrep(pv_filename,'.hdf5','_hub_bl_fine.csv')],[],1,[]);
        cas_coarse = read_paraview([pv_directory strrep(pv_filename,'.hdf5','_cas_bl_coarse.csv')],[],1,[]);
        cas_fine = read_paraview([pv_directory strrep(pv_filename,'.hdf5','_cas_bl_fine.csv')],[],1,[]);
        
%         ts_plot_surface(g,[0 0 0])
%         for i = 1:size(hub_coarse,3)
%             plot3(hub_coarse(:,1,i),hub_coarse(:,3,i),hub_coarse(:,2,i),'r-')
%         end
%         for i = 1:size(hub_fine,3)
%             plot3(hub_fine(:,1,i),hub_fine(:,3,i),hub_fine(:,2,i),'b-')
%         end
%         for i = 1:size(cas_coarse,3)
%             plot3(cas_coarse(:,1,i),cas_coarse(:,3,i),cas_coarse(:,2,i),'r-')
%         end
%         for i = 1:size(cas_fine,3)
%             plot3(cas_fine(:,1,i),cas_fine(:,3,i),cas_fine(:,2,i),'b-')
%         end
        
        % Assemble and re order
        hub = nan(max(size(hub_fine,1),size(hub_coarse,1)),3,size(hub_fine,3) + size(hub_coarse,3));
        cas = nan(max(size(cas_fine,1),size(cas_coarse,1)),3,size(cas_fine,3) + size(cas_coarse,3));
        
        hub(1:size(hub_fine,1),:,1:size(hub_fine,3)) = hub_fine;
        hub(1:size(hub_coarse,1),:,size(hub_fine,3)+1:size(hub_coarse,3)+size(hub_fine,3)) = hub_coarse;
        cas(1:size(cas_fine,1),:,1:size(cas_fine,3)) = cas_fine;
        cas(1:size(cas_coarse,1),:,size(cas_fine,3)+1:size(cas_coarse,3)+size(cas_fine,3)) = cas_coarse;
        
        [~, i] = sort(hub(1,3,:)); hub = hub(:,:,i);
        [~, i] = sort(cas(1,3,:)); cas = cas(:,:,i);
        
        % Intersect planes ensuring lines actually cross the inlet plane
        d_hub = diff(hub(1:2,:,:),1); d_cas = diff(cas(1:2,:,:),1);
        hub = cat(1,hub(1,:,:) - d_hub,hub);
        cas = cat(1,cas(1,:,:) - d_cas,cas);
%         squeeze(diff(hub(1:10,1,1:2),1))
        
        hub_in = paths_intersect_plane(permute(hub,[3 2 1]),Inlet);
        cas_in = paths_intersect_plane(permute(cas,[3 2 1]),Inlet);
        hub_out = paths_intersect_plane(permute(hub,[3 2 1]),Outlet);
        cas_out = paths_intersect_plane(permute(cas,[3 2 1]),Outlet);
        
        % Take a single pitch
        % Find middle of inlet or outlet cut and add half a pitch on each side
        [~,J_hub_in] = min(abs(mean(hub_in(:,2)) - Inlet.r(:,1)));
        [~,J_cas_in] = min(abs(mean(cas_in(:,2)) - Inlet.r(:,1)));
        [~,J_hub_out] = min(abs(mean(hub_out(:,2)) - Outlet.r(:,1)));
        [~,J_cas_out] = min(abs(mean(cas_out(:,2)) - Outlet.r(:,1)));
        t_hub_in = ( (Inlet.rt(J_hub_in,1) / Inlet.r(J_hub_in,1)) +...
            (Inlet.rt(J_hub_in,end) / Inlet.r(J_hub_in,end)) ) / 2;
        t_cas_in = ( (Inlet.rt(J_cas_in,1) / Inlet.r(J_cas_in,1)) +...
            (Inlet.rt(J_cas_in,end) / Inlet.r(J_cas_in,end)) ) / 2;
        t_hub_out = ( (Outlet.rt(J_hub_out,1) / Outlet.r(J_hub_out,1)) +...
            (Outlet.rt(J_hub_out,end) / Outlet.r(J_hub_out,end)) ) / 2;
        t_cas_out = ( (Inlet.rt(J_cas_out,1) / Inlet.r(J_cas_out,1)) +...
            (Outlet.rt(J_cas_out,end) / Outlet.r(J_cas_out,end)) ) / 2;
        
        nblade = double(Inlet.bv.nblade);
%         d = (diff(hub_in(:,1)).^2 + diff(hub_in(:,2)).^2 + diff(hub_in(:,3)).^2).^0.5;
%         [~, i] = min(d); t = hub_in(i,3) / hub_in(i,2);
        t_max = t_hub_in + pi / nblade; t_min = t_hub_in - pi / nblade; 
        t = hub_in(:,3) ./ hub_in(:,2); hub_in(t > t_max | t < t_min,:) = [];
%         
%         d = (diff(cas_in(:,1)).^2 + diff(cas_in(:,2)).^2 + diff(cas_in(:,3)).^2).^0.5;
%         [~, i] = min(d); t = cas_in(i,3) / cas_in(i,2);
        t_max = t_cas_in + pi / nblade; t_min = t_cas_in - pi / nblade; 
        t = cas_in(:,3) ./ cas_in(:,2); cas_in(t > t_max | t < t_min,:) = [];
%         
%         d = (diff(hub_out(:,1)).^2 + diff(hub_out(:,2)).^2 + diff(hub_out(:,3)).^2).^0.5;
%         [~, i] = min(d); t = hub_out(i,3) / hub_out(i,2);
        t_max = t_hub_out + pi / nblade; t_min = t_hub_out - pi / nblade; 
        t = hub_out(:,3) ./ hub_out(:,2); hub_out(t > t_max | t < t_min,:) = [];
%         
%         d = (diff(cas_out(:,1)).^2 + diff(cas_out(:,2)).^2 + diff(cas_out(:,3)).^2).^0.5;
%         [~, i] = min(d); t = cas_out(i,3) / cas_out(i,2);
        t_max = t_cas_out + pi / nblade; t_min = t_cas_out - pi / nblade; 
        t = cas_out(:,3) ./ cas_out(:,2); cas_out(t > t_max | t < t_min,:) = [];        
        
        % Repeat +/- a pitch
        hub_in = [hub_in(:,1:2) hub_in(:,3) - hub_in(:,2) * pi * 2 / nblade ; hub_in ; ...
            hub_in(:,1:2) hub_in(:,3) + hub_in(:,2) * pi * 2 / nblade];
        cas_in = [cas_in(:,1:2) cas_in(:,3) - cas_in(:,2) * pi * 2 / nblade ; cas_in ; ...
            cas_in(:,1:2) cas_in(:,3) + cas_in(:,2) * pi * 2 / nblade];
        hub_out = [hub_out(:,1:2) hub_out(:,3) - hub_out(:,2) * pi * 2 / nblade ; hub_out ; ...
            hub_out(:,1:2) hub_out(:,3) + hub_out(:,2) * pi * 2 / nblade];
        cas_out = [cas_out(:,1:2) cas_out(:,3) - cas_out(:,2) * pi * 2 / nblade ; cas_out ; ...
            cas_out(:,1:2) cas_out(:,3) + cas_out(:,2) * pi * 2 / nblade];
        
        % Reduce number of points
        d = 0; tol = 1e-5;
        while min(d) < tol
            d = (diff(hub_in(:,1)).^2 + diff(hub_in(:,2)).^2 + diff(hub_in(:,3)).^2).^0.5;
            j = find(d < tol);
            hub_in(j(1:2:end),:) = [];
        end

        d = 0;
        while min(d) < tol
            d = (diff(cas_in(:,1)).^2 + diff(cas_in(:,2)).^2 + diff(cas_in(:,3)).^2).^0.5;
            j = find(d < tol);
            cas_in(j(1:2:end),:) = [];
        end

        d = 0;
        while min(d) < tol
            d = (diff(hub_out(:,1)).^2 + diff(hub_out(:,2)).^2 + diff(hub_out(:,3)).^2).^0.5;
            j = find(d < tol);
            hub_out(j(1:2:end),:) = [];
        end

        d = 0;
        while min(d) < tol
            d = (diff(cas_out(:,1)).^2 + diff(cas_out(:,2)).^2 + diff(cas_out(:,3)).^2).^0.5;
            j = find(d < tol);
            cas_out(j(1:2:end),:) = [];
        end

        
        % Re order
%         [~, i] = sort(hub_in(:,3)); hub_in = hub_in(i,:);
%         [~, i] = sort(cas_in(:,3)); cas_in = cas_in(i,:);
%         [~, i] = sort(hub_out(:,3)); hub_out = hub_out(i,:);
%         [~, i] = sort(cas_out(:,3)); cas_out = cas_out(i,:);
        
        mid_in = [cas_in ; flipud(hub_in)];
        mid_out = [cas_out ; flipud(hub_out)];
        
        % Mass Averages
        plot_stuff = 0;
        [Hub.s2 Hub.mass_out Hub.area_out] = ts_mass_average_nonuniform(Outlet, 's', hub_out, 1,...
            plot_stuff);
        Hub.T2 = ts_mass_average_nonuniform(Outlet, 'T', hub_out, 1, plot_stuff);

        [Mid.s2 Mid.mass_out Mid.area_out] = ts_mass_average_nonuniform(Outlet, 's', mid_out, 2,...
            plot_stuff);
        Mid.T2 = ts_mass_average_nonuniform(Outlet, 'T', mid_out, 2, plot_stuff);    

        [Cas.s2 Cas.mass_out Cas.area_out] = ts_mass_average_nonuniform(Outlet, 's', cas_out, 3,...
            plot_stuff);
        Cas.T2 = ts_mass_average_nonuniform(Outlet, 'T', cas_out, 3, plot_stuff);
        
        [Hub.s1 Hub.mass_in Hub.area_in] = ts_mass_average_nonuniform(Inlet, 's', hub_in, 1,...
            plot_stuff);
        Hub.V1 = ts_mass_average_nonuniform(Inlet, 'V', hub_in, 1, plot_stuff);

        [Mid.s1 Mid.mass_in Mid.area_in] = ts_mass_average_nonuniform(Inlet, 's', mid_in, 2,...
            plot_stuff);
        Mid.V1 = ts_mass_average_nonuniform(Inlet, 'V', mid_in, 2, plot_stuff);    

        [Cas.s1 Cas.mass_in Cas.area_in] = ts_mass_average_nonuniform(Inlet, 's', cas_in, 3,...
            plot_stuff);
        Cas.V1 = ts_mass_average_nonuniform(Inlet, 'V', cas_in, 3, plot_stuff);
        
        [~, m_1, A_1] = ts_mass_average(Inlet,'ro',3);
        [~, m_2, A_2] = ts_mass_average(Outlet,'ro',3);
        
        % Loss Coefficients
        Zeta{n}.hub = (Hub.mass_in/m_1) * 2 * (Hub.s2 - Hub.s1) .* Hub.T2 ./ (Hub.V1.^2);
        Zeta{n}.mid = (Mid.mass_in/m_1) * 2 * (Mid.s2 - Mid.s1) .* Mid.T2 ./ (Mid.V1.^2);
        Zeta{n}.cas = (Cas.mass_in/m_1) * 2 * (Cas.s2 - Cas.s1) .* Cas.T2 ./ (Cas.V1.^2);
        Loss.hub = Zeta{n}.hub;  Loss.mid = Zeta{n}.mid;  Loss.cas = Zeta{n}.cas; 
        l.hub_in = hub_in; l.hub_out = hub_out; l.cas_in = cas_in; l.cas_out = cas_out;
        save([pv_directory strrep(pv_filename,'.hdf5','_hub_mid_cas_loss.mat')],'Hub','Mid',...
            'Cas','m_1','l','-v7.3');
    else
        load([pv_directory strrep(pv_filename,'.hdf5','_hub_mid_cas_loss.mat')])
        Zeta{n}.hub = (Hub.mass_in/m_1) * 2 * (Hub.s2 - Hub.s1) .* Hub.T2 ./ (Hub.V1.^2);
        Zeta{n}.mid = (Mid.mass_in/m_1) * 2 * (Mid.s2 - Mid.s1) .* Mid.T2 ./ (Mid.V1.^2);
        Zeta{n}.cas = (Cas.mass_in/m_1) * 2 * (Cas.s2 - Cas.s1) .* Cas.T2 ./ (Cas.V1.^2);   
    end
    
    % Plot wakes with streamlines
    plot_break_wakes_2(Inlet,Outlet,l,1,h(w))
    L{w} = l;
    
end

% Total loss
Total_Loss = zeros(size(f.number));
for n = 1:length(f.number)
    Inlet = f.Inlet{n}; %Outlet = f.Outlet{n};
    Outlet = f.cuts{n,2};
    
    T_2 = ts_mass_average(Outlet,'T',3);
    s_2 = ts_mass_average(Outlet,'s',3);
    V_1 = ts_mass_average(Inlet,'V',3);
    s_1 = ts_mass_average(Inlet,'s',3);    
    
    Total_Loss(n) = T_2 * (s_2 - s_1) / (0.5 * V_1^2);
end

% Loss loop breakdowns
[~, n_des] = min(abs(f.number - f.M * 10000));
Inc_des = f.Inc(n_des); Loss_des = f.Loss(n_des);
Inc = []; Hub_Loss = []; Mid_Loss = []; Cas_Loss = []; TF_Loss = []; 
m = 1;
for n = 1:length(Zeta)
    if isempty(Zeta{n}) == 0
        Inc(m) = f.Inc(n);
        Hub_Loss(m) = Zeta{n}.hub;
        Mid_Loss(m) = Zeta{n}.mid;
        Cas_Loss(m) = Zeta{n}.cas;
        TF_Loss(m) = Total_Loss(n);
        m = m+1;
    end
end

% Unscaled
figure(H.break);
hold on
plot(f.Inc - Inc_des, Total_Loss / Loss_des, 'k--');
plot(Inc - Inc_des, Hub_Loss / Loss_des, 'b--');
plot(Inc - Inc_des, (Hub_Loss + Mid_Loss) / Loss_des, 'r--');
plot(Inc - Inc_des, (Hub_Loss + Mid_Loss + Cas_Loss) / Loss_des, 'g--');

f.Hub_Loss.Unscaled = Hub_Loss;
f.Mid_Loss.Unscaled = Mid_Loss;
f.Cas_Loss.Unscaled = Cas_Loss;

% Scaled to total loss loop
frac = TF_Loss ./ (Hub_Loss + Mid_Loss + Cas_Loss);
Hub_Loss = Hub_Loss .* frac;
Mid_Loss = Mid_Loss .* frac;
Cas_Loss = Cas_Loss .* frac;

plot(f.Inc - Inc_des, Total_Loss / Loss_des, 'k.-');
plot(Inc - Inc_des, Hub_Loss / Loss_des, 'b.-');
plot(Inc - Inc_des, (Hub_Loss + Mid_Loss) / Loss_des, 'r.-');
plot(Inc - Inc_des, (Hub_Loss + Mid_Loss + Cas_Loss) / Loss_des, 'g.-');
v = axis; axis([v(1:2) 0 v(4)])
grid on; box on;

f.Hub_Loss.Scaled = Hub_Loss;
f.Mid_Loss.Scaled = Mid_Loss;
f.Cas_Loss.Scaled = Cas_Loss;
f.Inc_break = Inc;
f.Lines = L;

end


function f = surface_hub_limit(directory,pv_directory,f,H)
% Read in surface streamlines and pick most severe one limiting hub corner
% stall

load([directory strrep(f.name{1},'.hdf5','.mat')])
S.x = g{5}.x(91,1:11,1:25); S.r = g{5}.r(91,1:15,1:25); S.rt = g{5}.rt(91,1:11,1:25);
CLT = CLT_Blade2CLTCirc(g{5}.x(:,1,1),g{5}.rt(:,1,1),0); iLE = CLT.surfaces.iLE;
xrrt_hub_LE = [g{5}.x(iLE,1,1) g{5}.r(iLE,1,1) g{5}.rt(iLE,1,1)];
% x_1 = 0.1726; x_2 = 0.2; 
x_mid = 0.18;

[~,i_des] = min(abs(f.number - f.M * 10000));
numbers = f.number(i_des:end);

figure(H.dCp_limit);
subplot(2,1,1); hold on
subplot(2,1,2); hold on

dCp_max = nan(size(numbers)); inc_dCp_max = nan(size(numbers)); 
Cp_max = nan(size(numbers)); inc_Cp_max = nan(size(numbers));

figure(H.s_Cp)
hold on
%     subplot(2,1,1); hold on
%     subplot(2,1,2); hold on

l = cell(length(numbers),1);

for n = 1:length(numbers)
    pv_filename = strrep(f.name{f.number == numbers(n)},'.hdf5','_paraview.hdf5');
    surfline_name = strrep(pv_filename,'.hdf5','_hub_surflines0.csv');
    
    C = read_paraview([pv_directory surfline_name],0,1,S);
    
    p = struct('x',1);
    p.x = C.paths(:,1,:); p.r = C.paths(:,2,:); p.rt = C.paths(:,3,:);
    p.ro = C.primary(:,1,:); p.rovx = C.primary(:,2,:); p.rovr = C.primary(:,3,:); 
    p.rorvt = C.primary(:,4,:); p.roe = C.primary(:,5,:);
%         p.su_n_x = C.surf_normals(:,1,:); p.su_n_r = C.surf_normals(:,2,:);
%         p.su_n_rt = C.surf_normals(:,3,:); 
%         p.st_n_x = C.stream_normals(:,1,:); p.st_n_r = C.stream_normals(:,2,:);
%         p.st_n_rt = C.stream_normals(:,3,:); 
%         p.dPdx = C.dP(:,1,:); p.dPdr = C.dP(:,2,:); p.dPdrt = C.dP(:,3,:);
%         p.omega = C.omega;

    % Squeeze up the arrays and flip directions
    varnames = fieldnames(p);
    for v = 1:length(varnames)
        varname = varnames{v};
        p.(varname) = squeeze(p.(varname));
        for j = 1:size(p.(varname),2) % Flip paths front to back
            if isnan(p.(varname)(end,j)) == 1
                j_nan = find(isnan(p.(varname)(:,j)),1,'first') - 1;
            else
                j_nan = size(p.(varname),1);
            end
            p.(varname)(1:j_nan,j) = p.(varname)(j_nan:-1:1,j);
        end
    end
    
    % Truncate all paths starting from the point closest to the leading edge
    % on the fillet
    for j = 1:size(p.x,2)
        d = ((p.x(:,j) - xrrt_hub_LE(1)).^2 + (p.r(:,j) - xrrt_hub_LE(2)).^2 ...
        + (p.rt(:,j) - xrrt_hub_LE(3)).^2).^0.5;
        [~, i] = min(abs(d));
        for v = 1:length(varnames)
            varname = varnames{v};
            p.(varname)(1:end-i+1,j) = p.(varname)(i:end,j);
            p.(varname)(end-i+2:end,j) = nan;
        end
    end
            
    p = ts_secondary(p);
    P_1 = ts_mass_average(f.Inlet{f.number == numbers(n)},'P',3);
    Po_1 = ts_mass_average(f.Inlet{f.number == numbers(n)},'Po',3);
    P = p.P;
    
    s = ((diff(p.x,1).^2) + (diff(p.r,1).^2) + (diff(p.rt,1).^2)).^0.5;
    s = cat(1,zeros(1,size(p.x,2)),cumsum(s,1));
%     s = repmat(max(s,[],1),[size(s,1) 1]) - s;

    Cp = (P - P_1) ./ (Po_1 - P_1);

    % Record maximum Cp difference and maximum Cp
    Cp_maxs = max(-Cp,[],1); [x_max i_x_max] = max(p.x,[],1);
    Cp_tes = nan(1,size(Cp,2));
    for i = 1:size(Cp,2); 
        if x_max(i) > x_mid
            Cp_tes(i) = -Cp(i_x_max(i),i); 
        end
    end

    dCp = Cp_maxs - Cp_tes;

    [dCp_max(n) n_dmax] = max(dCp);
    [Cp_max(n) n_max] = max(Cp_maxs);
    inc_Cp_max(n) = f.Inc(f.number == numbers(n));
    inc_dCp_max(n) = f.Inc(f.number == numbers(n));

%         subplot(2,1,1);
%         plot(p.x,-Cp)
%         plot(p.x(:,n_max),-Cp(:,n_max),'k-','LineWidth',3)
%         plot(p.x(:,n_dmax),-Cp(:,n_dmax),'r-','LineWidth',3)
%         plot(s,-Cp)
    figure(H.s_Cp)
    plot(s(:,n_dmax),-Cp(:,n_dmax),'-','LineWidth',3,'Color',[n/length(numbers) 0 0])
%         axis([0.17 0.2 -0.6 0.8]) 

    % Calculate cross stream pressure gradients from dp calculations
%         V = cat(3,p.Vx,p.Vr,p.Vt);
%         su_n = cat(3,p.su_n_x,p.su_n_r,p.su_n_rt);
%         st_n = cat(3,p.st_n_x,p.st_n_r,p.st_n_rt);
%         dP = cat(3,p.dPdx,p.dPdr,p.dPdrt);
%         q = cat(3,su_n(:,:,2) .* V(:,:,3) - su_n(:,:,3) .* V(:,:,2),...
%             su_n(:,:,3) .* V(:,:,1) - su_n(:,:,1) .* V(:,:,3),...
%             su_n(:,:,1) .* V(:,:,2) - su_n(:,:,2) .* V(:,:,1));
%         
%         q = q ./ repmat(sum(q.^2,3).^0.5,[1 1 3]);
% %         
%         dPdq = sum(q .* dP,3);
%         
%         % From radius of curvatures
%         R = p.V ./ p.omega;
%         dPdn = p.ro .* (p.V).^2 ./ R;
%         dPdq = dPdq .* sum(st_n .* q,3);
%         
%         
%         subplot(2,1,2)
% %         plot(p.x,dPdq)
%         plot(p.x(:,n_dmax),dPdq(:,n_dmax),'r-','LineWidth',3)

    l{n}.number = numbers(n);
    l{n}.x = p.x(:,n_dmax);
    l{n}.r = p.r(:,n_dmax);
    l{n}.rt = p.rt(:,n_dmax);
    l{n}.s = s(:,n_dmax);
    l{n}.Cp = Cp(:,n_dmax);
    l{n}.ro = p.ro(:,n_dmax);
    l{n}.rovx = p.rovx(:,n_dmax);
    l{n}.rovr = p.rovr(:,n_dmax);
    l{n}.rorvt = p.rorvt(:,n_dmax);
    l{n}.roe = p.roe(:,n_dmax);
    
    figure(H.limit_line)
    plot3(l{n}.x,l{n}.rt,l{n}.r,'r-')
end

figure(H.dCp_limit)
subplot(2,1,1); plot(inc_Cp_max,Cp_max,'.-','Color',f.colour);
subplot(2,1,2); plot(inc_dCp_max,dCp_max,'.-','Color',f.colour);

% Record output streamlines and pressure gradient information
f.inc_dCp_max = inc_dCp_max;
f.dCp_max = dCp_max;
f.hub_limit = l;

end

function [] = plot_break_wakes(Inlet,Outlet,Loss,Surface,plot_entropy,n,d,H)

% Plot double passage wakes onto subplot
figure(H.break_wakes)
hfig = subplot(d,5,n);
t_mean_2 = ts_plot_wake(Inlet, Outlet, hfig, plot_entropy, 2, 1);

% Plot r lines onto subplots
t_hub_min = min(Outlet.t(Loss.J.hub_out_nearest,:)); t_hub_max = max(Outlet.t(Loss.J.hub_out_nearest,:));
t_hub_max = 2*t_hub_max - t_hub_min; t = linspace(t_hub_min,t_hub_max,200) - t_mean_2;
r = ones(1,200)*Loss.outlet.r_hub;
[y z] = pol2cart(t,r);
plot(z, y, 'k--')

t_cas_min = min(Outlet.t(Loss.J.cas_out_nearest,:)); t_cas_max = max(Outlet.t(Loss.J.cas_out_nearest,:));
t_cas_max = 2*t_cas_max - t_cas_min; t = linspace(t_cas_min,t_cas_max,200) - t_mean_2;
r = ones(1,200)*Loss.outlet.r_cas;
[y z] = pol2cart(t,r);
plot(z, y, 'k--')

% Plot dividing lines on subplots
line_hub = Surface.hub.line; line_cas = Surface.cas.line;
r_line_hub = line_hub(:,2); t_line_hub = line_hub(:,3)./line_hub(:,2);
r_line_hub = r_line_hub(t_line_hub > t_hub_min & t_line_hub < t_hub_max);
t_line_hub = t_line_hub(t_line_hub > t_hub_min & t_line_hub < t_hub_max) - t_mean_2;

r_line_cas = line_cas(:,2); t_line_cas = line_cas(:,3)./line_cas(:,2);
r_line_cas = r_line_cas(t_line_cas > t_cas_min & t_line_cas < t_cas_max);
t_line_cas = t_line_cas(t_line_cas > t_cas_min & t_line_cas < t_cas_max) - t_mean_2;

[y_hub z_hub] = pol2cart(t_line_hub, r_line_hub);
[y_cas z_cas] = pol2cart(t_line_cas, r_line_cas);

plot(z_hub,y_hub,'r.-')
plot(z_cas,y_cas,'r.-')

end

function [] = plot_break_wakes_2(Inlet,Outlet,l,plot_entropy,hfig)

% Plot double passage wakes onto subplot
t_mean_2 = ts_plot_wake(Inlet, Outlet, hfig, plot_entropy, 2, 1);

% Plot r lines onto subplots
t = l.hub_in(:,3) ./ l.hub_in(:,2) - t_mean_2;
r = l.hub_in(:,2);
[y z] = pol2cart(t,r);
plot(z, y, 'k--')

t = l.cas_in(:,3) ./ l.cas_in(:,2) - t_mean_2;
r = l.cas_in(:,2);
[y z] = pol2cart(t,r);
plot(z, y, 'k--')

% Plot dividing lines on subplots
line_hub = l.hub_out; line_cas = l.cas_out;
r_line_hub = line_hub(:,2); t_line_hub = line_hub(:,3)./line_hub(:,2);
t_line_hub = t_line_hub - t_mean_2;

r_line_cas = line_cas(:,2); t_line_cas = line_cas(:,3)./line_cas(:,2);
t_line_cas = t_line_cas - t_mean_2;

[y_hub z_hub] = pol2cart(t_line_hub, r_line_hub);
[y_cas z_cas] = pol2cart(t_line_cas, r_line_cas);

plot(z_hub,y_hub,'r.-')
plot(z_cas,y_cas,'r.-')

end
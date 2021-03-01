function [p,f,C] = bl_plot_run(dr,run_name,n_des,rows,r_secs,h,col,plot_stuff)
% Plot the current design output from a automated design function

% Make sure newest modified file is always read
A_hdf5 = [dir([dr.ts run_name '_avg.hdf5']) dir([dr.ts run_name '_avg.mat'])]; 
A_cut = dir([dr.ts run_name '_avg_cuts.mat']);
if numel(A_hdf5) == 1 && numel(A_cut) == 1 
    if A_hdf5(1).datenum > A_cut.datenum
        delete([dr.ts run_name '_avg_cuts.mat'])
    end
end

% Check if a cuts file already exists
if exist([dr.ts run_name '_avg_cuts.mat'],'file') == 0

    % Read in flow file and cut rows
    g = ts_read_hdf5([dr.ts run_name '_avg.hdf5']);
    load([dr.geom run_name '_xr_hub.mat']); load([dr.geom run_name '_xr_cas.mat']);
    C = ts_cut_rows(g,xr_hub,xr_cas);
    
    % Save cuts file
    save([dr.ts run_name '_avg_cuts.mat'],'C')
    
else
    
    % Load in cuts
    load([dr.ts run_name '_avg_cuts.mat'])
    
end

% Plot past TURBOSTREAM calculations
p = cell(length(C),1); f = cell(length(C),1);
for n = n_des
    [p{n},f{n}] = ts_plot_blade(C{n},r_secs(2:end-1),h.blade{n},col,plot_stuff);
    subtitle(rows{n});
end

% Extrapolate flow parameters to the endwalls
tarnames = fieldnames(p{1}); misnames = fieldnames(f{1});
for n = n_des
    
    % Extrapolate target parameters
    for v = 1:length(tarnames)
        if length(p{n}.(tarnames{v})) == length(r_secs) - 2
            p{n}.(tarnames{v}) = interp1(r_secs(2:end-1),p{n}.(tarnames{v}),r_secs,'linear','extrap');
        end
    end
    
    % Switch figure window
    figure(h.blade{n}); 
  
    % Plot extrapolated incidence
    subplot_tight(1,8,1,[0.08,0.01]); 
    plot(p{n}.psi_stag(1:2),r_secs(1:2),'--','color',col)
    plot(p{n}.psi_stag(end-1:end),r_secs(end-1:end),'--','color',col)
    
    % Position of peak suction
    subplot_tight(1,8,2,[0.08,0.01]); 
    plot(p{n}.s_Cp_max(1:2),r_secs(1:2),'--','color',col)
    plot(p{n}.s_Cp_max(end-1:end),r_secs(end-1:end),'--','color',col)
    
    % Shape factor area
    subplot_tight(1,8,3,[0.08,0.01]); 
    plot(p{n}.aHb_1(1:2),r_secs(1:2),'--','color',col)
    plot(p{n}.aHb_1(end-1:end),r_secs(end-1:end),'--','color',col)    
    plot(p{n}.aHb_2(1:2),r_secs(1:2),'--','color',col)
    plot(p{n}.aHb_2(end-1:end),r_secs(end-1:end),'--','color',col)
    
    % Trailing edge shape factor
    subplot_tight(1,8,4,[0.08,0.01]);
    plot(p{n}.Hb_te(1:2),r_secs(1:2),'--','color',col)
    plot(p{n}.Hb_te(end-1:end),r_secs(end-1:end),'--','color',col)    
    
    % Turning
    subplot_tight(1,8,5,[0.08,0.01]); 
    plot(p{n}.Alpha(1:2),r_secs(1:2),'--','color',col)
    plot(p{n}.Alpha(end-1:end),r_secs(end-1:end),'--','color',col)
    
    % Tangential average axial velocity gradient for vorticity calculation
    Vx = ts_mass_average(C{n}.inlet,'Vx',2); r = ts_mass_average(C{n}.inlet,'r_nondim',2);
    Vx = Vx / interp1(r,Vx,0.5); dVxdr = grad_mg(r,Vx);
    
    % Find and interpolate where the vorticity crosses the threshold 
    j1 = find(dVxdr < 0.2,1,'first'); j2 = find(dVxdr > -0.2,1,'last');
    r1 = interp1(dVxdr(j1-2:j1+2),r(j1-2:j1+2),0.2,'pchip');
    r2 = interp1(dVxdr(j2-2:j2+2),r(j2-2:j2+2),-0.2,'pchip');
    
    % Record endwall boundary layer depths
    p{n}.r_bl = [r1 r2];
    
    % Extrapolate MISES flow paramters with polynomial fit
%     r = linspace(r1,r2,100);
    for v = 1:length(misnames)
        
        % Interpolate at increased resolution
%         var_fit = interp1(r_secs(2:end-1),f{n}.(misnames{v}),r,'pchip');
        
        % Linear fit and extrapolation
%         f{n}.(misnames{v}) = polyval(polyfit(r,var_fit,2),r_secs);
        f{n}.(misnames{v}) = interp1(r_secs(2:end-1),f{n}.(misnames{v}),r_secs,'linear','extrap');
    end
    
end

end


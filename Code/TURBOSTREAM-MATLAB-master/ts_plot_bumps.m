function [F, H] = ts_plot_bumps(directory,file_prefix,M,flow_type,type,inlet,outlet,bump_prefix,mag_noise)
% A function to plot the output from many flow solutions with perturbed
% surfaces
% Choose type between 'Noise', 'Linearity' or 'Surface'

if exist('inlet','var') == 0
    inlet = 'Inlet';
end
if exist('outlet','var') == 0
    outlet = 'Outlet';
end

% Read in baseline file and calculate loss coefficient
load([directory file_prefix 'M0' num2str(100*M) '_t' num2str(flow_type) '_bump_base_avg_cut.mat']);
g_base_cut = g_cut;
g_base = ts_read_hdf5([directory file_prefix 'M0' num2str(100*M) '_t'...
    num2str(flow_type) '_bump_base_avg.hdf5']);

% Calculate loss in original baseline file
V1_base = ts_mass_average(ts_secondary(g_base_cut.(inlet)),'V',3);
T2_base = ts_mass_average(ts_secondary(g_base_cut.(outlet)),'T',3);
s1_base = ts_mass_average(ts_secondary(g_base_cut.(inlet)),'s',3);
s2_base = ts_mass_average(ts_secondary(g_base_cut.(outlet)),'s',3);
Zeta_base = double(2 * (s2_base - s1_base) .* T2_base ./ (V1_base .* V1_base));

if strcmp(type,'Noise') == 1
    disp('%% Plotting noise in baseline solution %%')
    % Plot noise in loss
    
    % Read in many baseline files
    A = dir([directory file_prefix 'M0' num2str(100*M) '_t'...
        num2str(flow_type) '_bump_base_*_avg_cut.mat']);
    g = cell(length(A),1); N = zeros(length(A),1);
    for n = 1:length(A)
        load([directory A(n).name]);
        g{n} = g_cut;
        B = regexp(A(n).name,'_','split');
        N(n) = str2double(B(7));
    end
    [~,I] = sort(N);
    g = g(I);

    % Calculate loss coefficients
    V1 = zeros(size(g)); T2 = zeros(size(g)); s1 = zeros(size(g)); s2 = zeros(size(g));
    for n = 1:length(g)
        V1(n) = ts_mass_average(ts_secondary(g{n}.(inlet)),'V',3);
        T2(n) = ts_mass_average(ts_secondary(g{n}.(outlet)),'T',3);
        s1(n) = ts_mass_average(ts_secondary(g{n}.(inlet)),'s',3);
        s2(n) = ts_mass_average(ts_secondary(g{n}.(outlet)),'s',3);
    end
    Zeta_list = 2 * (s2 - s1) .* T2 ./ (V1 .* V1);
    
    Zeta_ratio = Zeta_list ./ Zeta_base;
    
    % Calculate mean and plot distribution of loss ratios
    F.sig = std(Zeta_ratio);
    F.mu = mean(Zeta_ratio);
    F.Loss = Zeta_ratio;
    
    H.noise = figure();
    plot(Zeta_ratio,'k.')
    
    
elseif strcmp(type,'Linearity') == 1
    disp('%% Plotting linearity of bump heights %%')
    % Plot linearity with increasing bump height for a number of points
    
    % Load bump locations
    load([directory file_prefix 'Linearity_Bumps.mat']);
    
    for u = 1:size(u_bumps,1)
        u_bump = u_bumps(u,:);
        A = dir([directory file_prefix 'M0' num2str(100*M) '_t' num2str(flow_type)...
            '_bump_' num2str(u_bump(1)) '_' num2str(u_bump(2)) '_' num2str(u_bump(3))...
            '_' num2str(u_bump(4)) '*_output_avg_cut.mat']);    
    
        % Read output files
        g = cell(length(A),1);
        for n = 1:length(A)
            load([directory A(n).name]);
            g{n} = g_cut;
            
            B = regexp(A(n).name,'_','split');
            g{n}.bump.Lp = str2double(B(end-4));
            g{n}.bump.dp = str2double(B(end-3));
        end
        
        % Calculate losses
        V1 = zeros(size(g)); T2 = zeros(size(g)); s1 = zeros(size(g)); s2 = zeros(size(g));
        dp = zeros(size(g)); Lp = zeros(size(g));
        for n = 1:length(g)
            V1(n) = ts_mass_average(ts_secondary(g{n}.(inlet)),'V',3);
            T2(n) = ts_mass_average(ts_secondary(g{n}.(outlet)),'T',3);
            s1(n) = ts_mass_average(ts_secondary(g{n}.(inlet)),'s',3);
            s2(n) = ts_mass_average(ts_secondary(g{n}.(outlet)),'s',3);
            dp(n) = g{n}.bump.dp;
            Lp(n) = g{n}.bump.Lp;
        end
        Zeta_list = 2 * (s2 - s1) .* T2 ./ (V1 .* V1);
    
        % Plot output
        H(u).Linearity = figure();
        h = subplot(1,2,1);
        ts_plot_surface(g_base,[0 0 0],h);
        x = g_base{u_bump(1)+1}.x(u_bump(2)+1,u_bump(3)+1,u_bump(4)+1);
        r = g_base{u_bump(1)+1}.r(u_bump(2)+1,u_bump(3)+1,u_bump(4)+1);
        rt = g_base{u_bump(1)+1}.rt(u_bump(2)+1,u_bump(3)+1,u_bump(4)+1);
        plot3(x,rt,r,'r.','MarkerSize',40);
        
        Lp_max = max(Lp); Lp_min = min(Lp);
        subplot(1,2,2); hold on;
        for n = 1:length(g)
            plot(dp(n)/r_le * 100,Zeta_list(n)/Zeta_base,'.','Color',[(Lp(n) - Lp_min)/ (Lp_max - Lp_min) 0 0])
        end
        xlabel('d_p / r_{LE} x 100','FontSize',14); ylabel('\zeta / \zeta_{des}','FontSize',14)
        
        F{u}.Loss = Zeta_list/Zeta_base;
        F{u}.dp = dp/r_le * 100;
    end

elseif strcmp(type,'Surface') == 1
    
    load([directory file_prefix bump_prefix 'Bumps.mat']);
    load([directory file_prefix 'normal.mat']);
    d = double(g_base{1}.r(1,2,1) - g_base{1}.r(1,1,1));
    
    % Get filenames of all cut files
    A = dir([directory file_prefix 'M0' num2str(100*M) '_t' num2str(flow_type)...
        '_bump_*_output_avg_cut.mat']);

    % Find Lps which there are more than 20 files of
    A = {A.name}; B = regexp(A,'_','split'); B = B.'; B = vertcat(B{:});
    i_bump = find(strcmp(B(1,:),'bump')); i_Lp = i_bump + 5;
    Lps = str2double(B(:,i_Lp));
    
    u = unique(Lps); Lps = [u histc(Lps,u)];
    Lps = Lps(Lps(:,2) > 20); 
    
    f = 1;
    % Find dps of which there are more than 20 files
    for Lp = Lps.'
        A = dir([directory file_prefix 'M0' num2str(100*M) '_t' num2str(flow_type)...
            '_bump_*' num2str(Lp) '*_output_avg_cut.mat']);
        A = {A.name}; B = regexp(A,'_','split'); B = B.'; B = vertcat(B{:});
        dps = str2double(B(:,i_Lp + 1));
        u = unique(dps); dps = [u histc(dps,u)];    
        dps = dps(dps(:,2) > 20); 
        
        % Read in files and assign into a cell array
        for dp = dps.'
            if exist([directory file_prefix bump_prefix 'M0' num2str(100*M) '_t' num2str(flow_type)...
                    '_' num2str(Lp) '_' num2str(dp) '_Inlet_Outlet.mat'],'file') == 0
            
                % Define names in bumps data set
                filenames = cell(size(u_bumps,1),1);
                for n = 1:length(filenames)
                    u_bump = u_bumps(n,:);
                    filenames{n} = [file_prefix 'M0' num2str(100*M) '_t' num2str(flow_type)...
                        '_bump_' num2str(u_bump(1)) '_' num2str(u_bump(2)) '_' num2str(u_bump(3))...
                        '_' num2str(u_bump(4)) '_' num2str(Lp) '_' num2str(dp) '_output_avg_cut.mat'];
                end

                % Load in cut files if they exist
                g = cell(size(filenames));
                for n = 1:length(filenames)
                    if exist([directory filenames{n}],'file') ~= 0
                        load([directory filenames{n}]);

                        B = regexp(filenames{n},'_','split');
                        g{n}.bump.bid = str2double(B(i_Lp - 4));
                        g{n}.bump.I = str2double(B(i_Lp - 3));
                        g{n}.bump.J = str2double(B(i_Lp - 2));
                        g{n}.bump.K = str2double(B(i_Lp - 1));
                        g{n}.bump.Lp = str2double(B(i_Lp));
                        g{n}.bump.dp = str2double(B(i_Lp + 1));

                        g{n}.inlet = g_cut.(inlet);
                        g{n}.outlet = g_cut.(outlet);
                    end
                end
            
                g = g(~cellfun(@isempty, g));
                
                save([directory file_prefix bump_prefix 'M0' num2str(100*M) '_t' num2str(flow_type)...
                    '_' num2str(Lp) '_' num2str(dp) '_Inlet_Outlet.mat'],'g','-v7.3');

            else
                % Load all inlet and outlet planes from single mat file
                load([directory file_prefix bump_prefix 'M0' num2str(100*M) '_t' num2str(flow_type)...
                    '_' num2str(Lp) '_' num2str(dp) '_Inlet_Outlet.mat'])
            end
                            
            if length(g) > 5
                % Calculate loss ratios
                V1 = zeros(size(g)); T2 = zeros(size(g)); s1 = zeros(size(g)); s2 = zeros(size(g));
                for n = 1:length(g)
                    V1(n) = ts_mass_average(ts_secondary(g{n}.inlet),'V',3);
                    T2(n) = ts_mass_average(ts_secondary(g{n}.outlet),'T',3);
                    s1(n) = ts_mass_average(ts_secondary(g{n}.inlet),'s',3);
                    s2(n) = ts_mass_average(ts_secondary(g{n}.outlet),'s',3);
                end
                Zeta_list = double(2 * (s2 - s1) .* T2 ./ (V1 .* V1));

                % Initialise loss matrices
                Zeta = cell(size(g_base));
                for bid = ts_get_bids(g_base)
                    Zeta{bid+1} = ones(size(g_base{bid+1}.x));
                end

                % Input loss ratios of calculated points and extract
                % coordinates for interpolant
                x = zeros(size(Zeta_list)); r = zeros(size(Zeta_list)); 
                rt = zeros(size(Zeta_list)); rt1 = zeros(size(Zeta_list)); rt2 = zeros(size(Zeta_list));
                I = zeros(size(Zeta_list)); J = zeros(size(Zeta_list)); K = zeros(size(Zeta_list));
                bid = zeros(size(Zeta_list)); a = zeros(size(Zeta_list,1),3);
                for n = 1:length(g)
                    bid(n) = g{n}.bump.bid; I(n) = g{n}.bump.I;
                    J(n) = g{n}.bump.J; K(n) = g{n}.bump.K; nblade = double(g_base{bid(n)+1}.bv.nblade);

                    Zeta{bid(n)+1}(I(n)+1,J(n)+1,K(n)+1) = Zeta_list(n) / Zeta_base;
                    x(n) = double(g_base{bid(n)+1}.x(I(n)+1,J(n)+1,K(n)+1));
                    r(n) = double(g_base{bid(n)+1}.r(I(n)+1,J(n)+1,K(n)+1));
                    rt(n) = double(g_base{bid(n)+1}.rt(I(n)+1,J(n)+1,K(n)+1));
                    rt1(n) = double(g_base{bid(n)+1}.rt(I(n)+1,J(n)+1,K(n)+1)) + (2*pi*r(n)) / nblade;
                    rt2(n) = double(g_base{bid(n)+1}.rt(I(n)+1,J(n)+1,K(n)+1)) - (2*pi*r(n)) / nblade;
                    
                    a(n,:) = double(reshape(normal{bid(n)+1}(I(n)+1,J(n)+1,K(n)+1,:),[1 3]));
                end
                
                % Repeat by +/- a pitch 
                X = repmat(x,[3 1]); R = repmat(r,[3 1]); RT = [rt ; rt1 ; rt2];
                Z = repmat(Zeta_list/Zeta_base,[3 1]); A = repmat(a,[3 1]);
                
                % Repeat by +/- normals
                X = repmat(X,[3 1]) + [d * A(:,1) ; zeros(size(X)) ; - d * A(:,1)];
                R = repmat(R,[3 1]) + [d * A(:,2) ; zeros(size(R)) ; - d * A(:,2)];
                RT = repmat(RT,[3 1]) + [d * A(:,3) ; zeros(size(RT)) ; - d * A(:,3)];
                Z = repmat(Z,[3 1]);
                
                F_zeta = TriScatteredInterp(X,R,RT,Z);
                d = max(1 - min(Z), max(Z) -1);
                C = [1-d 1+d];

                % Record output
                F{f}.Loss = Zeta_list / Zeta_base;
                F{f}.bid = bid;
                F{f}.I = I;
                F{f}.J = J;
                F{f}.K = K;
                F{f}.dp = dp;
                F{f}.Lp = Lp;
                F{f}.x = x; F{f}.r = r; F{f}.rt = rt;

                % Find all walls and interpolate loss coefficient
                for bid = ts_get_bids(g_base)
                    [I,J,K] = ind2sub(size(g_base{bid+1}.mwall), find(g_base{bid+1}.mwall == 0));
                    I = I-1; J = J-1; K = K-1;

                    x_interp = zeros(length(I),1); r_interp = zeros(length(I),1); 
                    rt_interp = zeros(length(I),1);
                    for n = 1:length(I)
                        x_interp(n) = g_base{bid+1}.x(I(n)+1,J(n)+1,K(n)+1);
                        r_interp(n) = g_base{bid+1}.r(I(n)+1,J(n)+1,K(n)+1);
                        rt_interp(n) = g_base{bid+1}.rt(I(n)+1,J(n)+1,K(n)+1);
                    end

                    Zeta_interp = F_zeta(x_interp,r_interp,rt_interp);

                    for n = 1:length(Zeta_interp)
                        Zeta{bid+1}(I(n)+1,J(n)+1,K(n)+1) = Zeta_interp(n);
                    end
                end

                % Plot surfaces
                H.surface = figure();
                ts_plot_surface(g_base, [0 0 0], H.surface)
                ts_plot_surface(g_base, Zeta, H.surface, 0, 2)
                plot3(x,rt,r,'k.')
                colorbar
                caxis(C);
                
%                 F{f}.Zeta_matrix = Zeta;
                f = f+1;
            end
        end
    end
    
    Zeta = F{end}.Loss;
    if length(F) > 1 && length(F{1}.Loss) == length(F{2}.Loss)
        % Check correlation between runs of multiple bump heights
        dps = []; Zetas = []; bids = []; Is = []; Js = []; Ks = []; 
        for f = 1:length(F)
            Zeta = F{f}.Loss; Zetas = [Zetas Zeta];
            dp = F{f}.dp; dps = [dps ones(size(Zeta)) * dp];
            bid = F{f}.bid; bids = [bids bid];
            I = F{f}.I; Is = [Is I];
            J = F{f}.J; Js = [Js J];
            K = F{f}.K; Ks = [Ks K];
        end

        % Include zero height bump
        dps = [zeros(size(Zeta)) dps];
        Zetas = [ones(size(Zeta)) Zetas];
        bids = [bids(:,1) bids];
        Is = [Is(:,1) Is];
        Js = [Js(:,1) Js];
        Ks = [Ks(:,1) Ks];

        % Calculate correlation
        rho = zeros(size(Zetas,1),1);
        for n = 1:length(rho)
            rho(n) = corr(dps(n,:).',Zetas(n,:).','type','Pearson');
        end

        % Initialise correlation matrices
        Corr = cell(size(g_base));
        for bid = ts_get_bids(g_base)
            Corr{bid+1} = ones(size(g_base{bid+1}.x));
        end

        % Interpolate correlation elsewhere on walls
        F_corr = TriScatteredInterp(X,R,RT,repmat(abs(rho),[9 1]));

        for bid = ts_get_bids(g_base)
            [I,J,K] = ind2sub(size(g_base{bid+1}.mwall), find(g_base{bid+1}.mwall == 0));
            I = I-1; J = J-1; K = K-1;

            x_interp = zeros(length(I),1); r_interp = zeros(length(I),1); 
            rt_interp = zeros(length(I),1);
            for n = 1:length(I)
                x_interp(n) = g_base{bid+1}.x(I(n)+1,J(n)+1,K(n)+1);
                r_interp(n) = g_base{bid+1}.r(I(n)+1,J(n)+1,K(n)+1);
                rt_interp(n) = g_base{bid+1}.rt(I(n)+1,J(n)+1,K(n)+1);
            end

            Corr_interp = F_corr(x_interp,r_interp,rt_interp);

            for n = 1:length(Corr_interp)
                Corr{bid+1}(I(n)+1,J(n)+1,K(n)+1) = Corr_interp(n);
            end
        end

        % Plot correlation map
        H.corr = figure();
        ts_plot_surface(g_base, [0 0 0], H.surface)
        ts_plot_surface(g_base, Corr, H.surface, 0, 2)
        plot3(x,rt,r,'k.')
        colorbar

        % Calculate and plot sensitivities weighted by correlation and filtered by noise
%         Corr_filter = Corr; Noise_filter = F{end}.Zeta_matrix; Zeta = F{end}.Zeta_matrix;
%         Zeta_max = 1; Zeta_min = 1;
%         for bid = ts_get_bids(g_base)
%             Corr_filter{bid+1}(Corr_filter{bid+1} < 0.8) = 0;
%             Corr_filter{bid+1}(Corr_filter{bid+1} > 0.8) = 1;
% 
%             Noise_filter{bid+1}(Noise_filter{bid+1} < mag_noise) = 0;
%             Noise_filter{bid+1}(Noise_filter{bid+1} > mag_noise) = 1;
% 
%             Zeta{bid+1} = Zeta{bid+1} .* Noise_filter{bid+1} .* Corr_filter{bid+1};
%             Zeta{bid+1}(Zeta{bid+1} == 0) = 1;
%             Zeta_max = max(max(max(max(Zeta{bid+1}))),Zeta_max);
%             Zeta_min = min(min(min(min(Zeta{bid+1}))),Zeta_min);
%         end
        
        % Reject values of loss coefficient at bump locations based on
        % noise and linearity
        Corr_filter = abs(rho);
        Corr_filter(Corr_filter < 0.8) = 0;
        Corr_filter(Corr_filter > 0.8) = 1;
        
        Zeta = Zeta .* Corr_filter;
    end
    
    Noise_filter = Zeta;
    Noise_filter(Noise_filter < mag_noise) = 0;
    Noise_filter(Noise_filter > mag_noise) = 1;

    Zeta = Zeta .* Noise_filter;
    Zeta(Zeta == 0) = 1;

    % Record output
    F{end+1} = F{end};
    F{end}.Loss = Zeta;

    F_zeta = TriScatteredInterp(X,R,RT,repmat(Zeta,[9 1]));
    d = max(1 - min(Zeta), max(Zeta) -1);
    C = [1-d 1+d];

    % Initialise loss matrices
    Zeta = cell(size(g_base));
    for bid = ts_get_bids(g_base)
        Zeta{bid+1} = ones(size(g_base{bid+1}.x));
    end

    % Calculate new interpolated loss coefficients on all surfaces
    for bid = ts_get_bids(g_base)
        [I,J,K] = ind2sub(size(g_base{bid+1}.mwall), find(g_base{bid+1}.mwall == 0));
        I = I-1; J = J-1; K = K-1;

        x_interp = zeros(length(I),1); r_interp = zeros(length(I),1); 
        rt_interp = zeros(length(I),1);
        for n = 1:length(I)
            x_interp(n) = g_base{bid+1}.x(I(n)+1,J(n)+1,K(n)+1);
            r_interp(n) = g_base{bid+1}.r(I(n)+1,J(n)+1,K(n)+1);
            rt_interp(n) = g_base{bid+1}.rt(I(n)+1,J(n)+1,K(n)+1);
        end

        Zeta_interp = F_zeta(x_interp,r_interp,rt_interp);

        for n = 1:length(Zeta_interp)
            Zeta{bid+1}(I(n)+1,J(n)+1,K(n)+1) = Zeta_interp(n);
        end
    end

    % Plot Loss coefficients
    H.corrected = figure();
    ts_plot_surface(g_base, [0 0 0], H.surface)
    ts_plot_surface(g_base, Zeta, H.surface, 0, 2)
    plot3(x,rt,r,'k.')
    colorbar
    caxis(C)
end


end
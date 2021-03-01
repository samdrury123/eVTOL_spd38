function [H] = ts_plot_structured(g,property,IJK,type,range,bids,H)
% Plot contours of a variable throughout the domain
% Select res:       'Coarse'    'Fine'
% Select col:       'Bone'      'Jet'
% Select coords:    'Cart'      'Pol'
% Select view:      'Trav'      'B2B'       'Mer'
% Select mesh:      'on'        'off
% Select rep:       '%numeric'

if exist('type','var') == 0
    type.res = 'Fine'; type.col = 'Jet'; type.coords = 'Pol';...
        type.view = 'B2B'; type.rep = 1; type.mesh = 'off';
else
    if isfield(type,'res') == 0
        type.res = 'Fine';
    end
    if isfield(type,'col') == 0
        type.col = 'Jet';
    end
    if isfield(type,'coords') == 0
        type.coords = 'Pol';
    end
    if isfield(type,'view') == 0
        type.view = 'B2B';
    end    
    if isfield(type,'mesh') == 0
        type.mesh = 'off';
    end    
    if isfield(type,'rep') == 0
        type.rep = 1;
    end        
end

g = ts_secondary(g);
if exist('H','var') == 0
    H = figure();
else
    figure(H.fig)
    subplot(H.sub)
end
hold on

% Select blocks to loop over
if exist('bids','var') == 0
    bids = [];
end
if isempty(bids) == 1
    bids = ts_get_bids(g);
end

% Gather coordinate and property data
X = cell(size(g)); R = cell(size(g)); RT = cell(size(g));
Y = cell(size(g)); Z = cell(size(g)); P = cell(size(g)); nblade = zeros(size(g));
for bid = bids
    write_vars = 1;
    if IJK(1) < g{bid+1}.attribute.ni && IJK(1) > -1
        I = IJK(1);
    elseif IJK(1) >= g{bid+1}.attribute.ni
        write_vars = 0;
    else
        I = 0:g{bid+1}.attribute.ni-1;
    end
    if IJK(2) < g{bid+1}.attribute.nj && IJK(2) > -1
        J = IJK(2);
    elseif IJK(2) >= g{bid+1}.attribute.nj
        write_vars = 0;
    else
        J = 0:g{bid+1}.attribute.nj-1;
    end
    if IJK(3) < g{bid+1}.attribute.nk && IJK(3) > -1
        K = IJK(3);
    elseif IJK(3) >= g{bid+1}.attribute.nk
        write_vars = 0;
    else
        K = 0:g{bid+1}.attribute.nk-1;
    end    
    
    if write_vars == 1
        X{bid+1} = double(squeeze(g{bid+1}.x(I+1,J+1,K+1)));
        R{bid+1} = double(squeeze(g{bid+1}.r(I+1,J+1,K+1)));
        RT{bid+1} = double(squeeze(g{bid+1}.rt(I+1,J+1,K+1)));
        Y{bid+1} = double(squeeze(g{bid+1}.y(I+1,J+1,K+1)));
        Z{bid+1} = double(squeeze(g{bid+1}.z(I+1,J+1,K+1)));
        P{bid+1} = double(squeeze(g{bid+1}.(property)(I+1,J+1,K+1)));
        nblade(bid+1) = double(g{bid+1}.bv.nblade);
    end
end

% Remove empty blocks
X = X(~cellfun(@isempty, X));
R = R(~cellfun(@isempty, R));
RT = RT(~cellfun(@isempty, RT));
Y = Y(~cellfun(@isempty, Y));
Z = Z(~cellfun(@isempty, Z));
P = P(~cellfun(@isempty, P));
nblade = nblade(nblade ~= 0);

% Define range of plot
if exist('range','var') == 0
    range = [];
end
if isempty(range) == 1
    V_min = inf; V_max = -inf;
    for n = 1:length(P)
        V_min = min(min(min(P{n})),V_min);
        V_max = max(max(max(P{n})),V_max);
        if strcmp(type.res,'Coarse') == 1
            V = linspace(V_min,V_max,20);
        elseif strcmp(type.res,'Fine') == 1
            V = linspace(V_min,V_max,100);
        else
            disp('Specify resolution')        
        end
    end
elseif length(range) == 2
    V_min = range(1); V_max = range(2);
    if strcmp(type.res,'Coarse') == 1
        V = linspace(V_min,V_max,20);
    elseif strcmp(type.res,'Fine') == 1
        V = linspace(V_min,V_max,100);
    else
        disp('Specify resolution')
    end
else
    V = range;
end

% Repeat matrices for multiple passages if necessary
if type.rep > 1
    if length(P) == 1
        T = repmat(RT{n}(:,1:end-1)./R{n}(:,1:end-1),[1 type.rep]);
        dT = [];
        for m = 1:type.rep
            dT = [dT -ones(size(R{n}(:,1:end-1))) * (m-1) * 2 * pi / nblade(n)];
        end

        T = T + dT;
        
        if strcmp(type.view,'Trav') == 1
            T_mean = mean(mean(T));
            T = T - T_mean;
        end

        R{n} = repmat(R{n}(:,1:end-1),[1 type.rep]);
        X{n} = repmat(X{n}(:,1:end-1),[1 type.rep]);
        P{n} = repmat(P{n}(:,1:end-1),[1 type.rep]);

        RT{n} = T.*R{n};
        [Y{n} Z{n}] = pol2cart(T, R{n});
    else
        for n = 1:length(P)
            for m = 2:type.rep
                T = RT{n}./R{n};
                dT = -ones(size(R{n})) * (m-1) * 2 * pi / nblade(n);
                T = T + dT;
                
                R = [R ; R{n}];
                RT = [RT ; T.*R{n}];
                X = [X ; X{n}];
                P = [P ; P{n}];
                [y z] = pol2cart(T, R{n});
                Y = [Y ; y];
                Z = [Z ; z];
            end
        end                
    end
end        

% Select coordinates to plot over
if strcmp(type.view,'Trav') == 1
    if strcmp(type.coords,'Pol') == 1
        A = RT; B = R;
    elseif strcmp(type.coords,'Cart') == 1
        A = Z; B = Y;
    end
elseif strcmp(type.view,'B2B') == 1
    if strcmp(type.coords,'Pol') == 1
        A = X; B = RT;
    elseif strcmp(type.coords,'Cart') == 1
        A = X; B = Z;
    else
        disp('Specify coordinates')           
    end
elseif strcmp(type.view,'Mer') == 1
    if strcmp(type.coords,'Pol') == 1
        A = X; B = R;
    elseif strcmp(type.coords,'Cart') == 1
        A = X; B = Y;
    else
        disp('Specify coordinates')    
    end        
else
    disp('Specify view')
end

% Contour plot the data
for n = 1:length(P)
    
    contourf(A{n},B{n},P{n},V);
    
    % Plot mesh lines on top if required
    if strcmp(type.mesh,'on') == 1
        plot(A{n},B{n},'k-')
        plot(A{n}.',B{n}.','k-')
    elseif strcmp(type.mesh,'off') ~= 1
        disp('Specify mesh')
    end
end

% Select the colouring
if strcmp(type.col,'Bone') == 1
    colormap(flipud(bone)); 
elseif strcmp(type.col,'Jet') == 1
    colormap(jet)
else
    disp('Specify colourmap')
end

% Remove black lines from plot if fine contours used
if strcmp(type.res,'Fine') == 1 || strcmp(type.mesh,'on') == 1 ||...
    length(P) > 1
    shading('flat')
end

axis auto
axis equal
caxis([min(V) max(V)])
colorbar

end
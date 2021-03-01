function c = ts_unstructured_cut(g,x_cut,bids,plotnames)
% TS_UNSTRUCTURED_CUT  Triangulated cut plane at specified axial coordinate

% Use specified blocks or only consider blocks that contain the target x coordinate
if exist('bids','var') == 0 || isempty(bids) == 1
    bids = [];
    for bid = ts_get_bids(g)
        if max(g{bid+1}.x(:)) > x_cut && min(g{bid+1}.x(:)) < x_cut
            bids = [bids bid];
        end
    end
end

% Preallocate triangle arrays
fv = cell(size(g));

% Primary flow variables for interpolation
varnames = {'ro' 'rovx' 'rovr' 'rorvt' 'roe' 'trans_dyn_vis' 'phi'};

% Cut the primary flow variables at constant x
for bid = bids
    for v = 1:length(varnames)
        
        % Get coordinates and current variable
        x = g{bid+1}.x; r = g{bid+1}.r; rt = g{bid+1}.rt; 
        var = g{bid+1}.(varnames{v});
        
        % Take isosurface
        [faces,vertices] = isosurface(r,rt,var,x,x_cut,var);
        
        % Record variables and coordinates
        fv{bid+1}.faces = faces; 
        fv{bid+1}.x = ones(size(vertices,1),1) * x_cut;
        fv{bid+1}.r = vertices(:,1);
        fv{bid+1}.rt = vertices(:,2);
        fv{bid+1}.(varnames{v}) = vertices(:,3);
    end
end

% Get application and block variables
c.av = g{bids(1)+1}.av; c.bv = g{bids(1)+1}.bv;

% Include coordinates in variable names
varnames = [varnames {'x' 'r' 'rt'}];

% Preallocate variables
c.faces = [];
for v = 1:length(varnames)
    c.(varnames{v}) = [];
end

% Collate all blocks together into one triangulated cut
for bid = bids
    
    % Collate face indices
    ni = size(c.ro,1);
    c.faces = [c.faces ; fv{bid+1}.faces + ni];
    
    % Loop over primary flow variables and collate
    for v = 1:length(varnames)
        c.(varnames{v}) = [c.(varnames{v}) ; fv{bid+1}.(varnames{v})];
    end

end

% Plot specified variable on cut plane
if exist('plotnames','var') ~= 0 && isempty(plotnames) == 0
    
    % Arrange in cell array if just one variable chosen
    if ischar(plotnames)
        plotnames = {plotnames};
    end
    
    % Loop over all plot variables
    for a = 1:length(plotnames)
        figure; hold on;

        % Plot face patches of triangulated cut
        patch('Vertices',[c.z c.y c.(plotnames{a})],'Faces',c.faces,...
            'FaceVertexCData',c.(plotnames{a}),'FaceColor','interp','edgecolor','interp');
        
        % Set colours and axes
        axis equal; axis tight
        colormap(jet(64*4)); view(2);
    end
end

end

function [cut,average] = cut_averages(cut)

% Calculate the area of the face triangles using Heron's formula
for bid = 0:length(cut)-1
    if ~isempty(cut{bid+1})
        y = cut{bid+1}.y;
        z = cut{bid+1}.z;
        face = cut{bid+1}.faces;
        a = ((y(face(:,2)) - y(face(:,1))).^2 + (z(face(:,2)) - z(face(:,1))).^2).^0.5;
        b = ((y(face(:,3)) - y(face(:,2))).^2 + (z(face(:,3)) - z(face(:,2))).^2).^0.5;
        c = ((y(face(:,1)) - y(face(:,3))).^2 + (z(face(:,1)) - z(face(:,3))).^2).^0.5;
        d = sort([a b c],2,'descend');
        a = d(:,1);
        b = d(:,2);
        c = d(:,3);
        cut{bid+1}.area = 0.25*sqrt(4*a.*a.*b.*b - (a.*a + b.*b - c.*c).^2);
        
        if any(imag(cut{bid+1}.area)~=0)
            cut{bid+1}.area = real(cut{bid+1}.area);
        end
    end
end

% Individual block averages
for bid = 0:length(cut)-1
    
    if ~isempty(cut{bid+1})
        vars = fieldnames(cut{bid+1});

        face = cut{bid+1}.faces;
        f1 = face(:,1);
        f2 = face(:,2);
        f3 = face(:,3);
        
        % Face mass flow rates
        
        rovx = cut{bid+1}.rovx;
        rovx = (rovx(f1) + rovx(f2) + rovx(f2))/3;
        mdot = cut{bid+1}.area.*rovx;
        mdotTot = sum(mdot);
%         mdotTot = sum(mdot(mdot>0));
        cut{bid+1}.average.mdot = mdotTot;

        for j = 1:length(vars)
            var = cut{bid+1}.(vars{j});
            
            % Calculate face average values
            if size(var,1) > 1
                if length(var) >= max([f1 ; f2 ; f3])
                    varAvg = (var(f1) + var(f2) + var(f3))/3;
                    cut{bid+1}.average.(vars{j}) = (sum(varAvg.*mdot))/mdotTot;
                else
                    cut{bid+1}.average.(vars{j}) = nan;
                end
            end
        end
    end
end

% Overall mass-averages
for j = 1:length(vars)
    var = zeros(size(cut));
    mdot = zeros(size(cut));
    for bid = 0:length(cut)-1
        if ~isempty(cut{bid+1})  
            if size(cut{bid+1}.(vars{j}),1) > 1
                var(bid+1) = cut{bid+1}.average.(vars{j});
                mdot(bid+1) = cut{bid+1}.average.mdot;
            end
        end
    end
    average.(vars{j}) = sum(mdot.*var)/sum(mdot);
    average.mdot = sum(mdot);
end

end
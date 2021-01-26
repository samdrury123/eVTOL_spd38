function [d,c] = exp_fhp_squares(c,P_FHP,probe,Iota,Tau,coeffs,plot_stuff)
% Use a least squares method to interpret five hole probe data

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Default to Hodson & Dominy style coefficients
if exist('coeffs','var') == 0
    coeffs = 'hodson';
end

% Process five hole probe data with traditional coefficients
d = exp_fhp_calculate(c,P_FHP,probe,Iota,Tau,coeffs);

% Extract pressure from calibration map
c.P_side = squeeze(c.P_raw(:,:,end,4:7));
c.P_av = mean(c.P_side,3);
c.P_cen = c.P_raw(:,:,end,3);
c.Po = c.P_raw(:,:,end,1); 
c.P = c.P_raw(:,:,end,2);

% Calculate individual hole coefficients from calbration
if strcmp(coeffs,'hodson') == 1
    c.C = (c.P_side - repmat(c.P_cen,[1 1 4])) ./ repmat(c.P_cen - c.P_av,[1 1 4]);
elseif strcmp(coeffs,'curtis') == 1
    P_sort = sort(c.P_side,3);
    P_ref = (2/3) * c.P_cen + (1/3) * mean(P_sort(:,:,[3 4]),3) - mean(P_sort(:,:,[1 2]),3);
    c.C = (c.P_side - repmat(c.P_cen,[1 1 4])) ./ repmat(P_ref,[1 1 4]);
end

% Plot individual coefficients
if plot_stuff == 1
    figure('Position',[140 300 1700 400]);
    for n = 1:4
        subplot(1,4,n); mesh(c.Iota,c.Tau,c.C(:,:,n)); 
        xlabel('Iota'); ylabel('Tau'); zlabel(['Coeff #' num2str(n)]);
    end
end

% Reshape pressure data
if size(P_FHP,3) == 5
    [nj,nk,~] = size(P_FHP);
else
    nj = size(P_FHP,1);
    nk = 1;
end
P_FHP = reshape(P_FHP,nj*nk,5);

% Calculate coefficients from traverse
e.P_side = P_FHP(:,2:5);
e.P_cen = P_FHP(:,1);
if strcmp(coeffs,'hodson') == 1
    e.P_av = mean(e.P_side,2);
    e.C = (e.P_side - repmat(e.P_cen,[1 4])) ./ repmat(e.P_cen - e.P_av,[1 4]);
elseif strcmp(coeffs,'curtis') == 1
    P_sort = sort(e.P_side,2);
    P_ref = (2/3) * e.P_cen + (1/3) * mean(P_sort(:,[3 4]),2) - mean(P_sort(:,[1 2]),2);
    e.C = (e.P_side - repmat(e.P_cen,[1 4])) ./ repmat(P_ref,[1 4]);
end

% Open figure window
if plot_stuff == 1
    h.coeffs = figure(); hold on; grid on; box on; xlabel('Yaw Angle'); ylabel('Pitch Angle');
    Col = jet(size(e.C,1));
end

% Loop over all data points and calculate error
xy_int = nan(size(e.C,1),4,2);
for p = 1:size(e.C,1)
    
    % Calculate contour levels of each probe hole coefficient
    Ci = cell(4,1);
    for n = 1:4
        Ci{n} = contourcs(contourc(c.Iota(:,1),c.Tau(1,:),c.C(:,:,n)',e.C(p,n)*ones(1,2)));
    end
    
    % Calculate length of all contour line sections and choose longest continuous one
    xy_ci = cell(4,1);
    for n = 1:4
        if length(Ci{n}) > 1
            S = [];
            for m = 1:length(Ci{n})
                S = [S ; sum(diff(Ci{n}(m).X).^2 + diff(Ci{n}(m).Y).^2)];
            end
            [~,m] = sort(S);
            xy_ci{n} = [Ci{n}(m(end)).X' Ci{n}(m(end)).Y'];
        else
            xy_ci{n} = [Ci{n}.X' Ci{n}.Y'];
        end
    end
    c.xy_ci{p} = xy_ci;
    
    % Check contour lines have been found
    empty = 0;
    for n = 1:4
        if isempty(xy_ci{n}) == 1; empty = 1; end;
    end
    if empty == 1; continue; end;
    
    % Truncate contour lines to vicinity of point on calibration map
    dAlpha = 5; dBeta = 5;
    for n = 1:4
        q = xy_ci{n}(:,1) < d.Alpha(p)-Iota + dAlpha & xy_ci{n}(:,1) > d.Alpha(p)-Iota - dAlpha & ...
            xy_ci{n}(:,2) < d.Beta(p) + dBeta & xy_ci{n}(:,2) > d.Beta(p) - dBeta;
        xy_ci{n}(q == 0,:) = [];
    end
    
    % Calculate intersections between side hole coefficients
    xy_ci{5} = xy_ci{1};
    for n = 1:4
        if size(xy_ci{n},1) > 2 && size(xy_ci{n+1},1) > 2
            [x_temp,y_temp] = intersections(xy_ci{n}(:,1),xy_ci{n}(:,2),...
                xy_ci{n+1}(:,1),xy_ci{n+1}(:,2));
            if numel(x_temp) == 1 
                xy_int(p,n,1) = x_temp;
                xy_int(p,n,2) = y_temp;
            end
        else
            continue
        end
    end
    
    % Plot intersections of contour lines around calculated point on map
    if plot_stuff == 1
        plot(xy_int(p,[1:end 1],1),xy_int(p,[1:end 1],2),'-','Color',Col(p,:))
        plot(d.Alpha(p)-Iota,d.Beta(p),'.','Color',Col(p,:))
    end
end

% Calculate areas of error quadrilaterals
% dxy_1 = reshape(xy_int(:,1,:) - xy_int(:,3,:),[],2);
% dxy_2 = reshape(xy_int(:,2,:) - xy_int(:,4,:),[],2);
% A = 0.5 * abs(dxy_1(:,1) .* dxy_2(:,2) - dxy_2(:,1) .* dxy_1(:,2));
        
% Calculate average diagonal length of error quadrilaterals / polygons
dxy = xy_int - cat(3,repmat(reshape(d.Alpha-Iota,[],1),[1 4]),repmat(reshape(d.Beta,[],1),[1 4]));
A = zeros(size(e.C,1),1);
for p = 1:size(e.C,1)
    A(p) = mean(sum(dxy(p,isnan(dxy(p,:,1)) == 0,:).^2,3).^0.5,2);
end

% Reshape error array into original size
d.A = reshape(A,nj,nk);

% Reshape other arrays
varnames = fieldnames(d);
for v = 1:length(varnames)
    d.(varnames{v}) = reshape(d.(varnames{v}),[nj nk]);
end


end
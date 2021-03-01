function h = exp_fhp_plot(c,h,cols)
% EXP_FHP_PLOT  Plot a calibration map

% Generate new figure if not specified
if exist('h','var') == 0 || isempty(h) == 1
    h = figure('Position',[1 29 1920 985]);
end

% Generate colours if not specified
if exist('cols','var') == 0 || isempty(cols) == 1
    cols = lines(size(c.C_Alpha,3));
end

% Plot spiders webs of pitch and yaw coefficients
subplot(1,3,1); hold on; grid on; box on; axis equal;
for n = 1:size(c.C_Alpha,3)
    plot(c.C_Alpha(:,:,n),c.C_Beta(:,:,n),'-','Color',cols(n,:))
    plot(c.C_Alpha(:,:,n).',c.C_Beta(:,:,n).','-','Color',cols(n,:))
end
xlabel('Yaw Coefficient'); ylabel('Pitch Coefficient');

% Calculate limits
n_contours = 22;
v = [min(min(c.Iota)) max(max(c.Iota)) min(min(c.Tau)) max(max(c.Tau))];

% Plot individual coefficient contours
if size(c.C_Alpha,2) > 1 && size(c.C_Alpha,1) > 1
    for n = 1:size(c.C_Alpha,3)

        % Yaw coefficient
        subplot(2,3,2); hold on; grid on; box on;
        xlabel('Yaw Angle'); ylabel('Pitch Angle'); title('Yaw Coefficient');
        contour(c.Iota,c.Tau,c.C_Alpha(:,:,n),get_range(c.C_Alpha,n_contours),'LineColor',cols(n,:));
        axis(v);

        % Pitch coefficient
        subplot(2,3,3); hold on; grid on; box on;
        xlabel('Yaw Angle'); ylabel('Pitch Angle'); title('Pitch Coefficient');
        contour(c.Iota,c.Tau,c.C_Beta(:,:,n),get_range(c.C_Beta,n_contours),'LineColor',cols(n,:));
        axis(v);

        % Yaw coefficient
        subplot(2,3,5); hold on; grid on; box on;
        xlabel('Yaw Angle'); ylabel('Pitch Angle'); title('Stagnation Pressure Coefficient');
        contour(c.Iota,c.Tau,c.C_Po(:,:,n),get_range(c.C_Po,n_contours),'LineColor',cols(n,:));
        axis(v);

        % Yaw coefficient
        subplot(2,3,6); hold on; grid on; box on;
        xlabel('Yaw Angle'); ylabel('Pitch Angle'); title('Static Pressure Coefficient');
        contour(c.Iota,c.Tau,c.C_P(:,:,n),get_range(c.C_P,n_contours),'LineColor',cols(n,:));
        axis(v);
    end
end

end

function v = get_range(x,n_contours)
% Return a linear range of 1.5 standard deviations from the mean
x = x(:);
v = linspace(mean(x)-1.5*std(x),mean(x)+1.5*std(x),n_contours);
end
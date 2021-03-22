function a_new = add_box(h,add_grid)
% ADD_BOX  Overlay a box and grid over current plot

% Default to current figure window
if exist('h','var') == 0
    h = gcf;
end

% Default to add grid lines too
if exist('add_grid','var') == 0
    add_grid = 1;
end

% Switch to figure window and get axes handle
figure(h); a = gca;

% Overlay new axes
a_new = axes('Position',a.Position,'Color','none','XTick',a.XTick,'YTick',a.YTick,...
    'XTickLabel',[],'YTickLabel',[],'XLim',a.XLim,'YLim',a.YLim,'Box','on');

% Add grid lines
if add_grid == 1
    grid on;
end


end
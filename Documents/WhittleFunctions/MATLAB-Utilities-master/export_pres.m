function export_pres(h,filename,xy)
% EXPORT_PRES  Export a figure as eps format for publication in a paper or presentation
%
%   EXPORT_PRES(h,filename,xy)
%
%   h - figure handle to save to file
%   filename - string with name to save under
%   xy - vector of figure size in mm or a string to specify a default
%
%   xy can be a string:
%       'asme' to publish a 90 x 67.5 mm figure
%       'full' to publish a full width figure for a slide 200 x 125 mm
%       'half' to publsih a half width figure for a slide 100 x 125 mm

% Default to full presentation figure
if exist('xy','var') == 0
    xy = 'full';
end

% Default sizes, all x2 real size on the screen
if strcmp(xy,'asme') == 1
    set(h,'Position',[10 60 612 460])
elseif strcmp(xy,'full') == 1
    set(h,'Position',[10 60 883 527])
elseif strcmp(xy,'half') == 1
    set(h,'Position',[10 60 442 527])
elseif strcmp(xy,'poster') == 1
    set(h,'Position',[10 60 952 529]);
else
    set(h,'Position',[10 60 2*xy(1)/0.294 2*xy(2)/0.294])
end

% Get handles for all axes
h_axes = findall(h,'type','axes');

% Set axes colours to black
set(h_axes,'xcolor',[0 0 0],'ycolor',[0 0 0],'zcolor',[0 0 0])

% Minimum line thickness and marker size
h_line = findobj(h,'type','line');
for n = 1:length(h_line)
    set(h_line(n),'LineWidth',max(h_line(n).LineWidth,1));
%     if get(h_line(n),'LineWidth') == 0.5; set(h_line(n),'LineWdith',1); end;
    if strcmp(h_line(n).Marker,'.') == 1
        set(h_line(n),'MarkerSize',max(h_line(n).MarkerSize,8.5));
    end
%     if strcmp(h_line(n).Marker,'.') == 1 && h_line(n).MarkerSize == 10
%         set(h_line(n),'MarkerSize',max(h_line(n).MarkerSize,20));
%     end
end

% Axis label and tick label font and size
for n = 1:length(h_axes)
    h_label = get(h_axes,{'XLabel' 'YLabel' 'ZLabel'});
    for m = 1:numel(h_label)
        set(h_label{m},'FontSize',14,'FontName','Arial');
    end
    set(h_axes,'FontSize',12,'FontName','Arial');
end

% Line width of borders and grid lines
for n = 1:length(h_axes)
    set(h_axes(n),'LineWidth',1);
end

% Legend font size
h_leg = findobj(h,'Tag','legend');
set(h_leg,'FontSize',12,'FontName','Arial');

% White background
set(h,'Color',[1 1 1]);

% Export as eps format
if any(strcmp(filename(end-2:end),{'pdf' 'eps'})) == 1
%     export_fig(h,filename,'-painters');
    export_fig(filename,h)
else
%     export_fig(h,filename,'-m4');
    export_fig(filename,'-m4',h);
end

end


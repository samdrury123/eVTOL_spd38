function [h,c] = bl_spline_plot(b,h,col)
% BL_PLOT_SPLINE  Visualise a spline blade definition
%
%   [h,c] = BL_PLOT_SPLINE(b,h,C)
%   
%   b - spline data structure
%   h - optional figure handle passed to function
%   col - rgb colour vector
%   c - output struct of evaluated spline definition
%   

if exist('col','var') == 0 || isempty(col) == 1
    col = [0 0 0];
end
if exist('h','var') == 0 || isempty(h) == 1
    h = figure();
end

figure(h);
set(h,'Position',[1 82 1280 895])

% Evaluate all splines
if isnumeric(b.chi_le) == 0
    r_nondim = linspace(0,1,200);
    c = bl_spline_eval(b,r_nondim);
else
    c = b; r_nondim = b.r_nondim;
end

% Plot blade spline parameters
figure(h); n = 2; m = 7;

% Create new subplots or get handles of existing ones
if isempty(h.Children) == 1
    for o = 1:n*m; h_sub(o) = subplot(n,m,o); hold on; grid on; box on; end;
else
    for o = 1:n*m; h_sub(o) = h.Children(o); end;
end

plot(h_sub(1),c.chi_le,r_nondim,'-','Color',col)
xlabel(h_sub(1),'Leading Edge Angle')

plot(h_sub(2),c.dcam_le,r_nondim,'-','Color',col)
xlabel(h_sub(2),'Front Camber Gradient')

plot(h_sub(3),c.dcam_te,r_nondim,'-','Color',col)
xlabel(h_sub(3),'Rear Camber Gradient')

plot(h_sub(4),c.qcam,r_nondim,'-','Color',col)
xlabel(h_sub(4),'Camber Symmetry')

plot(h_sub(5),c.chi_te,r_nondim,'-','Color',col)
xlabel(h_sub(5),'Trailing Edge Angle')

plot(h_sub(6),c.tchord,r_nondim,'-','Color',col)
xlabel(h_sub(6),'True Chord')

plot(h_sub(7),c.rad_le,r_nondim,'-','Color',col)
xlabel(h_sub(7),'Leading Edge Radius')

plot(h_sub(8),c.thick_max,r_nondim,'-','Color',col)
xlabel(h_sub(8),'Max Thickness')

plot(h_sub(9),c.s_thick_max,r_nondim,'-','Color',col)
xlabel(h_sub(9),'Position Max Thickness')

plot(h_sub(10),c.rad_thick_max,r_nondim,'-','Color',col)
xlabel(h_sub(10),'Max Thickness Radius')

plot(h_sub(11),c.thick_te,r_nondim,'-','Color',col)
xlabel(h_sub(11),'Trailing Edge Thickness')

plot(h_sub(12),c.wedge_te,r_nondim,'-','Color',col)
xlabel(h_sub(12),'Trailing Edge Wedge')

plot(h_sub(13),c.sweep,r_nondim,'-','Color',col)
xlabel(h_sub(13),'True Sweep')

plot(h_sub(14),c.lean,r_nondim,'-','Color',col)
xlabel(h_sub(14),'True Lean')

% Set axes limits
a = findall(h,'type','axes');
for n = 1:length(a)
    axes(a(n)); axis tight; v = axis; axis([v(1:2) 0 1]);
end


end
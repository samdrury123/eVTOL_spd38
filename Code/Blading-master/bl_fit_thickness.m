function t = bl_fit_thickness(s_cl,S,thick_te,plot_stuff)
% BL_FIT_THICKNESS  Fit a thickness distribution smoothly in shape space
%
%   t = BL_FIT_THICKNESS(s,S,thick_te,plot_stuff)
%
%   s - non-dimensional coordinate along camberline
%   S - shape space parameter
%   thick_te - trailing edge thickness
%   plot_stuff - 0 or 1 for showing working

if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Calculate thickness
S(1) = 0; S(end) = -inf;
thick = S .* (s_cl.^0.5 .* (1 - s_cl)) + s_cl.* thick_te;
% z = S .* (s.^(1/3) .* (1 - s)) + s.* (thick_te / true_chord);
thick(1) = 0; thick(end) = 0;

% Cut out trailing edge
if thick_te > 0
    [~, curv] = grad_mg(s_cl,thick);
    dcurv = diff(curv.');
    i_te = find(abs(dcurv(2:end) - dcurv(1:end-1)) > 100 & s_cl(2:end-1) > 0.9,1,'first');
else
    i_te = length(s_cl);
end
s_max = s_cl(i_te)-0.005;

% Cut leading edge if inflexional in shape space
q = s_cl < 0.1 & s_cl > 0.005; 
s_temp = linspace(0.005,0.05,20); 
S_temp = interp1(s_cl,S,s_temp,'pchip');
[~,i] = findpeaks(S_temp,'MINPEAKDISTANCE',3,'SORTSTR','ascend');
if isempty(i) == 0
    s_min = s_temp(i(end));
else
    s_min = 0.02;
end

% Plot original thickness distribution
if plot_stuff == 1
    figure();
    subplot(2,1,1); hold on; grid on; box on;
    plot(s_cl,S,'r-');
    plot([s_min s_max],interp1(s_cl,S,[s_min s_max],'pchip'),'r.');
    q = s_cl > 0.02 & s_cl < 0.98;
    axis([0 1 min(S(q))-0.1 max(S(q))+0.1]);
    subplot(2,1,2); hold on; grid on; box on;
    plot(s_cl,thick,'r-');
    plot([s_min s_max],interp1(s_cl,thick,[s_min s_max],'pchip'),'r.');
end

%% Old fit types

% Fit quartic
% p = polyfit(s_fit,S_fit,4);

% 2 part cubic spline fit, split at 10% chord
% spline_fit = fnxtr(spap2([0 0 0 0 0.10 1 1 1 1],4,s,S),4);

% 2 part quartic spline fit, split at 15% chord
% spline_fit = spap2([0 0 0 0 0 0.15 1 1 1 1 1],5,s_cl,S);

% Uniform distribution in shape space before fit and cut out LE & TE
% s_fit = linspace(0.02,s_max-0.005,100); 
% S_fit = interp1(s_cl,S,s_fit,'pchip');

% Weighted distribution in shape space
% s_fit = [hyperbolic_bunch(50,0.02,0.05) * 0.48 + 0.02 linspace(0.51,1,50)];
% S_fit = interp1(s,S,s_fit,'pchip');

% 2 part cubic spline fit, split at 15% chord
% s_fit = linspace(s_min,s_max,100)'; S_fit = smooth(interp1(s_cl,S,s_fit,'pchip'));
% spline_fit = fnxtr(spap2([0.05 0.05 0.05 0.05 0.15 1 1 1 1],4,s_fit,S_fit),4);
% S = fnval(spline_fit,s_cl);

% Custom power fit up to 20% chord, polynomial fit on aft section
% q_0 = [0.07 -0.7 2 0.05 0];
% q_0 = [0.07 -0.7 0.05 0];
% s_fit = linspace(s_min,s_max,100)'; S_fit = interp1(s_cl,S,s_fit,'pchip');
% 
% if plot_stuff == 1
%     options = optimset('Display','iter'); 
% else
%     options = optimset('Display','off');
% end
% 
% [q,fval] = lsqcurvefit(@fun_thick,q_0,s_fit,S_fit);
% 
% S = fun_thick(q,s_cl);

% Automated curve fitting function for arbitary curves 
% s_fit = linspace(s_min,s_max,100)'; S_fit = interp1(s_cl,S,s_fit,'pchip');
% % ft = fittype('a*(x+0.03)^b+c*(x+0.03)^d+e*x^3','independent','x','dependent','y');
% % opts = fitoptions(ft); opts.StartPoint = [0.07 -0.7 2 0.05 0.6];
% ft = fittype('a*log(x) + b*(x+e)^c + d*x','independent','x','dependent','y');
% opts = fitoptions(ft); opts.Lower = [-Inf -Inf -Inf -Inf 0];
% opts.StartPoint = [-0.5 2.5 0.3 -0.9 0.02];
% f_1 = fit(s_fit,S_fit,ft,opts);
% S = f_1(s_cl);

%% 2 part cubic spline fit but stretch out leading edge in shape space

% Determine weighting for spline fit
s_fit = unique([linspace(s_min,s_max,101) linspace(0.1,0.45,47) linspace(0.05,0.1,11) ])';
S_fit = interp1(s_cl,S,s_fit,'pchip');

% Define cubic stretching function
s_split = 0.11; s_stretch = 0.08;
A = [0 0 0 1 ; s_split^3 s_split^2 s_split 1 ; 3*s_split^2 2*s_split 1 0 ; 6*s_split 2 0 0];
b = [s_stretch ; 0 ; 0 ; 0];
x = A\b;
s_fit(s_fit < s_split) = s_fit(s_fit < s_split) - polyval(x,s_fit(s_fit < s_split));

% Fit stretched shape space distribution with a two part spline 
s_j = 0.3;
spline_fit = fnxtr(spap2([-s_stretch *ones(1,4) s_j 1 1 1 1],4,s_fit,S_fit),4);

% Ensure spline goes through position of max thickness
[~,i] = max(thick); s_thick_max = s_cl(i);
spline_fit.coefs = spline_fit.coefs * interp1(s_cl,S,s_thick_max,'pchip') / fnval(spline_fit,s_thick_max);

% Evaluate shape space at stretched points
s_interp = s_cl; 
s_interp(s_interp < s_split) = s_interp(s_interp < s_split) - polyval(x,s_interp(s_interp < s_split));
S = fnval(spline_fit,s_interp);

% Plot stretched fit
if plot_stuff == 1
    subplot(2,1,1); 
    plot(s_fit,S_fit,'g-')
    plot(s_interp,S,'c-')
end

%% Calculate new thickness

% Transform thickness back from shape space
% t.S = polyval(p,s);
thick = S .* (s_cl.^0.5 .* (1 - s_cl)) + s_cl * thick_te;
% z = t.S .* (s.^(1/3) .* (1 - s)) + s.* (thick_te / true_chord);

% Calculate new thickness parameters
[~, i] = max(thick); s_temp = s_cl(i);
q = s_cl > s_temp - 0.05 & s_cl < s_temp + 0.05;
s_temp = linspace(s_temp - 0.05,s_temp + 0.05,1000).';
th = polyval(polyfit(s_cl(q),thick(q),4),s_temp);
[~, i] = max(th); s_thick_max = s_temp(i);

% Find trailing edge thickness and wedge angle
q = s_cl > 0.96 & s_cl < 0.98;
p = polyfit(s_cl(q),thick(q),2);
thick_te = polyval(p,1);
wedge_te = atand(-polyval(polyder(p),1));

% Leading edge radius of curvature
p = CircleFitByTaubin([flipud([s_cl(1:5) -thick(1:5)]) ; s_cl(1:5) thick(1:5)]);
rad_le = p(3);

% Radius of curvature at max thickness
[~,d2tds2] = grad_mg(s_cl,thick);
rad_thick_max = -1 / interp1(s_cl,d2tds2,s_thick_max,'pchip');

% Plot fitted thickness distribution
if plot_stuff == 1
    subplot(2,1,1);
    plot(s_cl,S,'b-');
    subplot(2,1,2);
    plot(s_cl,thick,'b-');
    legend('Original','Fit');
    axis([0 1 0 1.1]);
end

% Record parameters
t.s_cl = s_cl; t.thick = thick; t.S = S;
t.s_thick_max = s_thick_max; t.thick_te = thick_te; t.wedge_te = wedge_te;
t.rad_le = rad_le; t.rad_thick_max = rad_thick_max;

end

function S = fun_thick(q,s)
% Calculate square errors to curve fit shape space thickness distribution

% Define join
s_split = 0.03;

% Calculate first region of shape distribution by power fit
a = q(1); b = q(2); c = q(3); %d = q(4);
s_1 = s(s < s_split);
% S_1 = a * (s_1 + 0.005).^b + c * (s_1 + 0.005).^d;
% 
% % Calculate value, gradient and curvature at join
% S_split = a * (s_split + 0.005)^b + c * (s_split + 0.005)^d;
% dSds_split = a*b*(s_split+0.005)^(b-1) + c*d*(s_split+0.005)^(d-1);
% d2Sds2_split = a*b*(b-1)*(s_split+0.005)^(b-2) + c*d*(d-1)*(s_split+0.005)^(d-2);

s_shift = 0.005;
S_1 = a*(s_1 + 0.005).^b + c*s_1.^2;
S_split = a * (s_split + s_shift)^b + c * s_split^2;
dSds_split = a*b*(s_split+s_shift)^(b-1) + 2*c*s_split;
d2Sds2_split = a*b*(b-1)*(s_split+s_shift)^(b-2) + c;

% Calculate second region of shape distribution by cubic fit
a = q(4);
b = 0.5 * d2Sds2_split - 3*a*s_split;
c = dSds_split - 3*a*s_split^2 - 2*b*s_split;
d = S_split - a*s_split^3 - b*s_split^2 - c*s_split;
s_2 = s(s >= s_split);
S_2 = a*s_2.^3 + b*s_2.^2 + c*s_2 + d;

% Calculate error
% F = sum([S_1 ; S_2] - S).^2;

% Return new shape shape distribution
S = [S_1 ; S_2];

end
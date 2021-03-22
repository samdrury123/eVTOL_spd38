function t = bl_fit_camber(s_cl,cam,chi_le,chi_te,plot_stuff)
% BL_FIT_CAMBER  Fit a dimensional camber distribution smoothly
%
%   t = bl_fit_camber(s,cam,plot_stuff)
%  
%   s - non-dimensional coordinate along camberline
%   cam - non-dimensional camber value
%   chi_le - leading edge metal angle
%   chi_te - trailing edge metal angle
%   order - order of polynomial 3 for cubic or 4 for quartic
%   plot_stuff - 0 or 1 for showing working
%   t - fitted camber parameters

% Default to not show plots
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Plot original camber distribution
if plot_stuff == 1
    figure(); hold on; grid on; box on;
    plot(s_cl,cam,'r-');
end

% Fit with quartics
order = 4;
    
% Interpolate points uniformly along the non-dimensional camber line
r = dist_2d([s_cl cam],1); ni = 100; r_fit = linspace(0,1,ni); 
s_fit = interp1(r,s_cl,r_fit,'pchip'); cam_fit = interp1(r,cam,r_fit,'pchip');

% Fit a polynomial through the camber distribution
p = polyfit(s_fit,cam_fit,order); cam_fit = polyval(p,s_fit);

% Rescale camber distribution with new metal angles
chi = cam_fit * (chi_le - chi_te) + chi_te;
cam_new = (chi - chi(end)) / (chi(1) - chi(end));
p = polyfit(s_fit,cam_new,order);

% Calculate new camber distribution
t.s_cl = s_cl; t.cam = polyval(p,s_cl);

% Calculate new camber parameters
t.chi_le = chi(1); t.chi_te = chi(end);
t.dcam_le = polyval(polyder(p),0); t.dcam_te = polyval(polyder(p),1);
t.qcam = p(1) - 2 * p(2) + p(3);

% Plot fitted camberline
if plot_stuff == 1
    plot(t.s_cl,t.cam * (t.chi_le - t.chi_te) + t.chi_te,'b-');
    legend('Original','Fit');
end

% Record polynomial coefficients for non-dimensional camberline
t.p = p;

end
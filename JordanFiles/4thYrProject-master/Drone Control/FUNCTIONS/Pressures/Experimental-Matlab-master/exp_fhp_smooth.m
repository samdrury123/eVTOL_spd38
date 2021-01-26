function c = exp_fhp_smooth(c,plot_stuff)
% EXP_FHP_SMOOTH  Smooth and interpolate a calibration grid to improve results

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Plot original calibration map
if plot_stuff == 1
    cols = lines(2);
    h = exp_fhp_plot(c,[],cols(1,:));
end

% Two dimensional cubic interpolation of calibration variables
varnames = {'Iota' 'Tau' 'C_Po' 'C_P' 'C_Alpha' 'C_Beta'}; n_divide = 2;
for v = 1:length(varnames)
    c.(varnames{v}) = interp2(c.(varnames{v}),n_divide,'cubic');
end

% Smooth calibration variables
varnames = {'C_Po' 'C_P' 'C_Alpha' 'C_Beta'};
for v = 1:length(varnames)
    c.(varnames{v}) = smooth3(repmat(c.(varnames{v}),[1 1 2]),'box',[7 7 3]);
    c.(varnames{v}) = c.(varnames{v})(:,:,1);
end

% Plot new calibration map
if plot_stuff == 1
    exp_fhp_plot(c,h,cols(2,:));
end

end
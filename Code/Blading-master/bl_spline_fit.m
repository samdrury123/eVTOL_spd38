function b = bl_spline_fit(b)
% BL_SPLINE_FIT  Fit splines through numeric values of blade parameters

% Specify variables to evaluate
varnames = {'chi_le' 'dcam_le' 'dcam_te' 'qcam' 'chi_te' 'tchord'...
    'rad_le' 'thick_max' 's_thick_max' 'rad_thick_max' 'thick_te' 'wedge_te'...
    'sweep' 'lean'};

% Loop over all variables and fit numeric specifications
for v = 1:length(varnames)
    if isnumeric(b.(varnames{v}))
        
        % Find resolution of specification
        nj = length(b.(varnames{v}));
        
        % Fit spline reducing order by one point
        if isvector(b.(varnames{v})) == 0
            r_nondim = linspace(0,1,100);
            val = interp1(b.(varnames{v})(:,1),b.(varnames{v})(:,2),r_nondim,'pchip');
%             b.(varnames{v}) = spap2(nj-3,4,b.(varnames{v})(:,1),b.(varnames{v})(:,2));
            b.(varnames{v}) = spap2(nj-3,4,r_nondim,val);
        else
            b.(varnames{v}) = spap2(nj-3,4,linspace(0,1,nj),b.(varnames{v}));
        end
    end
end

end
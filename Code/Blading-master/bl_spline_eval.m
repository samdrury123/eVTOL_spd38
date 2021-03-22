function c = bl_spline_eval(b,r_nondim)
% BL_SPLINE_EVAL  Evaluate a spline definition of a blade at desired resolution

% Specify variables to evaluate
varnames = {'chi_le' 'dcam_le' 'dcam_te' 'qcam' 'chi_te' 'tchord'...
    'rad_le' 'thick_max' 's_thick_max' 'rad_thick_max' 'thick_te' 'wedge_te'...
    'sweep' 'lean' 'cam_md' 'chord' 'fchord'};

% Loop over all variables, extrapolate and evaluate spline
for v = 1:length(varnames)
    if isfield(b,varnames{v}) == 1
        if isstruct(b.(varnames{v})) == 1
            c.(varnames{v}) = fnval(fnxtr(b.(varnames{v}),b.(varnames{v}).order),r_nondim);
        elseif size(b.(varnames{v}),2) == 1
            if numel(b.(varnames{v})) == 1; b.(varnames{v}) = repmat(b.(varnames{v}),[2 1]); end
            c.(varnames{v}) = interp1(linspace(0,1,size(b.(varnames{v}),1)),b.(varnames{v}),...
                r_nondim','spline');
        elseif size(b.(varnames{v}),2) == 2
            c.(varnames{v}) = interp1(b.(varnames{v})(:,1),b.(varnames{v})(:,2),r_nondim','spline');
        end
    end
end

% Evaluate qcam for old definitions
if isfield(c,'qcam') == 0
    c.qcam = zeros(size(c.cam_md));
    for n = 1:length(r_nondim)
        p = [0 0 0 0 1 ; 1 1 1 1 1 ; 0 0 0 1 0 ; 4 3 2 1 0 ; 0.0625 0.125 0.25 0.5 1] \ ...
            [1 ; 0 ; c.dcam_le(n) ; c.dcam_te(n) ; c.cam_md(n)];
        c.qcam(n) = p(1) - 2 * p(2) + p(3);
    end
end

% Approximate true chord for old definitions
if isfield(c,'tchord') == 0
    c.tchord = c.chord;
end

% Reference coordinates for section
if isfield(b,'x_ref') == 1
    c.xy_cen = [b.x_ref 0];
end

% Record non-dimensional radius
c.r_nondim = r_nondim;

% Record non-dimensional axial coordinate
if isfield(b,'x_nondim') == 1
    c.x_nondim = b.x_nondim; 
end

% Calculate meridional parameters
if isfield(b,'r_le') == 1

    % Calculate hub and casing radii from mean radii and areas
    r_hub_le = b.r_le - b.A_le / (4 * pi * b.r_le);
    r_cas_le = b.r_le + b.A_le / (4 * pi * b.r_le);
    r_hub_te = b.r_te - b.A_te / (4 * pi * b.r_te);
    r_cas_te = b.r_te + b.A_te / (4 * pi * b.r_te);

    % Radii at leading and trailing edge
    c.r_le = r_nondim * (r_cas_le - r_hub_le) + r_hub_le;
    c.r_te = r_nondim * (r_cas_te - r_hub_te) + r_hub_te;
    
end

% Blade number
if isfield(b,'N') == 1
    c.N = b.N;
end

% Leading edge type
if isfield(b,'le') == 1
    c.le = b.le;
end

end
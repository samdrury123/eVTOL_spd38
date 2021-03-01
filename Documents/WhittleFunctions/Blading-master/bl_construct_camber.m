function t = bl_construct_camber(c,s_cl)
% BL_CONSTRUCT_CAMBER  Calculate a camberline from given camber parameters
%
%   t = BL_CONSTRUCT_CAMBER(c,s_cl)
%
%   c - input struct containing thickness parameters
%   s_cl - optional argument to specify chordwise points
%
%   t - output struct containing camber distribution
%
%   c must contain fields:
%       dcam_le - gradient of non-dimensional camberline at the leading edge
%       dcam_te - gradient of non-dimensional camberline at the trailing edge
%       qcam - parameter describing the kinky-ness of the camberline
%       chi_le - leading edge metal angle
%       chi_te - trailing edge metal angle
%       tchord - true chord distance along the camber line

% Default chordwise point distribution
if exist('s_cl','var') == 0
    ni = 1101; s_cl = linspace(-0.05,1.05,ni).';
end

%% Non-dimensional camber quartic by inverting matrix of boundary conditions

% Preallocate vectors and arrays
A = zeros(5,5); y = zeros(5,1);

% Y-intercept is 1
A(1,:) = [0 0 0 0 1]; y(1) = 1;

% X-intercept is 1
A(2,:) = [1 1 1 1 1]; y(2) = 0;

% Kinky-ness of the camberline
A(3,:) = [1 -2 1 0 0]; y(3) = c.qcam;

% Leading edge gradient
A(4,:) = [0 0 0 1 0]; y(4) = c.dcam_le;

% Trailing edge gradient
A(5,:) = [4 3 2 1 0]; y(5) = c.dcam_te;

% Invert to find polynomial coefficients
p = A \ y;


%% Evaluate and construct polynomial camber line

% Evaluate camberline at desired chordwise points
cam = polyval(p,s_cl);

% Dimensionalise camberline with inlet and outlet metal angles
chi = cam * (c.chi_le - c.chi_te) + c.chi_te;

% Find leading and trailing edge points
[~, i_le] = min(abs(s_cl - 0)); [~, i_te] = min(abs(s_cl - 1));

% Construct x-rt coordinates of camberline
x = [0 ; cumsum(0.5*(cosd(chi(2:end)) + cosd(chi(1:end-1))) .* diff(s_cl))]; 
rt = [0 ; cumsum(0.5*(sind(chi(2:end)) + sind(chi(1:end-1))) .* diff(s_cl))];

% Scale camberline by desired true chord
sf = c.tchord / sum((diff(x(i_le:i_te),1,1).^2 + diff(rt(i_le:i_te),1,1).^2).^0.5);
xy_cam = [x - x(i_le) rt - rt(i_le)] * sf; 

% Coordinates on straight chord line
xy_chord = [s_cl * (x(i_te) - x(i_le)) + x(i_le) s_cl * (rt(i_te) - rt(i_le)) + rt(i_le)] * sf;

% Record output
t.s = s_cl; t.cam = cam; t.chi = chi; t.xy_cam = xy_cam; t.xy_chord = xy_chord;
t.i_le = i_le; t.i_te = i_te;


end


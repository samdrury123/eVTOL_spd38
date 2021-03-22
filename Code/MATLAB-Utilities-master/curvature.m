function [k,varargout] = curvature(x,y,varargin)
% CURVATURE
%
%   Calculates the curvature of a curve in cartesian coords.
%
% Format:
%
%   "k = curvature(x,y);"
%


%% Defaults/ Program Settings

defaults.isClosed = false;
defaults.tol = 1e-6;
defaults.endVal = NaN;

%% Input Arguments

par = inputParser;
addRequired(par,'x',@isnumeric);
addRequired(par,'y',@isnumeric);
fields = fieldnames(defaults);
for ii = 1:length(fields)
    addParamValue(par,fields{ii},defaults.(fields{ii}));
end
parse(par,x,y,varargin{:});

% Housekeeping
[nr,nc] = size(x);

%% Surface Distance

dx = x(2:end,:) - x(1:end-1,:);
dy = y(2:end,:) - y(1:end-1,:);
ds = sqrt(dx.^2 + dy.^2);
s  = zeros(nr,nc);
s(2:end,:) = cumsum(ds,1);

%% Manages Closed Shapes

% A line is treated as closed if the start and end points are duplicated
isClosed = all( s(1,:)-s(end,:) < par.Results.tol );

% Duplicates the 2nd & penultimate points so that the curvature of the
% end-points is calculated.
if isClosed
    x = [x(end-1,:); x; x(2,:)];
    y = [y(end-1,:); y; y(2,:)];
    s = [s(end-1,:)-s(end); s; s(2)+s(end)];
end

%% Calculates Curvature

% Gradients of the bisectors
dx1 =  ( y(2:end-1,:) - y(1:end-2,:) );
dy1 = -( x(2:end-1,:) - x(1:end-2,:) );
dx2 =  ( y(3:end  ,:) - y(2:end-1,:) );
dy2 = -( x(3:end  ,:) - x(2:end-1,:) );

% Bisector points
x1 = 0.5.*( x(2:end-1,:) + x(1:end-2,:) );
y1 = 0.5.*( y(2:end-1,:) + y(1:end-2,:) );
x2 = 0.5.*( x(3:end  ,:) + x(2:end-1,:) );
y2 = 0.5.*( y(3:end  ,:) + y(2:end-1,:) );

% Finds the centre of the circle by solving matrix equations
det = -dx1.*dy2 + dx2.*dy1;
a11 = -dy2./det;
a12 =  dx2./det;
a21 = -dy1./det;
a22 =  dx1./det;
lam1 = a11.*(x2-x1) + a12.*(y2-y1);
lam2 = a21.*(x2-x1) + a22.*(y2-y1);
x0 = x1 + lam1.*dx1;
y0 = y1 + lam1.*dy1;

% Calculates curvature
dx = x0 - x(2:end-1);
dy = y0 - y(2:end-1);
r = sqrt(dx.^2 + dy.^2);
k = 1./r;

%% End points

if ~isClosed
    endVal = par.Results.endVal.*ones(1,nc);
    k  = [ endVal;  k; endVal];
    x0 = [ endVal; x0; endVal];
    y0 = [ endVal; y0; endVal];
else
    s = s(2:end-1,:);
end


%% Output Arguments

if nargout >= 2; varargout{1} = s ; end
if nargout >= 3; varargout{2} = x0; end
if nargout >= 4; varargout{3} = y0; end


%% ends function CURVATURE
end


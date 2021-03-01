function varargout = smoothSect(x,y);
%SMOOTHSECT

%program defaults
newQuadTolFac = 0.8;
curvExp = 1.0;
dsCurvMax = 0.1;
nNew = 500;

%removes duplicate points (i.e. repeated start/end points)
tmp = [x,y];
tmp = unique(tmp,'rows','stable');
x = tmp(:,1);
y = tmp(:,2);

%adds duplicated points back in - more robust this way (handles both closed
%and open polygons.
x = [x; x(1)];
y = [y; y(1)];


x = smooth(x,5);
y = smooth(y,5);

%distance vector
dx = x(2:end) - x(1:end-1);
dy = y(2:end) - y(1:end-1);
ds = sqrt( dx.*dx + dy.*dy );
s = [0; cumsum(ds)];

%converts to intrinsic coordinates
xtmp = [x(end-1); x; x(2)];
ytmp = [y(end-1); y; y(2)];
dx = xtmp(3:end) - xtmp(1:end-2);
dy = ytmp(3:end) - ytmp(1:end-2);
phi = atan2( dy, dx );

%finds axis crossing points and makes phi continuous.
npts = length(phi);
dphi = phi(2:end) - phi(1:end-1);
phiSign = sign( phi(2:end) .* phi(1:end-1) );
isNewQuad = phiSign < 0 & abs(dphi) > 2.*pi.*newQuadTolFac;
indNewQuad = find(isNewQuad)+1;
for ii = 1:length(indNewQuad)
    iiNewQuad = [indNewQuad(ii):npts]';
    phi(iiNewQuad) = phi(iiNewQuad) + 2.*pi;
end
phiAdded = 2.*pi.*length(indNewQuad);

%calculates curvatures
ds1 = s(end-1) - s(end);
ds2 = s(end)   + s(2);
stmp = [ ds1;   s;  ds2];
phi1 = phi(end-1) - phiAdded; 
phi2 = phi(2)     + phiAdded;
ptmp = [phi1; phi; phi2];
ds   = stmp(3:end) - stmp(1:end-2);
dphi = ptmp(3:end) - ptmp(1:end-2);
curv = dphi./ds;


%reinterpolates blade
% dphi = abs(phi(2:end) - phi(1:end-1));
% dphi = dphi.^curvExp;
% phiTmp = [0;cumsum(dphi)];
% phiTmp = phiTmp./phiTmp(end);
% phiNew = [0:nNew-1]'./(nNew-1);
% pNew = interp1( phiTmp, phi, phiNew,'cubic');
% sNew = interp1( phiTmp,   s, phiNew,'cubic');
% xNew = interp1( phiTmp,   x, phiNew,'cubic');
% yNew = interp1( phiTmp,   y, phiNew,'cubic');

%recalculates curvature
% phiTmp = [pNew(end); pNew; pNew(1)];
% sTmp = [ sNew(end-1)-sNew(end); sNew; sNew(end) + sNew(2)];
% dphi = phiTmp(3:end) - phiTmp(1:end-2);
% ds   =   sTmp(3:end) -   sTmp(1:end-2);
% cNew = dphi./ds;

sNew = [0:nNew-1]'./(nNew-1);
sTmp = s./s(end);
pNew = interp1( sTmp, phi, sNew,'cubic');
xNew = interp1( sTmp,   x, sNew,'cubic');
yNew = interp1( sTmp,   y, sNew,'cubic');
sNew = interp1( sTmp,   s, sNew,'cubic');

if nargout >= 1
    varargout{1} = [xNew,yNew];
end
if nargout >= 2
    varargout{1} = xNew;
    varargout{2} = yNew;
end


% 
subplot(2,2,1)
plot(s,curv,'bx-');grid on

subplot(2,2,2)
plot(xNew,yNew,'r+-',x,y,'bx-');
axis equal
grid on

subplot(2,2,3)
plot(s,phi,'bx-');
xlabel('s')
ylabel('\phi')
grid on

subplot(2,2,4)
plot(cumsum(abs(curv).^0.5),s,'b-',curv,s,'r-');
grid on
% 
% 
% 
% %plots
% 
% 
% 

%ends function smooth sect
end
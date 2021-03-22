function xyz = cutMesh(wrl,cutDir,cutVal)
%CUTMESH

%defaults
dsLim = 0.025;

%generates cutting plane
switch lower(cutDir)
    case {'x'}
        plane.origin = [cutVal, 0, 0];
        plane.vect1 = [0, 1, 0];
        plane.vect2 = [0, 0, 1];
        
    case {'y'}
        plane.origin = [0, cutVal, 0];
        plane.vect1 = [1, 0, 0];
        plane.vect2 = [0, 0, 1];       
    
    case {'z'}
        plane.origin = [0,0 ,cutVal];
        plane.vect1 = [1, 0, 0];
        plane.vect2 = [0, 1, 0];     
    
    otherwise
        errstr = ['Cut direction "',cutDir,'" not recognised.'];
        error(errstr)
end

%unpacks x,y & z coords
x1 = wrl.pts(:,1);
y1 = wrl.pts(:,2);
z1 = wrl.pts(:,3);
x1 = x1(wrl.knx);
y1 = y1(wrl.knx);
z1 = z1(wrl.knx);

x2 = [x1(:,2:end),x1(:,1)];
y2 = [y1(:,2:end),y1(:,1)];
z2 = [z1(:,2:end),z1(:,1)];

%calculates distance from cut plane
s1 = distancePointPlane( wrl.pts, [plane.origin, plane.vect1, plane.vect2]);
s1 = s1( wrl.knx);
s2 = [s1(:,2:end),s1(:,1)];

%look for sign changes
cutPlane = sign( s1 .* s2);
cutPlane = cutPlane <= 0;

%only keeps vectors which cut the specified plane
x1  = x1(cutPlane);
y1  = y1(cutPlane);
z1  = z1(cutPlane);
s1  = s1(cutPlane);
x2  = x2(cutPlane);
y2  = y2(cutPlane);
z2  = z2(cutPlane);
s2  = s2(cutPlane);

%interpolates to intersection points
ds = s2 - s1;
dxds = (x2-x1)./ds;
dyds = (y2-y1)./ds;
dzds = (z2-z1)./ds;
xyz(:,1) = -dxds.*s1 + x1;
xyz(:,2) = -dyds.*s1 + y1;
xyz(:,3) = -dzds.*s1 + z1;

%removes duplicate points
xyz = unique(xyz,'rows');

%reorders points using a shortest distance method
[npts,nc] = size(xyz);
xyz0 = xyz(1,:);
xyzTmp = xyz(2:end,:);
iipt = 1;
for ii = 1:npts-1
    %finds least distance
    dxyz = xyzTmp - repmat( xyz0, npts-ii, 1);
    ds = sqrt( sum(dxyz.*dxyz,2) );
    [dsMin,iiMin] = min(ds);

    %updates coordinate arrays (only if min dist is > limit set above).
    if dsMin >= dsLim
        iipt = iipt + 1;
        xyz(iipt,:) = xyzTmp(iiMin,:);
        xyz0 = xyzTmp(iiMin,:);
    end
    
    %removes point from temp coordinate array
    iiGtMin = [1:iiMin-1,iiMin+1:npts-ii]';
    xyzTmp = xyzTmp(iiGtMin,:);
end
xyz = xyz(1:iipt,:);
xyz = [xyz;xyz(1,:)];

%smooths points
% for ii = 1:3
%     xyz(:,ii) = smooth(xyz(:,ii),5);
% end

%ends function CUTMESH
end
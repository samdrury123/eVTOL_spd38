function g = ts_add_duct(g)
% Add an upstream duct to a blade row to improve the inlet boundary
% conditions

ER = 1.2;

% Get inlet data
for bid = ts_get_bids(g)
    for pid = ts_get_pids(g{bid+1})
        if g{bid+1}.patch{pid+1}.attribute.kind == 0
            bid_in = bid;
            pid_in = pid;
        end
    end
end

% Determine coordinates and spacings on inlet plane
x2 = g{bid_in+1}.x(1,:,:);
r2 = g{bid_in+1}.r(1,:,:);
t2 = g{bid_in+1}.rt(1,:,:) ./ r2;
dx2 = g{bid_in+1}.x(2,:,:) - x2;
dri2 = g{bid_in+1}.r(2,:,:) - r2;
drj2 = diff(r2(:,2));

% Determine length required to straighten the inlet
l = 2 * (max(max(x2)) - min(min(x2)));

% Determine resolution
x1 = max(max(x2)) - l;
dx_max = max(max(x2)) - x1;
dx2_min = min(min(dx2));

B = 1 - dx_max * (1 - ER) / dx2_min;
ni = double(ceil(log(B) / log(ER)));

if ni < 9
    ni = 9;
    l = 1.5 * l;
    x1 = max(max(x2)) - l;
end

nj = g{bid_in+1}.attribute.nj;
nk = g{bid_in+1}.attribute.nk;

% Determine spacings
sr = repmat( ( r2-repmat(min(r2,[],2),[1 nj 1]) ) ./ (repmat(max(r2,[],2),[1 nj 1]) - repmat(min(r2,[],2),[1 nj 1])), [ni 1 1]);

% Calculate coordinates
% Axial coordinates
era = 0.8 * ones(1,nj,nk);
erb = 1.5 * ones(1,nj,nk);
error = ones(1,nj,nk);

% Get expansion ratio by bisection method
while max(max(error)) > 1e-9
    er = 0.5 * (era + erb);
    dx = dx2 .* (1 - er.^(ni-1)) ./ (1 - er);
    diff_dx = dx - (x2 - x1);

    era(diff_dx < 0) = er(diff_dx < 0);
    erb(diff_dx > 0) = er(diff_dx > 0);
    error = abs(diff_dx);
end

% Calculate x coordinates from expansion ratios
n = repmat(reshape(ni-1:-1:0,[ni 1 1]),[1 nj nk]);
er = repmat(er, [ni 1 1]);
dx = repmat(dx2,[ni 1 1]) .* (1 - er.^n) ./ (1 - er);
x = repmat(x2,[ni 1 1]) - dx;

% Rescale x for errors in the geometric progression
% x = (x - repmat(x(end,:,:),[ni 1 1])) .* (repmat(x(end,:,:),[ni 1 1]) - x1) ./ repmat(x(end,:,:) - x(1,:,:),[ni 1 1]) + repmat(x(end,:,:),[ni 1 1]);

% Radial Coordinates
% Define exponential hub and casing lines
x2h = x2(1,1,:); r2h = r2(1,1,:); 
x2c = x2(1,end,:); r2c = r2(1,end,:); 
drdx2h = dri2(1,1,:) ./ dx2(1,1,:);
drdx2c = dri2(1,end,:) ./ dx2(1,end,:);

% Set scaling and calculate constants for casing
ac = 0.1 * (r2c - r2h);
y1c = 12;

B = abs(drdx2c .* (x2c - x1) ./ ac);

y2 = -1:0.0001:y1c;
f = (y1c(1) - y2) .* exp(-y2) - B(1);
f(abs(f) < -10) = -10;

[~, i] = min(f(1:end-1).*f(2:end));
y2c = y2(i);

% Calculate r coordinates on the casing
xc = x(:,end,:);
ac = repmat(ac,[ni 1 1]);
x2c = repmat(x2c,[ni 1 1]);
y2c = repmat(y2c,[ni 1 nk]);
r2c = repmat(r2c,[ni 1 1]);

if min(drdx2c) < 0
    z2c = (1 - exp(-y2c));
    rc = ac .* (1 - exp(-((x2c - xc).*(y1c - y2c)./(x2c - x1) + y2c)) - z2c) + r2c;
else
    z2c = (exp(-y2c) - 1);
    rc = ac .* (exp(-((x2c - xc).*(y1c - y2c)./(x2c - x1) + y2c)) - 1 - z2c) + r2c;
end    

% Set scaling and calculate constants for the hub
r2c = r2(1,end,:); 
ah = 0.1 * (r2c - r2h);
y1h = 12;

B = abs(drdx2h .* (x2h - x1) ./ ah);

y2 = -1:0.0001:y1h;
f = (y1h(1) - y2) .* exp(-y2) - B(1);
f(abs(f) < -10) = -10;

[~, i] = min(f(1:end-1).*f(2:end));
y2h = y2(i);

% Calculate r coordinates on the casing
xh = x(:,1,:);
ah = repmat(ah,[ni 1 1]);
x2h = repmat(x2h,[ni 1 1]);
y2h = repmat(y2h,[ni 1 nk]);
r2h = repmat(r2h,[ni 1 1]);

if min(drdx2h) < 0
    z2h = (1 - exp(-y2h));
    rh = ah .* (1 - exp(-((x2h - xh).*(y1h - y2h)./(x2h - x1) + y2h)) - z2h) + r2h;
else
    z2h = (exp(-y2h) - 1);
    rh = ah .* (exp(-((x2h - xh).*(y1h - y2h)./(x2h - x1) + y2h)) - 1 - z2h) + r2h;
end    

% Calculate r coordinates in whole block
r = sr .* (repmat(rc,[1 nj 1]) - repmat(rh,[1 nj 1])) + repmat(rh,[1 nj 1]);

% Tangential coordinates
t = repmat(t2,[ni 1 1]);

% Patch on to end of ts data structure
% Block properties
bid_duct = length(g);
b.attribute.bid = bid_duct;
b.attribute.ni = ni;
b.attribute.nj = nj;
b.attribute.nk = nk;
b.attribute.np = 4;
b.attribute.nc = 0;
b.attribute.ncl = 0;
b.attribute.procid = g{bid_in+1}.attribute.procid;
b.attribute.threadid = 0;
b.x = x;
b.r = r;
b.rt = r .* t;
b.bv = g{bid_in+1}.bv;
b.av = g{1}.av;

g = [g ; b];

% Set patch attributes and properties
% Move inlet patch
g{end}.patch{1} = g{bid_in+1}.patch{pid_in+1};
g{end}.patch{1}.attribute.bid = bid_duct;
g{end}.patch{1}.attribute.pid = 0;

% Set periodic to old inlet patch
g{bid_in+1}.patch{pid_in+1}.attribute.kind = 5;
g{bid_in+1}.patch{pid_in+1}.attribute.idir = 0;
g{bid_in+1}.patch{pid_in+1}.attribute.jdir = 1;
g{bid_in+1}.patch{pid_in+1}.attribute.kdir = 2;
g{bid_in+1}.patch{pid_in+1}.attribute.nxbid = bid_duct;
g{bid_in+1}.patch{pid_in+1}.attribute.nxpid = 1;
g{bid_in+1}.patch{pid_in+1} = rmfield(g{bid_in+1}.patch{pid_in+1},...
    {'pitch' 'pstag' 'rfin' 'sfinlet' 'tstag' 'yaw'});

% Periodic on new block
p1.bid = bid_duct;
p1.pid = 1;
p1.kind = 5;
p1.ist = ni-1;
p1.ien = ni;
p1.jst = 0;
p1.jen = nj;
p1.kst = 0;
p1.ken = nk;
p1.idir = 0;
p1.jdir = 1;
p1.kdir = 2;
p1.nxpid = pid_in;
p1.nxbid = bid_in;
g{end}.patch{2,1}.attribute = p1;

% Periodic patches accross the duct block
bid0 = bid_duct;
bid1 = bid_duct;

Q0 = [0 g{bid0+1}.attribute.ni ; 0 g{bid0+1}.attribute.nj ;...
    0 1];
dir0 = [0 1 2];
Q1 = [0 g{bid1+1}.attribute.ni ; 0 g{bid1+1}.attribute.nj ; g{bid1+1}.attribute.nk-1 g{bid1+1}.attribute.nk];
dir1 = [0 1 2];

g = ts_add_periodic_patch(bid0, bid1, Q0, Q1, dir0, dir1, g, 5);

% Set initial guess
if isfield(g{bid_in+1},'rovx') == 1
    g{bid_duct+1}.ro = repmat(g{bid_in+1}.ro(1,:,:),[ni 1 1]);
    g{bid_duct+1}.rovx = repmat(g{bid_in+1}.rovx(1,:,:),[ni 1 1]);
    g{bid_duct+1}.rovr = repmat(g{bid_in+1}.rovr(1,:,:),[ni 1 1]);
    g{bid_duct+1}.rorvt = repmat(g{bid_in+1}.rorvt(1,:,:),[ni 1 1]);
    g{bid_duct+1}.roe = repmat(g{bid_in+1}.roe(1,:,:),[ni 1 1]);
    g{bid_duct+1}.trans_dyn_vis = repmat(g{bid_in+1}.trans_dyn_vis(1,:,:),[ni 1 1]);
end

g{bid_duct+1}.bv.tstagin = g{bid_in+1}.bv.tstagin;
g{bid_duct+1}.bv.tstagout = g{bid_in+1}.bv.tstagout;
g{bid_duct+1}.bv.pstatin = g{bid_in+1}.bv.pstatin;
g{bid_duct+1}.bv.pstatout = g{bid_in+1}.bv.pstatout;
g{bid_duct+1}.bv.vgridin = g{bid_in+1}.bv.vgridin;
g{bid_duct+1}.bv.vgridout = g{bid_in+1}.bv.vgridout;

end
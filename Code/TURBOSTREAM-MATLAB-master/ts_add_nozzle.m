function g = ts_add_nozzle(g, M_2, length_ratio, mixing_plane, outlet, AR)
% Generates a downstream nozzle mesh

if exist('length_ratio','var') == 0
    length_ratio = 1;
end

if exist('mixing_plane','var') == 0
    mixing_plane = 1;
end

% Set defaults
if exist('AR','var') == 0
    AR = 400; % Aspect ratio at exit
end
ER = 1.2; % Expansion ratio

% Set area ratio
for bid = ts_get_bids(g)
    for pid = ts_get_pids(g{bid+1})
        if g{bid+1}.patch{pid+1}.attribute.kind == 1
            bid_exit = bid;
            pid_exit = pid;
        end
    end
end
if exist('outlet','var') == 0 || isempty(outlet) == 1
    % Get exit data
    ien = g{bid_exit+1}.patch{pid_exit+1}.attribute.ien;
    P_1 = g{bid_exit+1}.patch{pid_exit+1}.pout;
    outlet = ts_secondary(ts_structured_cut(g,bid_exit,ien,ien,1,'en',1,'en')); 
    M_1 = ts_mass_average(outlet,'M',3); Alpha_1 = ts_mass_average(outlet,'Alpha',3);
    M_1 = M_1 * cosd(Alpha_1);
    r_an = [g{bid_exit+1}.r(end,1,1) g{bid_exit+1}.r(end,end,1)];
else
    M_1 = outlet.M;
    P_1 = outlet.P;
    r_an = outlet.r_an;
    if isfield(outlet,'area_ratio') == 1
        area_ratio = outlet.area_ratio;
    end
end

% Determine area ratio and back pressure from compressible relations
if exist('area_ratio','var') == 0
    M_func_2 = 1.4/(0.4^0.5) * M_2 * (1 + 0.2 * M_2^2)^(-3);
    M_func_1 = 1.4/(0.4^0.5) * M_1 * (1 + 0.2 * M_1^2)^(-3);
    area_ratio = M_func_1 / M_func_2;
end

M_stat_2 = 1.4/(0.4^0.5) * M_2 * (1 + 0.2 * M_2^2)^0.5;
M_stat_1 = 1.4/(0.4^0.5) * M_1 * (1 + 0.2 * M_1^2)^0.5;
Pr = M_stat_1 / (M_stat_2 * area_ratio);
P_2 = Pr * P_1;
pout = P_2;

% Blade chord and span to set nozzle length
chords = ts_find_chord(g,r_an); chord = chords(end);
span = r_an(2) - r_an(1);

%==========================================================================
% Resolution, and geometrical parameters for the nozzle
%==========================================================================

% Upstream of exit patch
x0 = squeeze(g{bid_exit+1}.x(end-1,:,:)); 
r0 = squeeze(g{bid_exit+1}.r(end-1,:,:)); 

% Exit patch
x1 = squeeze(g{bid_exit+1}.x(end,:,:)); 
r1 = squeeze(g{bid_exit+1}.r(end,:,:)); 
rt1 = squeeze(g{bid_exit+1}.rt(end,:,:));
t1 = rt1 ./ r1;

% Nozzle exit
% length_total = 0.3*D*length_ratio;
% length_total = length_ratio * chord;
length_total = length_ratio * (chord * span)^0.5;
x2 = length_total + min(min(x1));

% Determine number of axial points in the nozzle mesh
% Section based on maximum expansion ratio then a section with maximum
% aspect ratio in order to minimise ni
dri1 = r1 - r0;
drj1 = diff(r1(:,1));
dx1 = x1 - x0;
dx1_min = min(dx1(:,1));
x1_min = min(x1(:,1));
x1_max = max(x1(:,1));
x_max = max(max(x2 - x1));

dr2 = drj1 * area_ratio^0.5;
dx2_max = min(dr2) * AR;
dx2 = dx2_max;

% Points in expanding section
ni1 = double(ceil((log(dx2_max) - log(dx1_min))/log(ER)));

% Points in uniform section
% x12 = dx1_min * (1 - ER^ni1+1) / (1 - ER) + x1_min;
x12 = dx1_min * (1 - ER^ni1+1) / (1 - ER) + x1_max;

if ni1 < 4
    ni1 = 1;
    x12 = max(max(x1));
end

ni2 = double(ceil((x2 - x12) / dx2_max));

ni = ni1 + ni2;

nj = g{bid_exit+1}.attribute.nj;

if mixing_plane == 1
    nk = 9;
    dx1 = repmat(dx1(:,1),1,nk);
    x1 = repmat(x1(:,1),1,nk);
    drj1 = repmat(drj1(:,1),1,nk);
    dri1 = repmat(dri1(:,1),1,nk);
    r1 = repmat(r1(:,1),1,nk);
else
    nk = g{bid_exit+1}.attribute.nk;
end

% Clustering of points in radial direction
sr = repmat(reshape((r1-repmat(min(r1),nj,1))./(repmat(max(r1),nj,1)-repmat(min(r1),nj,1)),1,nj,nk),[ni 1 1]);

%==========================================================================
% Extend the domain downstream
%==========================================================================

% Axial coordinates
dx2 = (x2 - x12) / ni2;
dx = x12 - x1;

% Find changes in spacing to determine boundary lines
d2x1 = diff(dx1(:,1),1);
[~,j_split] = max(d2x1);
J = [1 j_split nj];

if ni1 > 1
    er = repmat(0.5:0.00001:ER^2,[ni1 1]);
    n = repmat((0:ni1-1).',[1 size(er,2)]);
    A = zeros(length(J),1);
    for p = 1:length(J)
        j = J(p);
        f = sum(dx1(j,1) .* er.^(n-1), 1) - dx(j,1);
        f = f(1:end-1) .* f(2:end);
        [~, i] = min(f);
        A(p) = er(1,i);
    end

    % Determine clustering on boundary lines
    er = A;
    n = 0:ni1-1;
    xa = zeros(ni1,nj,nk);
    for p = 1:length(J)
        j = J(p);
        S = dx1(J(p),1)*(1 - er(p).^n)./(1 - er(p));
        S = S/max(S);
        xa(:,j,:) = repmat(reshape(S,[ni1 1 1]),[1 1 nk]) .* repmat(reshape(dx(j,:),[1 1 nk]),[ni1 1 1]) +...
            repmat(reshape(x1(j,:),[1 1 nk]),[ni1 1 1]);
    end
    

    % Fill in gaps with bilinear interpolation
    for p = 1:length(J)-1
        nj1 = J(p+1) - J(p) + 1;

        sp2 = (xa(:,J(p+1),:) - repmat(reshape(x1(J(p+1),:),[1 1 nk]),[ni1 1 1])) ./ repmat(reshape(dx(J(p+1),:),[1 1 nk]),[ni1 1 1]);
        sp1 = (xa(:,J(p),:) - repmat(reshape(x1(J(p),:),[1 1 nk]),[ni1 1 1])) ./ repmat(reshape(dx(J(p),:),[1 1 nk]),[ni1 1 1]);

        sx = repmat((sp2 - sp1),[1 nj1 1]) .* (sr(1:ni1,J(p):J(p+1),:) - repmat(sr(1:ni1,J(p),:),[1 nj1 1])) ./...
            (repmat(sr(1:ni1,J(p+1),:) - sr(1:ni1,J(p),:),[1 nj1 1])) + repmat(sp1,[1 nj1 1]);

        xa(:,J(p):J(p+1),:) = sx .* repmat(reshape(dx(J(p):J(p+1),:),[1 nj1 nk]),[ni1 1 1]) + repmat(reshape(x1(J(p):J(p+1),:),[1 nj1 nk]),[ni1 1 1]);
    end
else
    xa = x12 * ones(1,nj,nk);
end

xb = cumsum(dx2*ones(ni2,nj,nk),1) + x12;

x = cat(1, xa, xb);

% Radial coordinates - define hub & casing lines

% Gradients and known coordinates
drdx1h = repmat(reshape(dri1(1,:) ./ dx1(1,:),1,1,[]),[ni,1,1]);
drdx1c = repmat(reshape(dri1(end,:) ./ dx1(end,:),1,1,[]),[ni,1,1]); 
x1h = repmat(reshape(x1(1,:),1,1,[]),[ni,1,1]); x1c = repmat(reshape(x1(end,:),1,1,[]),[ni,1,1]); 
r1h = repmat(reshape(r1(1,:),1,1,[]),[ni,1,1]); r1c = repmat(reshape(r1(end,:),1,1,[]),[ni,1,1]); 

% Calculate scaling for exponential casing line
ac = 0.1 * (r1c - r1h);
y2c = 12;

B = abs(drdx1c .* (x2 - x1c) ./ ac);

y1 = -1:0.0001:y2c;
f = (y2c(1) - y1) .* exp(-y1) - B(1);
f(abs(f) < -10) = -10;

[~, i] = min(f(1:end-1).*f(2:end));
y1c = y1(i);

xc = x(:,end,:);

if min(drdx1c) > 0
    z1c = (1 - exp(-y1c));
    rc = ac .* (1 - exp(-((xc - x1c).*(y2c - y1c)./(x2 - x1c) + y1c)) - z1c) + r1c;
else
    z1c = (exp(-y1c) - 1);
    rc = ac .* (exp(-((xc - x1c).*(y2c - y1c)./(x2 - x1c) + y1c)) - 1 - z1c) + r1c;
end    

% Calculate scaling for tanh hub line
r2c = repmat(rc(end,1,:),[ni 1 1]);
r2h = (r2c.^2 - area_ratio*(r1c.^2 - r1h.^2)).^0.5;
y2h = 3; z2h = tanh(y2h);
B = drdx1h .* (x2 - x1h) ./ (r2h - r1h);

y1 = -20:0.0001:y2h;
f = (r2h(1) - r1h(1)) .* (y2h(1) - y1) .* ((sech(y1)).^2) ./ ((z2h(1) - tanh(y1)) .* (x2(1) - x1h(1)));
f(abs(f) > 10) = 10;
f = f - drdx1h(1);

[~, i] = min(f(1:end-1).*f(2:end));
y1h = y1(i);

if min(min(drdx1h)) < 0 || min(min(y1h)) < -y2h
    y1h = y1h./y1h * -y2h;
end
z1h = tanh(y1h);

xh = x(:,1,:);
rh = ((r2h - r1h) ./ (z2h - z1h)) .* (tanh( (xh - x1h).*(y2h-y1h)./(x2 - x1h) + y1h ) - z1h) + r1h;

% Fill in radial coordinates between hub and casing lines
r = sr .* repmat(rc - rh, [1 nj 1]) + repmat(rh, [1 nj 1]);

% Tangential coordinates
if mixing_plane == 1
    t = repmat(reshape(linspace(min(t1(1,:)),max(t1(1,:)),nk),1,1,nk),[ni nj 1]);
else
    t = repmat(reshape(t1,[1 nj nk]),[ni 1 1]);
end


%==========================================================================
% Patch on to end of TS data structure
%==========================================================================

% Add on the block
b.attribute.bid = -1;
b.attribute.ni = ni;
b.attribute.nj = nj;
b.attribute.nk = nk;
b.attribute.np = 4;
b.attribute.nc = 0;
b.attribute.ncl = 0;
b.attribute.procid = g{bid_exit+1}.attribute.procid;
b.attribute.threadid = 0;
b.x = x;
b.r = r;
b.rt = r .* t;
b.bv = g{bid_exit+1}.bv;
b.av = g{1}.av;

g = ts_add_block(b,g);

% Patch it appropriately
bid_nozzle = size(g,1)-1;
% Exit patch
g{end}.patch{1}.attribute.bid = bid_nozzle;
g{end}.patch{1}.attribute.pid = 0;
g{end}.patch{1}.attribute.kind = 1;
g{end}.patch{1}.attribute.ist = ni-1;
g{end}.patch{1}.attribute.ien = ni;
g{end}.patch{1}.attribute.jst = 0;
g{end}.patch{1}.attribute.jen = nj;
g{end}.patch{1}.attribute.kst = 0;
g{end}.patch{1}.attribute.ken = nk;

g{end}.patch{1}.pout = pout;
g{end}.patch{1}.ipout = 3;
g{end}.patch{1}.throttle_type = 0;
g{end}.patch{1}.throttle_target = 0;
g{end}.patch{1}.throttle_k0 = 0;
g{end}.patch{1}.throttle_k1 = 0;
g{end}.patch{1}.throttle_k2 = 0;

% Mixing plane or periodic patch to old exit block
bid0 = bid_exit;
pid0 = pid_exit;
bid1 = bid_nozzle;
pid1 = size(g{bid_nozzle+1}.patch,1);

p1.bid = bid1;
p1.pid = pid1;
if mixing_plane == 1
    p1.kind = 2;
else
    p1.kind = 5;
end
p1.ist = 0;
p1.ien = 1;
p1.jst = 0;
p1.jen = nj;
p1.kst = 0;
p1.ken = nk;
p1.idir = 0;
p1.jdir = 1;
p1.kdir = 2;
p1.nxpid = pid0;
p1.nxbid = bid0;

g{bid_nozzle+1}.patch{pid1+1,1}.attribute = p1;
g{bid_nozzle+1}.attribute.np = g{bid_nozzle+1}.attribute.np + 1;

p0.bid = bid0;
p0.pid = pid0;
if mixing_plane == 1
    p0.kind = 2;
else
    p0.kind = 5;
end
p0.ist = g{bid_exit+1}.attribute.ni-1;
p0.ien = g{bid_exit+1}.attribute.ni;
p0.jst = 0;
p0.jen = nj;
p0.kst = 0;
p0.ken = g{bid_exit+1}.attribute.nk;
p0.idir = 0;
p0.jdir = 1;
p0.kdir = 2;
p0.nxpid = pid1;
p0.nxbid = bid1;

g{bid_exit+1}.patch{pid_exit+1,1} = [];
g{bid_exit+1}.patch{pid_exit+1,1}.attribute = p0;

% Periodic patches
bid0 = bid_nozzle;
bid1 = bid_nozzle;

Q0 = [0 g{bid0+1}.attribute.ni ; 0 g{bid0+1}.attribute.nj ;...
    0 1];
dir0 = [0 1 2];
Q1 = [0 g{bid1+1}.attribute.ni ; 0 g{bid1+1}.attribute.nj ; g{bid1+1}.attribute.nk-1 g{bid1+1}.attribute.nk];
dir1 = [0 1 2];
g = ts_add_periodic_patch(bid0, bid1, Q0, Q1, dir0, dir1, g, 5);

for bid = ts_get_bids(g)
    g{bid+1}.attribute.np = length(ts_get_pids(g{bid+1}));
end

% Set initial guess
if isfield(g{bid_exit+1},'rovx') == 1
    g{bid_nozzle+1}.ro = ones(size(g{bid_nozzle+1}.x)) * g{bid_exit+1}.ro(end,round(nj/2),1);
    g{bid_nozzle+1}.rovx = ones(size(g{bid_nozzle+1}.x)) * g{bid_exit+1}.rovx(end,round(nj/2),1);
    g{bid_nozzle+1}.rovr = ones(size(g{bid_nozzle+1}.x)) * g{bid_exit+1}.rovr(end,round(nj/2),1);
    g{bid_nozzle+1}.rorvt = ones(size(g{bid_nozzle+1}.x)) * g{bid_exit+1}.rorvt(end,round(nj/2),1);
    g{bid_nozzle+1}.roe = ones(size(g{bid_nozzle+1}.x)) * g{bid_exit+1}.roe(end,round(nj/2),1);
    g{bid_nozzle+1}.trans_dyn_vis = ones(size(g{bid_nozzle+1}.x)) * g{bid_exit+1}.trans_dyn_vis(end,round(nj/2),1);
end

g{bid_nozzle+1}.bv.tstagin = g{bid_exit+1}.bv.tstagin;
g{bid_nozzle+1}.bv.tstagout = g{bid_exit+1}.bv.tstagout;
g{bid_nozzle+1}.bv.pstatin = g{bid_exit+1}.bv.pstatin;
g{bid_nozzle+1}.bv.pstatout = g{bid_exit+1}.bv.pstatout;
g{bid_nozzle+1}.bv.vgridin = g{bid_exit+1}.bv.vgridin;
g{bid_nozzle+1}.bv.vgridout = g{bid_exit+1}.bv.vgridout;

end


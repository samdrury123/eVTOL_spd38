function Q = ts_check_dir(p0,p1)
% Checks the 4 corners of two patches for coincidence

Q0a = [p0.x(1,1) ; p0.r(1,1) ; p0.rt(1,1)];
Q0b = [p0.x(1,end) ; p0.r(1,end) ; p0.rt(1,end)];
Q0c = [p0.x(end,end) ; p0.r(end,end) ; p0.rt(end,end)];
Q0d = [p0.x(end,1) ; p0.r(end,1) ; p0.rt(end,1)];
    
if p0.ni == 1
    Q1a = corner([-1 1 1], p0, p1);
    Q1b = corner([-1 1 2], p0, p1);
    Q1c = corner([-1 2 2], p0, p1);
    Q1d = corner([-1 2 1], p0, p1);
elseif p0.nj == 1
    Q1a = corner([1 -1 1], p0, p1);
    Q1b = corner([1 -1 2], p0, p1);
    Q1c = corner([2 -1 2], p0, p1);
    Q1d = corner([2 -1 1], p0, p1);
else
    Q1a = corner([1 1 -1], p0, p1);
    Q1b = corner([1 2 -1], p0, p1);
    Q1c = corner([2 2 -1], p0, p1);
    Q1d = corner([2 1 -1], p0, p1);        
end

Q = zeros(4,1);
if p0.nxbid ~= p0.bid && p0.kind ~= 2
    Q(1) = sum((Q0a.^2 - Q1a.^2).^0.5);
    Q(2) = sum((Q0b.^2 - Q1b.^2).^0.5);
    Q(3) = sum((Q0c.^2 - Q1c.^2).^0.5);
    Q(4) = sum((Q0d.^2 - Q1d.^2).^0.5);
end

end

function Q1 = corner(ijk0,p0,p1)

i0 = ijk0(1);
j0 = ijk0(2);
k0 = ijk0(3);

idir = p0.idir;
jdir = p0.jdir;
kdir = p0.kdir;

A = [1 2];

i1 = -1;
j1 = -1;
k1 = -1;

% Deal with idir
if idir == 0
    i1 = i0;
elseif idir == 1
    j1 = i0;
elseif idir == 2
    k1 = i0;
elseif idir == 3
    i1 = A(A~=i0);
elseif idir == 4
    j1 = A(A~=i0);
elseif idir == 5
    k1 = A(A~=i0);
else
    
end

% Deal with jdir
if jdir == 0
    i1 = j0;
elseif jdir == 1
    j1 = j0;
elseif jdir == 2
    k1 = j0;
elseif jdir == 3
    i1 = A(A~=j0);
elseif jdir == 4
    j1 = A(A~=j0);
elseif jdir == 5
    k1 = A(A~=j0);
else
    
end

% Deal with kdir
if kdir == 0
    i1 = k0;
elseif kdir == 1
    j1 = k0;
elseif kdir == 2
    k1 = k0;
elseif kdir == 3
    i1 = A(A~=k0);
elseif kdir == 4
    j1 = A(A~=k0);
elseif kdir == 5
    k1 = A(A~=k0);
else
    
end
    
% Remove singular row
ijk1 = [i1 j1 k1];
ijk1(ijk1==-1) = [];

if ijk1(1) == 2
    ijk1(1) = size(p1.x,1);
end

if ijk1(2) == 2
    ijk1(2) = size(p1.x,2);
end

Q1 = [p1.x(ijk1(1),ijk1(2)) ; p1.r(ijk1(1),ijk1(2)) ; p1.rt(ijk1(1),ijk1(2))];


end

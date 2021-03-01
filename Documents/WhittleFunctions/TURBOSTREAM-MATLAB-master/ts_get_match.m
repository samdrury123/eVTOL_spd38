function p1 = ts_get_match(p0,p1,b1)
% TS_GET_MATCH  Pull out aligned coordinates of a matching periodic patch

% Indices of matching patch
I = p1.ist+1:p1.ien; J = p1.jst+1:p1.jen; K = p1.kst+1:p1.ken;
p1.ni = p1.ien - p1.ist; p1.nj = p1.jen - p1.jst; p1.nk = p1.ken - p1.kst;

% Generate array of indices
ni = b1.attribute.ni; nj = b1.attribute.nj; nk = b1.attribute.nk; 
b1.I = repmat(reshape(1:ni,[ni 1 1]),[1 nj nk]);
b1.J = repmat(reshape(1:nj,[1 nj 1]),[ni 1 nk]);
b1.K = repmat(reshape(1:nk,[1 1 nk]),[ni nj 1]);
p1.I = squeeze(b1.I(I,J,K)); p1.J = squeeze(b1.J(I,J,K)); p1.K = squeeze(b1.K(I,J,K));

% Extract coordinates from matching block
p1.x = squeeze(b1.x(I,J,K)); p1.r = squeeze(b1.r(I,J,K)); p1.rt = squeeze(b1.rt(I,J,K));
p1.t = p1.rt./p1.r;

% Directions of matching coordinates to on current patch
ds = [p0.idir p0.jdir p0.kdir];

% Set singular dimension based on size of current patch
if p0.ni == 1; ds(1) = -1; end;
if p0.nj == 1; ds(2) = -1; end;
if p0.nk == 1; ds(3) = -1; end;

% Initialise swapping and flipping vectors
or = nan(1,3); fl = nan(1,3);

% Deal with directions
for n = 1:3
    if ds(n) == 0
        or(n) = 1; fl(1) = 0;
    elseif ds(n) == 1
        or(n) = 2; fl(2) = 0;
    elseif ds(n) == 2
        or(n) = 3; fl(3) = 0;
    elseif ds(n) == 3
        or(n) = 1; fl(1) = 1;
    elseif ds(n) == 4
        or(n) = 2; fl(2) = 1;
    elseif ds(n) == 5
        or(n) = 3; fl(3) = 1;
    end
end

% Remove singular dimension
or = or(isnan(or) == 0); fl = fl(isnan(fl) == 0);
or = (or - min(or)) ./ (max(or) - min(or)) + 1;

% Flip the coordinate arrays
varnames = {'x' 'r' 'rt' 't' 'I' 'J' 'K'};
for n = 1:2
    if fl(n) == 1
        for v = 1:length(varnames)
            p1.(varnames{v}) = flip(p1.(varnames{v}),n);
        end
    end
end

% Permute the coordinate arrays
for v = 1:length(varnames)
    p1.(varnames{v}) = permute(p1.(varnames{v}),or);
end


end


function [f,v] = triangulate_grid(xyz,f,v)
% Triangulate a 2D structured mesh and append if f and v are specified already

% Empty f and v array if not specified
if exist('f','var') == 0 || isempty(f) == 1
    f = []; f_max = 0;
else
    f_max = max(f(:));
end
if exist('v','var') == 0
    v = [];
end

% Size of array
[ni,nj,~] = size(xyz);

% Reshape all vertices into a list
v = [v ; reshape(xyz,[ni * nj 3])];

% Index number in array format
n = reshape(1:ni*nj,[ni nj]);

% Vertex numbers for bottom right triangles
f_1 = [reshape(n(1:end-1,1:end-1),[],1) reshape(n(2:end,1:end-1),[],1) reshape(n(2:end,2:end),[],1)];

% Vertex numbers for top left triangles
f_2 = [reshape(n(1:end-1,1:end-1),[],1) reshape(n(2:end,2:end),[],1) reshape(n(1:end-1,2:end),[],1)];

% Return list of all vertex numbers
f = [f ; f_1 + f_max ; f_2 + f_max];


end
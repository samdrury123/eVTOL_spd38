function [vnew, fnew]= stlSlimVerts(v, f, tol)
% PATCHSLIM removes duplicate vertices in surface meshes.
% 
% This function finds and removes duplicate vertices.
%
% USAGE: [v, f]=patchslim(v, f)
%
% Where v is the vertex list and f is the face list specifying vertex
% connectivity.
%
% v contains the vertices for all triangles [3*n x 3].
% f contains the vertex lists defining each triangle face [n x 3].
%
% This will reduce the size of typical v matrix by about a factor of 6.
%
% For more information see:
%  http://www.esmonde-white.com/home/diversions/matlab-program-for-loading-stl-files
%
% Francis Esmonde-White, May 2010

% Default to zero tolerance
if exist('tol','var') == 0
    tol = 0;
end

if ~exist('v','var')
    error('The vertex list (v) must be specified.');
end
if ~exist('f','var')
    error('The vertex connectivity of the triangle faces (f) must be specified.');
end

% Round vertices to nearest tolerance value if specfied
if tol > 0
    v = round(v / tol) * tol;
end

[vnew, indexm, indexn] =  unique(v, 'rows');
fnew = indexn(f);

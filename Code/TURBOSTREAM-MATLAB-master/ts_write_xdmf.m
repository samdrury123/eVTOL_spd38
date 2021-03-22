function [] = ts_write_xdmf(g,filename,bp,bids)
% TS_WRITE_XDMF  Write a file to interpret TURBOSTREAM data in PARAVIEW
%
%   [] = TS_WRITE_XDMF(g,filename,bp,bids)
%
%   g - cell array of TURBOSTREAM solution
%   filename - string of .xdmf filename to write to
%   bp - limitations on block properties to write in cell array format
%   bids - vector of block ids to include or leave empty for all

if exist('bp','var') == 0
    attributes = {'ro' 'Float' 4; 'rovx' 'Float' 4; 'rovr' 'Float' 4; 'rorvt' 'Float' 4;...
        'roe' 'Float' 4; 'mwall' 'Char' 1; 'tlength' 'Float' 4; 'trans_dyn_vis' 'Float' 4;...
        'xlength' 'Float' 4;  'misc' 'Float' 4; 'dwallsq' 'Float' 4; 'phi' 'Float' 4; ...
        'yplus' 'Float' 4; 'egen' 'Float' 4; 'tdamp' 'Float' 4; 'fac_st0' 'Float' 4};
else
    attributes = bp;
end

if exist('bids','var') == 0
    bids = ts_get_bids(g);
end

% Check whether mesh is in cartesian or polar coordinates
if isfield(g{1},'rt') == 1
    geometry = { 'x' 'Float' ; 'r' 'Float' ; 'rt' 'Float'};
else
    geometry = { 'x' 'Float' ; 'y' 'Float' ; 'z' 'Float'};
end

% Remove velocities from output
velocities = { 'Vx' 'Float' 4 ; 'Vr' 'Float' 4 ; 'Vt_rel' 'Float' 4};
for a = size(attributes,1):-1:1
    for v = 1:size(velocities,1)
        if strcmp(attributes{a,1},velocities{v,1}) == 1
            attributes(a,:) = [];
        end
    end
end

% Convert filename to xdmf if not already
if isempty(regexp(filename,'.xdmf', 'once')) == 1
    xdmf_filename = strrep(filename,'.hdf5','.xdmf');
else
    xdmf_filename = filename;
    filename = strrep(filename,'.xdmf','.hdf5');
end

% Break the filename into parts
[~,filename,~] = fileparts(filename);
filename = [filename '.hdf5'];

% Open the file for printing
fid = fopen(xdmf_filename,'w');

% Print headers
fprintf(fid,'<?xml version=''1.0'' ?>\n');
fprintf(fid,'<!DOCTYPE Xdmf SYSTEM ''Xdmf.dtd'' []>\n');
fprintf(fid,'<Xdmf>\n');
fprintf(fid,'\t<Domain>\n');
fprintf(fid,'\t\t<Grid Name=''grid'' GridType=''Collection'' CollectionType=''Spatial''>\n');

% Print block info
for bid = bids
    fprintf(fid,'\t\t\t<Grid Name=''block%d'' GridType=''Uniform''>\n',bid);
    ni = g{bid+1}.attribute.ni;
    nj = g{bid+1}.attribute.nj;
    nk = g{bid+1}.attribute.nk;
    fprintf(fid,'\t\t\t\t<Topology Type=''3DSMESH'' Dimensions=''%d %d %d''/>\n', nk, nj, ni);
    
    % Print mesh
    fprintf(fid,'\t\t\t\t<Geometry Type=''X_Y_Z''>\n');
    for n = 1:size(geometry,1)
        coord = geometry{n,1};
        fprintf(fid,['\t\t\t\t\t<DataItem DataType=''Float'' Precision=''4'' Dimensions=''%d %d %d'' '...
            'Format=''HDF''>\n'], nk, nj, ni);
        fprintf(fid,'\t\t\t\t\t\t%s:/block%d/%s_bp\n',filename,bid,coord);
        fprintf(fid,'\t\t\t\t\t</DataItem>\n');
    end
    fprintf(fid,'\t\t\t\t</Geometry>\n');
    
    % Print block properties
    varnames = fieldnames(g{bid+1});
    for v = 1:length(varnames)
        varname = varnames{v};
        for a = 1:size(attributes,1)
            attribute = attributes{a,1};
            if strcmp(attribute,varname) == 1
                fprintf(fid,['\t\t\t\t<Attribute Name=''%s'' AttributeType=''Scalar'' '...
                    'Center=''Node''>\n'],attribute);
                fprintf(fid,['\t\t\t\t\t<DataItem Dimensions=''%d %d %d'' '...
                    'NumberType=''%s'' Precision=''%d'' Format=''HDF''>\n'],nk, nj, ni,...
                    attributes{a,2}, attributes{a,3});
                fprintf(fid,'\t\t\t\t\t\t%s:/block%d/%s_bp\n',filename,bid,attribute);
                fprintf(fid,'\t\t\t\t\t</DataItem>\n');
                fprintf(fid,'\t\t\t\t</Attribute>\n');
            end
        end
    end
    
    % Print velocity tensor
    if isfield(g{1},'Vx') == 1
        fprintf(fid,['\t\t\t\t<Attribute Name=''%s'' AttributeType=''Vector'' '...
            'Center=''Node''>\n'],'V');
        fprintf(fid,['\t\t\t\t\t<DataItem Dimensions=''%d %d %d %d'' '...
            'NumberType=''%s'' Precision=''%d'' Format=''HDF''>\n'], nk, nj, ni, 3,...
            'Float', 4);
        fprintf(fid,'\t\t\t\t\t\t%s:/block%d/%s_bp\n',filename,bid,'V');
        fprintf(fid,'\t\t\t\t\t</DataItem>\n');
        fprintf(fid,'\t\t\t\t</Attribute>\n');
    end
    
    fprintf(fid,'\t\t\t</Grid>\n');
end

fprintf(fid,'\t\t</Grid>\n');
fprintf(fid,'\t</Domain>\n');
fprintf(fid,'</Xdmf>\n');
fprintf(fid,'\n');

fclose(fid);


end

    
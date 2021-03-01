function [] = ts_write_hdf5(g,filename,bp,bids)
% TS_WRITE_HDF5  Write a heirarchical file to store TURBOSTREAM data
%
%   [] = TS_WRITE_HDF5(g,filename,bp,bids)
%
%   g - cell array of TURBOSTREAM solution
%   filename - string of .hdf5 filename to write to
%   bp - limitations on block properties to write in cell array format
%   bids - vector of block ids to include or leave empty for all

% Delete file if it already exists
if exist(filename,'file') ~= 0
    delete(filename);
end

% Define required block properties if not specified
if exist('bp','var') == 0
    bp = { 'x' 'fl' ; 'r' 'fl' ; 'rt' 'fl' ; 'ro' 'fl' ; 'rovx' 'fl' ; 'rovr' 'fl' ; 'rorvt' 'fl' ;...
        'roe' 'fl'; 'mwall' 'in8'; 'tlength' 'fl'; 'trans_dyn_vis' 'fl'; 'xlength' 'fl';...
        'misc' 'fl'; 'dwallsq' 'fl' ; 'phi' 'fl' ; 'yplus' 'fl' ; 'egen' 'fl' ; 'tdamp' 'fl' ; ...
        'fac_st0' 'fl'};
    ts_type; pv_mode = 0;
else
    pv_mode = 1;
end

% Get all block ids if not specifiec
if exist('bids','var') == 0
    bids = ts_get_bids(g);
end

% Reset some problems with patches before writing the file
if pv_mode == 0
    for bid = ts_get_bids(g)
        for pid = ts_get_pids(g,bid)
            
            % Reset probe append flag
            if g{bid+1}.patch{pid+1}.attribute.kind == 8
                g{bid+1}.patch{pid+1}.probe_append = 0;
            end
            
            % Delete sliding information for non-sliding patches
            if g{bid+1}.patch{pid+1}.attribute.kind ~= 3
                if isfield(g{bid+1}.patch{pid+1},'slide_nxbid') == 1
                    g{bid+1}.patch{pid+1} = rmfield(g{bid+1}.patch{pid+1},'slide_nxbid');
                end
                if isfield(g{bid+1}.patch{pid+1},'slide_nxpid') == 1
                    g{bid+1}.patch{pid+1} = rmfield(g{bid+1}.patch{pid+1},'slide_nxpid');
                end
            end
            
        end
    end
end


%% Write the file

% Open the file
fid = H5F.create(filename,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');

% Create group creation property lists
gcpl = H5P.create('H5P_LINK_CREATE');
H5P.set_create_intermediate_group(gcpl, 1);

% Close and reopen the file
H5F.close(fid); fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');

% Write application attributes
write_atts(fid,'/',{'nb','ntb'},{length(bids),0});

% Write application variables
if pv_mode == 0
    [varnames,varvalues,varclass] = get_vars(g{1}.av,av,'_av');
    write_vars(fid,'/',varnames,varvalues,1,varclass)
end

% Loop over all blocks
for bid = bids

    % Block group
    bg = ['/block' num2str(bid)];
    gid = H5G.create(fid,bg,gcpl,'H5P_DEFAULT','H5P_DEFAULT');
    H5G.close(gid);
    
    % Block size
    ni = double(g{bid+1}.attribute.ni); nj = double(g{bid+1}.attribute.nj); 
    nk = double(g{bid+1}.attribute.nk);
    
    % Write block properties
    [varnames,varvalues,varclass] = get_vars(g{bid+1},bp,'_bp');
    write_vars(fid,bg,varnames,varvalues,[ni nj nk],varclass)
    
    % Write block velocity tensor
    if pv_mode == 1 && isfield(g{bid+1},'Vx') == 1
        
        % Assemble matrix of velocities
        V = single(zeros(3,ni*nj*nk)); V(1,:) = g{bid+1}.Vx(:); V(2,:) = g{bid+1}.Vr(:); 
        if g{bid+1}.bv.rpm == 0; V(3,:) = g{bid+1}.Vt(:); else; V(3,:) = g{bid+1}.Vt_rel(:); end; 
        
        % Write the velocites
        write_vars(fid,bg,{'V_bp'},{V},[3 ni nj nk],{'fl'})
    end
    
    % Write variables only required for running in TURBOSTREAM
    if pv_mode == 0
    
        % Write block variables
        [varnames,varvalues,varclass] = get_vars(g{bid+1}.bv,bv,'_bv');
        write_vars(fid,bg,varnames,varvalues,1,varclass)        

        % Write block attributes
        [attnames,attvalues] = get_vars(g{bid+1}.attribute,ba,'');
        write_atts(fid,bg,attnames,attvalues)

        % Loop over all patches
        for pid = ts_get_pids(g,bid)

            % Patch group
            pg = [bg '/patch' num2str(pid)];
            gid = H5G.create(fid,pg,gcpl,'H5P_DEFAULT','H5P_DEFAULT');
            H5G.close(gid);
            
            % Write patch properties
            [varnames,varvalues,varclass] = get_vars(g{bid+1}.patch{pid+1},pp,'_pp');
            if isempty(varvalues) == 0
                write_vars(fid,pg,varnames,varvalues,size(varvalues{1}),varclass);
            end
            
            % Write patch variables
            [varnames,varvalues,varclass] = get_vars(g{bid+1}.patch{pid+1},pv,'_pv');
            if isempty(varvalues) == 0
                varnames = strrep(varnames,'id_pv','id');
                write_vars(fid,pg,varnames,varvalues,size(varvalues{1}),varclass);
            end
            
            % Write patch attributes
            [attnames,attvalues] = get_vars(g{bid+1}.patch{pid+1}.attribute,pa,'');
            write_atts(fid,pg,attnames,attvalues)
            
        end
    end
end    

% Close the file
H5P.close(gcpl); H5F.close(fid);


end


function [varnames,varvalues,varclass] = get_vars(s,v,suffix)
% Get variable names, values and class from a struct

% Check which variables are actually present in the struct
v_act = fieldnames(s);
q = strcmp(repmat(v_act,[1 size(v,1)]),repmat(v(:,1)',[length(v_act) 1]));
v = v(any(q,1),:); 

% Get all values, names and classes
varvalues = cell(size(v,1),1); varnames = cell(size(v,1),1); varclass = cell(size(v,1),1);
for n = 1:size(v,1)
    varnames{n} = [v{n,1} suffix]; 
    varvalues{n} = s.(v{n,1});
    varclass{n} = v{n,2};
end

% Convert to specified class
for n = 1:length(varvalues)
    if strcmp(varclass{n},'fl') == 1
        varvalues{n} = single(varvalues{n}); 
    elseif strcmp(varclass{n},'in') == 1
        varvalues{n} = int32(varvalues{n}); 
    elseif strcmp(varclass{n},'in8') == 1
        varvalues{n} = uint8(varvalues{n}); 
    end    
end


end


function write_vars(fid,oname,varnames,varvalues,varsize,varclass)
% Write a list of variables into the new file

% Class types to look for
cl_types = {'in' 'H5T_NATIVE_INT' ; 'fl' 'H5T_NATIVE_FLOAT' ; 'in8' 'H5T_NATIVE_UCHAR'};

% Open the object
oid = H5O.open(fid,oname,'H5P_DEFAULT');

% Create space identifier
dims = fliplr(varsize);
sid = H5S.create_simple(numel(varsize),dims,dims);

% Write all classes
for m = 1:size(cl_types,1)

    % Find variables with the current class
    i = find(strcmp(varclass,cl_types{m,1}) == 1);
    
    % Create class identifier
    tid = H5T.copy(cl_types{m,2});

    % Write the variables
    for n = 1:length(i)
        vid = H5D.create(oid,varnames{i(n)},tid,sid,'H5P_DEFAULT');
        H5D.write(vid,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',varvalues{i(n)}) 
        H5D.close(vid);
    end

    % Close class identifier
    H5T.close(tid);
    
end

% Close all identifiers
H5O.close(oid); H5S.close(sid);


end


function write_atts(fid,oname,attnames,attvalues)
% Write a list of attributes into the new file

% Open the object
oid = H5O.open(fid,oname,'H5P_DEFAULT');

% Create attribute identifiers
tid = H5T.copy('H5T_NATIVE_INT');
sid = H5S.create_simple(1,1,1);

% Write the attributes
for n = 1:length(attnames)
    aid = H5A.create(oid,attnames{n},tid,sid,'H5P_DEFAULT');
    H5A.write(aid,tid,int32(attvalues{n}));
    H5A.close(aid);
end

% Close all identifiers
H5O.close(oid); H5T.close(tid); H5S.close(sid);


end


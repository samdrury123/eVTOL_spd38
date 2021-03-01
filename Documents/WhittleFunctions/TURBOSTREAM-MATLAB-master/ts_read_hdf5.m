function g = ts_read_h5(filename,read_patches,in_parallel)
% TS_READ_H5 Reads in data from a TURBOSTREAM HDF5 file
% Author Alejandro Castillo feb 2020, based on ts_read_hdf5 from James Taylor


% Default to read patches
if exist('read_patches','var') == 0
    read_patches = 1;
end

% Default to not read in parallel
if exist('in_parallel','var') == 0
    in_parallel = 0;
end

% Make sure newest modified file is always read
A_hdf5 = dir(filename); A_mat = dir(strrep(filename,'.hdf5','.mat'));
if numel(A_hdf5) == 1 && numel(A_mat) == 1 
    if A_hdf5.datenum > A_mat.datenum
        delete(strrep(filename,'.hdf5','.mat'))
    end
end

% Necessary block names for regular perfect gas computions
bp_names = {'x' 'r' 'rt' 'mwall' 'phi' 'ro' 'rovx' 'rovr' 'rorvt' 'roe' 'tlength'...
    'trans_dyn_vis' 'xlength' 'yplus' 'tdamp' 'fac_st0'};

% Set number of parallel workers
if in_parallel == 0
    n_parfor = 0;
else
    n_parfor = 8;
    p=gcp('nocreate');
    if isempty(p)
        parpool(n_parfor);
    else
        n_parfor = p.NumWorkers;
    end
end

% Check if matlab equivalent exists, if not read hdf5 file
if exist(strrep(filename,'.hdf5','.mat'),'file') == 0
    
    % Read structure of HDF5 file
    info = h5info(filename);
    nGroups = size(info.Groups,1);
    nBlocks = nGroups;
    nGlobalDatasets = size(info.Datasets,1);

    % Read in application variables
    for i = 1:nGlobalDatasets
        varName=info.Datasets(i,1).Name;
        varNameFull=['/' varName];
        varName=varName(1:end-3);
        appVars.(varName) = h5read(filename,varNameFull);
    end

    % Read in the variables and coordinates from every block
    g = cell(nBlocks,1); iblock = zeros(nBlocks,1);
%     parfor (i = 1:nBlocks,n_parfor)
    for i = 1:nBlocks

        % Extract block number, not ordered correctly in HDF5 file
        blockStr = info.Groups(i,1).Name;
        nBlock = str2double(blockStr(7:end))+1;
        iblock(i) = nBlock;
        
        % Read in all block properties and variables
        nGroupDatasets = size(info.Groups(i,1).Datasets,1);
        for j = 1:nGroupDatasets
            varNameFull=info.Groups(i,1).Datasets(j,1).Name;
            varName=varNameFull(1:end-3);
            g_temp.av = appVars;
            
            
            % Check if the flow is a necessary block property or block variable
            if strcmp(varNameFull(end-2:end),'_bv') == 1 || any(strcmp(varName,bp_names)) == 1
                g_temp.(varName) = h5read(filename,['/block' int2str(nBlock-1) '/' varNameFull]);

                % Move the block variables into a separate sub-structure
                if size(g_temp.(varName),1) == 1 && ...
                        size(g_temp.(varName),2) == 1 && ...
                        size(g_temp.(varName),3) == 1
                    g_temp.bv.(varName) = g_temp.(varName);
                    g_temp = rmfield(g_temp, varName);

                else
                    % Reshape the block properties so they follow the I, J, K convention
                    [K,J,I] = size(g_temp.(varName));
                    g_temp.(varName) = reshape(g_temp.(varName), I, J, K);
                end
            end
        end

        % Read block attributes
        nGroupAttributes = size(info.Groups(i,1).Attributes,1);
        for j = 1:nGroupAttributes
            attName = info.Groups(i,1).Attributes(j,1).Name;
            g_temp.attribute.(attName) = info.Groups(i,1).Attributes(j,1).Value;
        end
        
        % Read patches
        if exist('read_patches','var') && read_patches
            np = g_temp.attribute.np;
            ncl = g_temp.attribute.ncl;
            g_temp.patch = cell(np,1);
            for j = ncl+1:np+ncl
                groupName = info.Groups(i,1).Groups(j,1).Name;
                pName = info.Groups(i,1).Groups(j,1).Name;
                pid = h5readatt(filename,groupName,'pid');
                for k = 1:size(info.Groups(i,1).Groups(j,1).Attributes,1)
                    attName = info.Groups(i,1).Groups(j,1).Attributes(k,1).Name;
                    g_temp.patch{pid+1}.attribute.(attName) = info.Groups(i,1).Groups(j,1).Attributes(k,1).Value;
                end

                for k = 1:size(info.Groups(i,1).Groups(j,1).Datasets,1)
                    varName = info.Groups(i,1).Groups(j,1).Datasets(k,1).Name;
                    varNameFull = [pName '/' varName];
                    if ~strcmp(varName,'slide_nxbid') && ~strcmp(varName,'slide_nxpid')
                        varName = varName(1:end-3);
                    end
                    g_temp.patch{pid+1}.(varName) = h5read(filename, varNameFull);
                end
            end

            % Read clists
            for j = 1:ncl
                groupName = info.Groups(i,1).Groups(j,1).Name;
                clid = h5readatt(filename,groupName,'clid');
                for k = 1:size(info.Groups(i,1).Groups(j,1).Attributes,1)
                    attName = info.Groups(i,1).Groups(j,1).Attributes(k,1).Name;
                    g_temp.clist{clid+1}.attribute.(attName) = info.Groups(i,1).Groups(j,1).Attributes(k,1).Value;
                end

                for k = 1:size(info.Groups(i,1).Groups(j,1).Datasets,1)
                    varName = info.Groups(i,1).Groups(j,1).Datasets(k,1).Name;
                    varNameFull = [groupName '/' varName];
                    g_temp.clist{clid+1}.(varName) = h5read(filename, varNameFull);
                end
            end

        end
        
        % Record current block
        g{i} = g_temp;

    end
    
    % Rearrange block order
    [~,i] = sort(iblock);
    g = g(i);

    % Save matlab file and delete hdf5
    if read_patches == 1
        save(strrep(filename,'.hdf5','.mat'),'g','-v7.3')
        delete(filename);
    end
    
else
    
    % Load data direct from matlab file
    g = load(strrep(filename,'.hdf5','.mat')); g = g.g;
    
end


end


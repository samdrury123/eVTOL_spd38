function G = ts_read_probes(directory,filestring)
% Reads in a set of TURBOSTREAM probe files for all timesteps

% Find all probe files
F = dir([directory filestring '_probe_*.*']);
G = cell(length(F),1);

% Get nsteps and sort
nstep = zeros(length(F),1);
for f = 1:length(F)
    A = regexp(F(f).name,'_','split');
    nstep(f) = str2double(A{end}(1:end-5));
end
[nstep,i] = sort(nstep);
F = F(i);

% Necessary block names for regular perfect gas computions
bp_names = {'x' 'r' 'rt' 'mwall' 'phi' 'ro' 'rovx' 'rovr' 'rorvt' 'roe' 'tlength'...
    'trans_dyn_vis' 'xlength' 'yplus' 'tdamp'};

% Loop over all probes and read data
for f = 1:length(F)    
    filename = [directory F(f).name];
    if exist(strrep(filename,'.hdf5','.mat'),'file') == 0
        
        % Open file and check size
        info = hdf5info(filename);
        nGroups = size(info.GroupHierarchy.Groups,2);
        nBlocks = nGroups;

        % Open the file
        fid = H5F.open(filename);

        % Loop over each entry and record in cell array
        g = cell(nBlocks,1);
        for i = 1:nBlocks    
            
            % Get block and patch data on original probe patch
            blockStr = info.GroupHierarchy.Groups(1,i).Name;
            A = regexp(blockStr,'_','split');
            bid = str2double(A{1}(7:end));
            pid = str2double(A{2}(6:end));

            % Loop over each variable and store into cell array
            nGroupDatasets = size(info.GroupHierarchy.Groups(1,i).Datasets,2);
            for j = 1:nGroupDatasets
                
                % Read the data
                varNameFull = info.GroupHierarchy.Groups(1,i).Datasets(1,j).Name;
                dset_id = H5D.open(fid,varNameFull);
                A = regexp(varNameFull,'/','split');
                varName = strrep(A{end},'_bp','');
                data = H5D.read(dset_id);

                % Determine whether data is block properties, variables or attributes
                if numel(data) == 1
                    if strcmp(varName(1),'n') == 0
                        g{i}.bv.(varName) = data;
                    else
                        g{i}.attribute.(varName) = data;
                    end
                else
                    if any(strcmp(varName,bp_names)) == 1
                        [K, J, I] = size(data);
                        g{i}.(varName) = squeeze(reshape(data,[I J K]));
                    end
                end
            end

            % Record data location and timestep
            g{i}.attribute.bid = bid;
            g{i}.attribute.pid = pid;
            g{i}.attribute.nstep = nstep(f);

        end
        
        % Save mat file and store in cell array
        save(strrep(filename,'.hdf5','.mat'),'g');
        delete(filename);
        G{f} = g;

    else
        
        % Load directly from mat file
        load(strrep(filename,'.hdf5','.mat'));
        G{f} = g;
        
    end
end


end


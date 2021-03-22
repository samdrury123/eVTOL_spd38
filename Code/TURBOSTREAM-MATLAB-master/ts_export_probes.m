function [] = ts_export_probes(directory,file_prefix,bp,frequency,nstep_cycle,...
    conv_cart,p,nsteps,t_pitch,skip_mat,domain)
% Write all probes at all timesteps into a paraview compatible hdf5 file
% 23Jan17: edited by hv256 so patch xdmf files are compatible with xdmf reader from paraview to read in 2D probe patches 

if exist('skip_mat','var') == 0
    skip_mat = 0;
end

if exist('domain','var') == 0
    domain = 0;
end

if exist('p','var') == 0
    p = [];
end

if strcmp(bp,'HighSpeed') == 1
     bp = {'M_rel', 'sfunc', 'Alpha_rel'};
end
if strcmp(bp,'LowSpeed') == 1
     bp = {'P', 'Po_rel', 'Alpha_rel'};
end
if strcmp(bp,'Track') == 1
    bp = {'ro','roe'};
end

if exist('nsteps','var') == 0
    % Get nsteps from available hdf5 files
    F = dir([directory file_prefix '_probe_*.hdf5']);

    % Get nsteps and sort
    nsteps = zeros(length(F),1);
    for f = 1:length(F)
        A = regexp(F(f).name,'_','split');
        nsteps(f) = str2double(A{end}(1:end-5));
    end
    [nsteps,i] = sort(nsteps);
    F = F(i);
else
    % Create filenames based on specified nsteps
%     F = cell(length(nsteps),1);
    for f = 1:length(nsteps)
        F(f).name = [file_prefix '_probe_' num2str(nsteps(f)) '.hdf5'];
    end
end

% Calculate ellapsed time for each nstep
times = nsteps / (nstep_cycle  * frequency);

% Save all probes as mat files for speed of reading later
if skip_mat == 0
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
                varNameFull = info.GroupHierarchy.Groups(1,i).Datasets(1,j).Name;

                dset_id = H5D.open(fid,varNameFull);

                A = regexp(varNameFull,'/','split');
                varName = strrep(A{end},'_bp','');

    %             data = hdf5read(filename,varNameFull);
                data = H5D.read(dset_id);

                % Determine whether data is block properties, variables or
                % attributes
                if numel(data) == 1
                    if strcmp(varName(1),'n') == 0
                        g{i}.bv.(varName) = data;
                    else
                        g{i}.attribute.(varName) = data;
                    end
                else
                    [K, J, I] = size(data);
                    g{i}.(varName) = squeeze(reshape(data,[I J K]));
                end
            end

            % Record data location and timestep
            g{i}.attribute.bid = bid;
            g{i}.attribute.pid = pid;
            g{i}.attribute.nstep = nsteps(f);
            g{i}.attribute.time = times(f);

        end
        
        % Save as a mat file
        save(strrep(filename,'.hdf5','.mat'),'g');
        delete(filename)
        disp(F(f).name)
    end
end
end

% Determine cell array for hdf5 file
bp_hdf5_master = {'ro' 'fl' ; 'rovx' 'fl' ; 'rovr' 'fl' ; 'rorvt' 'fl' ; 'roe' 'fl';...
    'mwall' 'in8'; 'tlength' 'fl'; 'trans_dyn_vis' 'fl'; 'xlength' 'fl';  'misc' 'fl';...
    'dwallsq' 'fl' ; 'phi' 'fl' ; 'yplus' 'fl' ; 'egen' 'fl' ; 'tdamp' 'fl';...
    'n_x' 'fl' ; 'n_r' 'fl' ; 'n_rt' 'fl'};

bp_hdf5 = cell(length(bp),2);
for n = 1:length(bp)
    if sum(strcmp(bp(n),bp_hdf5_master(:,1))) > 0
        bp_hdf5{n,1} = bp{n};
        bp_hdf5{n,2} = bp_hdf5_master{strcmp(bp(n),bp_hdf5_master(:,1)),2};
    else
        bp_hdf5{n,1} = bp{n};
        bp_hdf5{n,2} = 'fl';
    end
end

% Determine cell array for xdmf file
bp_xdmf_master = {'ro' 'Float' 4; 'rovx' 'Float' 4; 'rovr' 'Float' 4; 'rorvt' 'Float' 4;...
    'roe' 'Float' 4; 'mwall' 'Char' 1; 'tlength' 'Float' 4; 'trans_dyn_vis' 'Float' 4;...
    'xlength' 'Float' 4;  'misc' 'Float' 4; 'dwallsq' 'Float' 4;...
    'phi' 'Float' 4; 'yplus' 'Float' 4; 'egen' 'Float' 4; 'tdamp' 'Float' 4;...
    'n_x' 'Float' 4; 'n_r' 'Float' 4; 'n_rt' 'Float' 4;};

bp_xdmf = cell(length(bp),3);
for n = 1:length(bp)
    if sum(strcmp(bp(n),bp_xdmf_master(:,1))) > 0
        bp_xdmf{n,1} = bp{n};
        bp_xdmf{n,2} = bp_xdmf_master{strcmp(bp(n),bp_xdmf_master(:,1)),2};
        bp_xdmf{n,3} = bp_xdmf_master{strcmp(bp(n),bp_xdmf_master(:,1)),3};
    else
        bp_xdmf{n,1} = bp{n};
        bp_xdmf{n,2} = 'Float';
        bp_xdmf{n,3} = 4;
    end
end

% Cut up probes into domains as specified by index if present
if isempty(p) == 1
    for f = 1:length(nsteps)
        load([directory file_prefix '_probe_' num2str(nsteps(f)) '.mat'])
        g = ts_secondary(g);
        filename = [directory file_prefix '_probe_' num2str(nsteps(f)) '.hdf5'];
        delete(filename);
        write_hdf5(g,filename,bp_hdf5,conv_cart,domain);
        write_xdmf(g,filename,bp_xdmf,times(f),nsteps(f),conv_cart);
    end
else
    load([directory file_prefix '_probe_' num2str(nsteps(1)) '.mat'])
    disp([file_prefix '_probe_' num2str(nsteps(1)) '.mat'])
    % Get probe ids for each probe set from p
    varnames = fieldnames(p);
    for v = 1:length(varnames)
        varname = varnames{v};

        bids = p.(varname)(:,1); pids = p.(varname)(:,2);
        prid_set.(varname) = [];
        for n = 1:length(bids)
            for prid = 0:length(g)-1
                if g{prid+1}.attribute.bid == bids(n) && g{prid+1}.attribute.pid == pids(n)
                    prid_set.(varname) = [prid_set.(varname) prid];
                end
            end
        end
    end
    
    % Remove empty set names
    for v = length(varnames):-1:1
        if isempty(prid_set.(varnames{v})) == 1
            varnames(v) = [];
        end
    end

    % Determine zero theta position for rotating domain in each probe set
    t_0 = zeros(length(varnames),1);
    for v = 1:length(varnames)
        % Get probe set
        varname = varnames{v};
        g_temp = ts_secondary(g(prid_set.(varname)+1));
        
        % Get rotating surfaces
        prids_rotating = [];
        for prid = 0:size(p.(varname),1)-1
            if p.(varname)(prid+1,3) ~= 0
                prids_rotating = [prids_rotating prid];
            end
        end

        % Determine mid theta coordinate
        if isempty(prids_rotating) == 0
            t_min = inf; t_max = -inf;
            for prid = prids_rotating
                t_min = min(min(min(g_temp{prid+1}.t)),t_min);
                t_max = max(max(max(g_temp{prid+1}.t)),t_max);
            end
            t_0(v) = 0.5 * (t_min + t_max);
        end
    end
    
    for f = 1:length(nsteps)
        
        % Check what files exist already
        file_flag = zeros(length(varnames),2);
        for v = 1:length(varnames)
            varname = varnames{v};
            if exist([directory file_prefix '_probe_' num2str(nsteps(f)) '_'...
                    varname '.hdf5'],'file') ~= 0
                file_flag(v,1) = 1;
            end
            if exist([directory file_prefix '_probe_' num2str(nsteps(f)) '_'...
                    varname '.xdmf'],'file') ~= 0
                file_flag(v,2) = 1;
            end
            
        end
        
        % Read in probe at current step
        if sum(sum(file_flag == 0)) > 0
            load([directory file_prefix '_probe_' num2str(nsteps(f)) '.mat'])
            disp([file_prefix '_probe_' num2str(nsteps(f)) '.mat'])
            g = ts_secondary(g);

            % Loop over all sets, extract relevant probes then write hdf5 and
            % xdmf files
            for v = 1:length(varnames)
                varname = varnames{v};
                filename = [directory file_prefix '_probe_' num2str(nsteps(f)) '_' varname '.hdf5'];
                if exist(filename,'file') == 0 || exist(strrep(filename,'.hdf5','.xdmf'),'file') == 0
    %                 delete(filename);
                    g_temp = g(prid_set.(varname)+1);

                    % Analyse in stationary or rotor relative frame
                    if domain == 0
                        g_temp = shift_rotors(p.(varname),g_temp,t_pitch);
                    else
                        g_temp = freeze_rotors(p.(varname),g_temp,t_pitch,t_0(v));
                    end

                    % Determine whether to use xyz or xrrt coordinates
                    if strcmp(varname(1),'R') == 1 && strcmp(varname(2),'_') == 0
                        bp_hdf5_temp = [{ 'x' 'fl' ; 'r' 'fl' ; 'rt' 'fl'} ; bp_hdf5];
                        g_temp = flatten_probes(g_temp);
                        write_hdf5(g_temp,filename,bp_hdf5_temp,0,domain);
                        write_xdmf(g_temp,filename,bp_xdmf,times(f),nsteps(f),0);
                    else
                        bp_hdf5_temp = [{ 'x' 'fl' ; 'y' 'fl' ; 'z' 'fl'} ; bp_hdf5];
                        for n = 1:length(g_temp)
                            [g_temp{n}.y, g_temp{n}.z] = pol2cart(g_temp{n}.rt./g_temp{n}.r,...
                                g_temp{n}.r);
                        end
                        write_hdf5(g_temp,filename,bp_hdf5_temp,1,domain);
                        write_xdmf(g_temp,filename,bp_xdmf,times(f),nsteps(f),1);
                    end

                end
            end
        end
    end
    % Write xdmf file referring to all timesteps
    for v = 1:length(varnames)
        varname = varnames{v};
        fid = fopen([directory file_prefix '_probes_' num2str(min(nsteps)) '_' ...
            num2str(max(nsteps)) '_' varname '.xdmf'],'w');
        fprintf(fid,'<?xml version=''1.0'' ?>\n');
        fprintf(fid,'<Xdmf xmlns:xi=''http://www.w3.org/2001/XInclude'' Version=''2.0''>\n');
        fprintf(fid,'\t<Domain>\n');
        fprintf(fid,'\t\t<Grid GridType=''Collection'' CollectionType=''Temporal''>\n');
        for f = 1:length(nsteps)
            filename = [file_prefix '_probe_' num2str(nsteps(f)) '_' varname '.xdmf'];
            fprintf(fid,'\t\t\t<xi:include href=''%s'' xpointer=''xpointer(//Xdmf/Domain/Grid)''/>\n',filename);
        end
        fprintf(fid,'\t\t</Grid>\n');
        fprintf(fid,'\t</Domain>\n');
        fprintf(fid,'</Xdmf>\n');
    end
end


end

function g = freeze_rotors(p,g,t_pitch,t_0)
% Freeze rotors at zeros nstep coordinates, shift all stationary rows to
% best line up with them

% Find which patches are rotating
prids_rotating = []; prids_stationary = [];
for prid = 0:size(p,1)-1
    if p(prid+1,3) ~= 0
        prids_rotating = [prids_rotating prid];
        rpm = p(prid+1,3);
    else
        prids_stationary = [prids_stationary prid];
    end
end

% Get average theta for rotating patches
t_min = inf; t_max = -inf;
for prid = prids_rotating
    t_min = min(min(min(g{prid+1}.t)),t_min);
    t_max = max(max(max(g{prid+1}.t)),t_max);
end
if isinf(t_min) == 0
    t_mean_rotating = 0.5 * (t_min + t_max);
else
    t_mean_rotating = 0;
end

% Get difference between this and start of run
t_nstep = t_0 - t_mean_rotating;

% % Determine how far rotors have moved since start of run
% t_cycle = 2 * pi * rpm / (60 * frequency);
% t_nstep = t_cycle * nstep / nstep_cycle;

% Shift all patches, rotating back to the start and stationary by the same
for prid = 0:length(g)-1
    g{prid+1}.rt = g{prid+1}.rt + t_nstep * g{prid+1}.r;
    g{prid+1}.t = g{prid+1}.rt ./ g{prid+1}.r;
end

% Get average theta for both rotating and stationary patches
t_min = inf; t_max = -inf;
for prid = prids_stationary
    t_min = min(min(min(g{prid+1}.t)),t_min);
    t_max = max(max(max(g{prid+1}.t)),t_max);
end
if isinf(t_min) == 0
    t_mean_stationary = 0.5 * (t_min + t_max);
else
    t_mean_stationary = 0;
end

% Shift stationary patches to line up best
dt = t_mean_stationary - t_0;
dt = t_pitch * round(dt / t_pitch);
for prid = prids_stationary
    g{prid+1}.rt = g{prid+1}.rt - dt * g{prid+1}.r;
end

end

function g = shift_rotors(p,g,t_pitch)
% Shift all patches with non-zero rpm to line up best with stationary ones

% Find which patches are rotating
prids_rotating = []; prids_stationary = [];
for prid = 0:size(p,1)-1
    if p(prid+1,3) ~= 0
        prids_rotating = [prids_rotating prid];
    else
        prids_stationary = [prids_stationary prid];
    end
end

% Get maximum and minimum t for stationary patches
t_min = inf; t_max = -inf;
for prid = prids_stationary
    t_min = min(min(min(g{prid+1}.t)),t_min);
    t_max = max(max(max(g{prid+1}.t)),t_max);
end
if isinf(t_min) == 0
    t_mean_stationary = 0.5 * (t_min + t_max);
else
    t_mean_stationary = 0;
end

% Get maximum and minimum rt for rotating patches
t_min = inf; t_max = -inf;
for prid = prids_rotating
    t_min = min(min(min(g{prid+1}.t)),t_min);
    t_max = max(max(max(g{prid+1}.t)),t_max);
end
if isinf(t_min) == 0
    t_mean_rotating = 0.5 * (t_min + t_max);
else
    t_mean_rotating = 0;
end

% Add on shift in rotating blocks to bring in line with stationary blocks
dt = t_mean_stationary - t_mean_rotating;
dt = t_pitch * round(dt / t_pitch);
for prid = prids_rotating
    g{prid+1}.rt = g{prid+1}.rt + dt * g{prid+1}.r;
end

end

function g = flatten_probes(g)
% Project all patches onto a constant radius

% Determine new mean radius
r_min = inf; r_max = -inf;
for prid = 0:length(g)-1
    r_min = min(min(min(g{prid+1}.r)),r_min);
    r_max = max(max(max(g{prid+1}.r)),r_max);
end
r_mean = 0.5 * (r_min + r_max);
r_mean = round(r_mean * 10000) / 10000;

% Project onto new radius
for prid = 0:length(g)-1
    t = g{prid+1}.rt ./ g{prid+1}.r;
    g{prid+1}.rt = r_mean .* t;
    g{prid+1}.r = r_mean * ones(size(g{prid+1}.r));
end

end

function [] = write_hdf5(g,filename,bp,conv_cart,V_rel)
% Create new file 
fid = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');

% Create group creation property lists
gcpl = H5P.create ('H5P_LINK_CREATE');
H5P.set_create_intermediate_group (gcpl, 1);
H5F.close(fid);

% Get the file ID for the functions that write straight to the HDF5 file
flags = 'H5F_ACC_RDWR';
fapl = 'H5P_DEFAULT';
file_id = H5F.open(filename,flags,fapl);
c_file_id = onCleanup(@()H5F.close(file_id));

% Write block properties
for prid = 0:length(g)-1
    bstr = ['/probe' num2str(prid) '/'];
    bp_names = fieldnames(g{prid+1});
    for i = 1:length(bp_names)
        if sum(strcmp(bp_names{i}, bp(:,1))) == 1
            
            [ni nj nk] = size(g{prid+1}.(bp_names{i}));
            current_bp = g{prid+1}.(bp_names{i});

            j = find(strcmp(bp_names{i}, bp(:,1)) == 1);
            if strcmp(bp(j,2),'fl') == 1
                datatype = 'single';
                current_bp = single(current_bp);
            elseif strcmp(bp(j,2),'in8') == 1
                datatype = 'uint8';
                current_bp = uint8(current_bp);
            elseif strcmp(bp(j,2),'in') == 1
                datatype = 'int32';
                current_bp = int32(current_bp);
            else
                disp('Error - Datatype is missing')
                disp(bp_names{i});
            end

            h5createTS(filename, file_id, [bstr bp_names{i} '_bp'], size(current_bp), 'Datatype', datatype);
            h5writeTS(filename, file_id, [bstr bp_names{i} '_bp'], current_bp);
        end
    end
end

% Write block velocities
for prid = 0:length(g)-1
    bstr = ['/probe' num2str(prid) '/'];

    [ni, nj, nk] = size(g{prid+1}.Vx);
    if V_rel == 1
        if conv_cart == 0
            current_bp = reshape(cat(1,reshape(g{prid+1}.Vx,1,[]),reshape(g{prid+1}.Vr,1,[]),...
                reshape(g{prid+1}.Vt_rel,1,[])),[nk nj ni 3]);
        else
            current_bp = reshape(cat(1,reshape(g{prid+1}.Vx,1,[]),reshape(g{prid+1}.Vy_rel,1,[]),...
                reshape(g{prid+1}.Vz_rel,1,[])),[nk nj ni 3]);
        end
    else
        if conv_cart == 0
            current_bp = reshape(cat(1,reshape(g{prid+1}.Vx,1,[]),reshape(g{prid+1}.Vr,1,[]),...
                reshape(g{prid+1}.Vt,1,[])),[nk nj ni 3]);
        else
            current_bp = reshape(cat(1,reshape(g{prid+1}.Vx,1,[]),reshape(g{prid+1}.Vy,1,[]),...
                reshape(g{prid+1}.Vz,1,[])),[nk nj ni 3]);
        end
    end
    datatype = 'single';
    h5createTS(filename, file_id, [bstr 'V_bp'], size(current_bp), 'Datatype', datatype);
    h5writeTS(filename, file_id, [bstr 'V_bp'], current_bp);
end

% % Write block variables
% bv_names = {'Pref' 'Tref' 'rpm' 'R' 'cp' 'cv' 'ga' 'mu'};
% for prid = 0:length(g)-1
%     bstr = ['/probe' num2str(prid) '/'];
%     for i = 1:length(bv_names)
%         if isfield(g{prid+1}.bv,bv_names{i})
%         current_bv = single(g{prid+1}.bv.(bv_names{i}));
% 
%         datatype = 'single';
%         h5createTS(filename, file_id, [bstr bv_names{i} '_bv'], size(current_bv), 'Datatype', datatype);
%         h5writeTS(filename, file_id, [bstr bv_names{i} '_bv'], current_bv);
%         end
%     end
% end

end


function [] = write_xdmf(g,filename,bp,time,nstep,conv_cart)
% Write xdmf file for paraview

attributes = bp;
if conv_cart == 0
    geometry = { 'x' 'Float' ; 'r' 'Float' ; 'rt' 'Float'};
    velocities = { 'Vx' 'Float' 4 ; 'Vr' 'Float' 4 ; 'Vt_rel' 'Float' 4};
else
    geometry = { 'x' 'Float' ; 'y' 'Float' ; 'z' 'Float'};
    velocities = { 'Vx' 'Float' 4 ; 'Vy' 'Float' 4 ; 'Vz' 'Float' 4};
end
for a = size(attributes,1):-1:1
    for v = 1:size(velocities,1)
        if strcmp(attributes{a,1},velocities{v,1}) == 1
            attributes(a,:) = [];
        end
    end
end
bv_names = {'Pref' 'Tref' 'rpm' 'R' 'cp' 'cv' 'ga' 'mu'};
    
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

% remove directory info if attached to hdf5 filename
A = regexp(filename,'/');
if isempty(A) == 0
    filename = filename(A(end)+1:end);
end

% Open the file for printing
fid = fopen(xdmf_filename,'w');

% Print headers
fprintf(fid,'<?xml version=''1.0'' ?>\n');
fprintf(fid,'<!DOCTYPE Xdmf SYSTEM ''Xdmf.dtd'' []>\n');
fprintf(fid,'<Xdmf>\n');
fprintf(fid,'\t<Domain>\n');

fprintf(fid,'\t\t<Grid Name=''step%d'' GridType=''Collection'' CollectionType=''Spatial''>\n',nstep);
fprintf(fid,'\t\t<Time Type=''Single'' Value = ''%10.8f''/>\n',time);

% Print block info
for prid = 0:length(g)-1
    fprintf(fid,'\t\t\t<Grid Name=''probe%d'' GridType=''Uniform''>\n',prid);

    [ni, nj, nk] = size(g{prid+1}.x);
%     fprintf(fid,'\t\t\t\t<Topology Type=''3DSMESH'' Dimensions=''%d %d %d''/>\n', nk, nj, ni);
    fprintf(fid,'\t\t\t\t<Topology Type=''2DSMESH'' Dimensions=''%d %d''/>\n', nj, ni);

    % Print mesh
    fprintf(fid,'\t\t\t\t<Geometry Type=''X_Y_Z''>\n');
    for n = 1:size(geometry,1)
        coord = geometry{n,1};
        fprintf(fid,['\t\t\t\t\t<DataItem DataType=''Float'' Precision=''4'' Dimensions=''%d %d'' '...
            'Format=''HDF''>\n'], nj, ni);
        fprintf(fid,'\t\t\t\t\t\t%s:/probe%d/%s_bp\n',filename,prid,coord);
        fprintf(fid,'\t\t\t\t\t</DataItem>\n');
    end
    fprintf(fid,'\t\t\t\t</Geometry>\n');

    % Print block properties
    varnames = fieldnames(g{prid+1});
    for v = 1:length(varnames)
        varname = varnames{v};
        for a = 1:size(attributes,1)
            attribute = attributes{a,1};
            if strcmp(attribute,varname) == 1
                fprintf(fid,['\t\t\t\t<Attribute Name=''%s'' AttributeType=''Scalar'' '...
                    'Center=''Node''>\n'],attribute);
                fprintf(fid,['\t\t\t\t\t<DataItem Dimensions=''%d %d'' '...
                    'NumberType=''%s'' Precision=''%d'' Format=''HDF''>\n'], nj, ni,...
                    attributes{a,2}, attributes{a,3});
                fprintf(fid,'\t\t\t\t\t\t%s:/probe%d/%s_bp\n',filename,prid,attribute);
                fprintf(fid,'\t\t\t\t\t</DataItem>\n');
                fprintf(fid,'\t\t\t\t</Attribute>\n');
            end
        end
    end

    % Print velocity vector
    fprintf(fid,['\t\t\t\t<Attribute Name=''%s'' AttributeType=''Vector'' '...
        'Center=''Node''>\n'],'V');
    fprintf(fid,['\t\t\t\t\t<DataItem Dimensions=''%d %d %d'' '...
        'NumberType=''%s'' Precision=''%d'' Format=''HDF''>\n'], nj, ni, 3,...
        'Float', 4);        
    fprintf(fid,'\t\t\t\t\t\t%s:/probe%d/%s_bp\n',filename,prid,'V');
    fprintf(fid,'\t\t\t\t\t</DataItem>\n');
    fprintf(fid,'\t\t\t\t</Attribute>\n');
   
%     % Print block variables
%     for v = 1:length(bv_names)
%         varname = bv_names{v};
%         fprintf(fid,['\t\t\t\t<Attribute Name=''%s'' AttributeType=''Scalar'' '...
%             'Center=''Grid''>\n'],varname);
%         fprintf(fid,['\t\t\t\t\t<DataItem Dimensions=''1 1'' '...
%             'NumberType=''Float'' Precision=''4'' Format=''HDF''>\n']);
%         fprintf(fid,'\t\t\t\t\t\t%s:/probe%d/%s_bv\n',filename,prid,varname);
%         fprintf(fid,'\t\t\t\t\t</DataItem>\n');
%         fprintf(fid,'\t\t\t\t</Attribute>\n');
%     end
    
    fprintf(fid,'\t\t\t</Grid>\n');
end
fprintf(fid,'\t\t</Grid>\n');

% fprintf(fid,'\t\t</Grid>\n');
fprintf(fid,'\t</Domain>\n');
fprintf(fid,'</Xdmf>\n');
fprintf(fid,'\n');

fclose(fid);


end


function h5writeTS(Filename,file_id,Dataset,Data,varargin)
%H5WRITE Write to HDF5 dataset.
%   H5WRITE(FILENAME,DATASETNAME,DATA) writes to an entire dataset.
%
%   H5WRITE(FILENAME,DATASETNAME,DATA,START,COUNT) writes a subset of
%   data.  START is the index of the first element to be written and is
%   one-based.  COUNT defines how many elements to write along each
%   dimension.  An extendible dataset will be extended along any unlimited
%   dimensions if necessary.
%
%   H5WRITE(FILENAME,DATASETNAME,DATA,START,COUNT,STRIDE) writes a
%   hyperslab of data.  STRIDE is the inter-element spacing along each
%   dimension.   STRIDE always defaults to a vector of ones if not
%   supplied.
%
%   Only floating point and integer datasets are supported.
%
%   Example:  Write to an entire dataset.
%       h5create('myfile.h5','/DS1',[10 20]);
%       h5disp('myfile.h5');
%       mydata = rand(10,20);
%       h5writeTS('myfile.h5', '/DS1', mydata);
%
%   Example:  Write a hyperslab to the last 5-by-7 block of a dataset.
%       h5create('myfile.h5','/DS2',[10 20]);
%       h5disp('myfile.h5');
%       mydata = rand(5,7);
%       h5writeTS('myfile.h5','/DS2',mydata,[6 14],[5 7]);
%
%   Example:  Append to an unlimited dataset.
%       h5create('myfile.h5','/DS3',[20 Inf],'ChunkSize',[5 5]);
%       h5disp('myfile.h5');
%       for j = 1:10
%            data = j*ones(20,1);
%            start = [1 j];
%            count = [20 1];
%            h5writeTS('myfile.h5','/DS3',data,start,count);
%       end
%       h5disp('myfile.h5');
%
%   See also H5CREATE, H5DISP, H5READ, H5WRITEATT.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2011/05/17 02:24:56 $

p.Results.Filename = Filename;
p.Results.Dataset = Dataset;
p.Results.Data = Data;
p.Results.start = [];
p.Results.count = [];
p.Results.stride = [];

% p = inputParser;
% p.addRequired('Filename',@(x)ischar(x));
% p.addRequired('Dataset',@ischar);
% p.addRequired('Data',@isnumeric);
% p.addOptional('start',[], @(x) isa(x,'double') && isrow(x) && ~any(x<=0));
% p.addOptional('count',[], @(x) isa(x,'double') && isrow(x) && ~any(x<=0));
% p.addOptional('stride',[],@(x) isa(x,'double') && isrow(x) && ~any(x<=0));
%
% p.parse(Filename,Dataset,Data,varargin{:});
%
% validate_write_arguments(p.Results);

dapl = 'H5P_DEFAULT';
dxpl = 'H5P_DEFAULT';

% if ~exist(Filename,'file')
%     error(message('MATLAB:imagesci:h5write:fileDoesNotExist', Filename, Dataset));
% end
% if ~H5F.is_hdf5(Filename)
%     error(message('MATLAB:imagesci:h5write:notHDF5', Filename));
% end

% file_id = H5F.open(Filename,flags,fapl);
% c_file_id = onCleanup(@()H5F.close(file_id));

try
    dataset_id   = H5D.open(file_id,Dataset,dapl);
catch me
    if strcmp(me.identifier,'MATLAB:imagesci:hdf5lib:libraryError')
        error(message('MATLAB:imagesci:h5write:datasetDoesNotExist', Dataset));
    else
        rethrow(me);
    end
end
c_dataset_id = onCleanup(@()H5D.close(dataset_id));


datatype_id   = H5D.get_type(dataset_id);
c_datatype_id = onCleanup(@()H5T.close(datatype_id));

space_id   = H5D.get_space(dataset_id);
c_space_id = onCleanup(@()H5S.close(space_id));

[offset,count,stride] = determine_indexing(space_id,p.Results);

filespace_id   = get_filespace(dataset_id,offset,count,stride);
c_filespace_id = onCleanup(@()H5S.close(filespace_id));

memspace_id = H5S.create_simple(numel(count),fliplr(count),[]);
c_memspace_id = onCleanup(@()H5S.close(memspace_id));


H5D.write(dataset_id,'H5ML_DEFAULT',memspace_id,filespace_id,dxpl,Data);

return;

end

%--------------------------------------------------------------------------
function [offset,count,stride] = determine_indexing(space_id,options)

[offset,count,stride] = deal(options.start,options.count,options.stride);

[ndims,dims] = H5S.get_simple_extent_dims(space_id);
dims = fliplr(dims);

sz = size(options.Data);
% Check the length of the size vector against the rank of the dataset's
% dataspace.  If the size vector lags, it is probably because of trailing
% singletons that MATLAB automatically trims.  We want to augment the size
% vector in that case.
if numel(sz) < ndims
    sz = [sz ones(1,(ndims-numel(sz)))];
end

if isempty(offset)
    % The user gave us no information, so we assume we are writing the
    % entire dataset.
    offset = zeros(1,ndims);
    count = dims;
    
    
    stride = ones(1,ndims);
    
    if numel(options.Data) ~= prod(count)
        error(message('MATLAB:imagesci:h5write:fullDatasetDataMismatch'));
    end
    
    if ndims > 1
        % 1D datasets get a free pass.
        if any(sz - count)
            error(message('MATLAB:imagesci:h5write:fullDatasetDimsMismatch'));
        end
    end
    return
    
elseif isempty(stride)
    
    % The low-level interface requires zero-based indexing.
    offset = offset-1;
    
    % Supply a default stride.
    stride = ones(1,numel(offset));
    
else
    
    % Stride was provided.
    % The low-level interface requires zero-based indexing.
    offset = offset-1;
    
end

% Special error checking if START and COUNT are provided.
if (numel(options.Data) ~= prod(count))
    error(message('MATLAB:imagesci:h5write:datasetCountMismatch'));
end
if ndims > 1
    if any(sz - count)
        error(message('MATLAB:imagesci:h5write:datasetDimsMismatch'));
    end
end


return
end

%--------------------------------------------------------------------------
function filespace_id = get_filespace(dataset_id,start,count,stride)

filespace_id = H5D.get_space(dataset_id);

% Get the size of the dataset so we know if we have to extend it.
[~, spaceDims,maxDims] = H5S.get_simple_extent_dims(filespace_id);
spaceDims = fliplr(spaceDims);
maxDims = fliplr(maxDims);

boundsEnd = start + (count-1).*stride;

% Define the hyperslab selection.
H5S.select_hyperslab(filespace_id, 'H5S_SELECT_SET', ...
    fliplr(start), fliplr(stride), ...
    fliplr(count), ones(1,length(start)));


% Do any of the bounding box dimensions exceed the current dimensions
% of the dataset?
dataset_extents_exceeded = (boundsEnd > (spaceDims-1));
if ~any(dataset_extents_exceeded)
    % No, no need to extend.
    return;
end

% The hyperslab selection does exceed the current dataset size, so now
% we need to figure out if we can extend.
extend_idx = dataset_extents_exceeded;
if any(maxDims(extend_idx)>-1)
    H5S.close(filespace_id);
    error(message('MATLAB:imagesci:h5write:cannotExtend', H5I.get_name( dataset_id ), H5I.get_name( dataset_id )));
end

% Figure out the new dimensions, extend the dataset, and return the new
% space ID after selecting the hyperslab.
new_dims = max(spaceDims,boundsEnd+1);
H5S.close(filespace_id);
H5D.extend(dataset_id,fliplr(new_dims));


% Retrieve the new file space ID and make the hyperslab selection.
filespace_id = H5D.get_space(dataset_id);
H5S.select_hyperslab(filespace_id, 'H5S_SELECT_SET', ...
    fliplr(start), fliplr(stride), ...
    fliplr(count), ones(1,length(start)));

return

end

%--------------------------------------------------------------------------
function dataspace_id = createDataspaceId(attvalue)
% Setup the dataspace ID.  This just depends on how many elements the
% attribute actually has.

if isempty(attvalue)
    dataspace_id = H5S.create('H5S_NULL');
    return;
elseif ischar(attvalue)
    if isrow(attvalue)
        dataspace_id = H5S.create('H5S_SCALAR');
        return
    else
        error(message('MATLAB:imagesci:h5writeatt:badStringSize'));
    end
else
    if ( ndims(attvalue) == 2 ) && ( any(size(attvalue) ==1) )
        rank = 1;
        dims = numel(attvalue);
    else
        % attribute is a "real" 2D value.
        rank = ndims(attvalue);
        dims = fliplr(size(attvalue));
    end
end
dataspace_id = H5S.create_simple(rank,dims,dims);

end

%--------------------------------------------------------------------------
function datatype_id = createDatatypeId ( attvalue )
% We need to choose an appropriate HDF5 datatype based upon the attribute
% data.
switch class(attvalue)
    case 'single'
        datatype_id = H5T.copy('H5T_NATIVE_FLOAT');
    case 'int32'
        datatype_id = H5T.copy('H5T_NATIVE_INT');
    case 'uint8'
        datatype_id = H5T.copy('H5T_NATIVE_UCHAR');
    case 'double'
        datatype_id = H5T.copy('H5T_NATIVE_DOUBLE');        
    case 'int64'
        datatype_id = H5T.copy('H5T_NATIVE_LLONG');
    case 'uint64'
        datatype_id = H5T.copy('H5T_NATIVE_ULLONG');
    case 'uint32'
        datatype_id = H5T.copy('H5T_NATIVE_UINT');
    case 'int16'
        datatype_id = H5T.copy('H5T_NATIVE_SHORT');
    case 'uint16'
        datatype_id = H5T.copy('H5T_NATIVE_USHORT');
    case 'int8'
        datatype_id = H5T.copy('H5T_NATIVE_SCHAR');
    case 'char'
        datatype_id = H5T.copy('H5T_C_S1');
        if ~isempty(attvalue)
            % Don't do this when working with empty strings.
            H5T.set_size(datatype_id,numel(attvalue));
        end
        H5T.set_strpad(datatype_id,'H5T_STR_NULLTERM');
    otherwise
        error(message('MATLAB:imagesci:h5writeatt:unsupportedAttributeDatatype', class( attvalue )));
end
return

end

function h5createTS(Filename,file_id,Dataset,Size,varargin)
%H5CREATE  Create HDF5 dataset.
%   H5CREATE(FILENAME,DATASETNAME,SIZE,Param1,Value1, ...) creates an HDF5
%   dataset with name DATASETNAME and with extents given by SIZE in the
%   file given by FILENAME.  If DATASETNAME is a full path name, all
%   intermediate groups are created if they don't already exist.  If
%   FILENAME does not already exist, it is created.
%
%   Elements of SIZE should be Inf in order to specify an unlimited extent.
%
%   Parameter Value Pairs
%   ---------------------
%       'Datatype'   - May be one of 'double', 'single', 
%                      'uint64', 'int64', 'uint32', 'int32', 'uint16', 
%                      'int16', 'uint8', or 'int8'.  Defaults to 'double'.
%       'ChunkSize'  - Defines chunking layout.  Default is not chunked.
%       'Deflate'    - Defines gzip compression level (0-9).  Default is 
%                      no compression.
%       'FillValue'  - Defines the fill value for numeric datasets.
%       'Fletcher32' - Turns on the Fletcher32 checksum filter.  Default 
%                      value is false.
%       'Shuffle'    - Turns on the Shuffle filter.  Default value is
%                      false.
%
%   Example:  create a fixed-size 100x200 dataset.
%       h5create('myfile.h5','/myDataset1',[100 200]);
%       h5disp('myfile.h5');
%
%   Example:  create a single precision 1000x2000 dataset with a chunk size
%   of 50x80.  Apply the highest level of compression possible.
%       h5create('myfile.h5','/myDataset2',[1000 2000], 'Datatype','single', ...
%                'ChunkSize',[50 80],'Deflate',9);
%       h5disp('myfile.h5');
%
%   Example:  create a two-dimensional dataset that is unlimited along the
%   second extent.
%       h5create('myfile.h5','/myDataset3',[200 Inf],'ChunkSize',[20 20]);
%       h5disp('myfile.h5');
%
%   See also:  h5read, h5write, h5info, h5disp.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2011/05/17 02:24:50 $


p = inputParser;
p.addRequired('Filename',@(x)ischar(x));
p.addRequired('Dataset', @(x)ischar(x) & ~isempty(x));  

p.addRequired('Size', @(x) isa(x,'double') && isrow(x) && ~any(x<0));
p.addParamValue('Datatype','double', @(x) ischar(x) && isrow(x) );
p.addParamValue('ChunkSize', [], ...
    @(x) isa(x,'double') && isrow(x) && ~any(x<0));
p.addParamValue('Deflate', [], ...
    @(x) isa(x,'double') && isscalar(x) && (x>=0) && (x<=9));                                      
p.addParamValue('FillValue',[], @(x) isnumeric(x) && isscalar(x));
p.addParamValue('Fletcher32',false, ...
    @(x) (isa(x,'double') || islogical(x)) && isscalar(x));
p.addParamValue('Shuffle',false, ...
    @(x) (isa(x,'double') || islogical(x)) && isscalar(x));

p.parse(Filename,Dataset,Size,varargin{:});
options = validate_options(p.Results);
options.file_id = file_id;
create_dataset(options);

return

end
%--------------------------------------------------------------------------
function options = validate_options(options)

% Give a better error message than the low level library would.
if ~isempty(options.Deflate) && isempty(options.ChunkSize)
    error(message('MATLAB:imagesci:h5create:deflateRequiresChunking'));
end

% Setup Extendable.  Either
options.Extendable = false(1,numel(options.Size));
options.Extendable(isinf(options.Size)) = true;
options.Extendable(options.Size == 0) = true;

% Force Shuffle and Fletcher32 options to be logical.
if isnumeric(options.Fletcher32)
    options.Fletcher32 = logical(options.Fletcher32);
end  
if isnumeric(options.Shuffle)
    options.Shuffle = logical(options.Shuffle);
end  
if (options.Fletcher32 || options.Shuffle) && isempty(options.ChunkSize)
     error(message('MATLAB:imagesci:h5create:filterRequiresChunking')); 
end
  
if ~isempty(options.FillValue) && ~strcmp(options.Datatype,class(options.FillValue))
    error(message('MATLAB:imagesci:h5create:datasetFillValueMismatch', class( options.FillValue ), options.Datatype));
end

% Make sure that chunk size does not exceed dataset size.  After that,
% reset any Infs to zero before continuing.  The initial size of an
% unlimited extent must be zero to begin with.

if ~isempty(options.ChunkSize)
    if numel(options.ChunkSize) ~= numel(options.Size)
        error(message('MATLAB:imagesci:h5create:chunkSizeDatasetSizeMismatch'));
    end
    if any((options.ChunkSize - options.Size) > 0)
        error(message('MATLAB:imagesci:h5create:chunkSizeLargerThanDataset'));
    end
end
options.Size(isinf(options.Size)) = 0;

if ( ~isempty(options.Extendable) ) 
    if any(options.Extendable) && isempty(options.ChunkSize)
        error(message('MATLAB:imagesci:h5create:extendibleRequiresChunking'));   
    end
end


% If the file exists, check that it is an HDF5 file.
% if exist(options.Filename,'file')
%     fid = fopen(options.Filename);
%     fullpath = fopen(fid);
%     fclose(fid);
%     if ~H5F.is_hdf5(fullpath)
%         error(message('MATLAB:imagesci:h5create:notHDF5', options.Filename));
%     end
%     options.Filename = fullpath;
% end

if options.Dataset(1) ~= '/'
    error(message('MATLAB:imagesci:h5create:notFullPathName'));
end

end

%--------------------------------------------------------------------------
function create_dataset(options)

% if exist(options.Filename,'file')
%     fid = H5F.open(options.Filename,'H5F_ACC_RDWR','H5P_DEFAULT');
%     file_was_created = false;
% else
%     fid = H5F.create(options.Filename,'H5F_ACC_TRUNC','H5P_DEFAULT', ...
%         'H5P_DEFAULT');
%     file_was_created = true;
% end



% Does the dataset already exist?
% try
%     dset = H5D.open(options.file_id, options.Dataset);
%     H5D.close(dset);
%     error(message('MATLAB:imagesci:h5create:datasetAlreadyExists', options.Dataset));
% catch me
%     if strcmp(me.identifier,'MATLAB:imagesci:h5create:datasetAlreadyExists')
%         rethrow(me)
%     end
% end

try    
    switch(options.Datatype)
        case 'single'
            datatype = 'H5T_NATIVE_FLOAT';
        case 'uint32'
            datatype = 'H5T_NATIVE_UINT';  
        case 'uint8'
            datatype = 'H5T_NATIVE_UCHAR';
        case 'uint64'
            datatype = 'H5T_NATIVE_UINT64';
        case 'double'
            datatype = 'H5T_NATIVE_DOUBLE';
        case 'int64'
            datatype = 'H5T_NATIVE_INT64';
        case 'int32'
            datatype = 'H5T_NATIVE_INT';
        case 'uint16'
            datatype = 'H5T_NATIVE_USHORT';
        case 'int16'
            datatype = 'H5T_NATIVE_SHORT';
        case 'int8'
            datatype = 'H5T_NATIVE_CHAR';
        otherwise
            error(message('MATLAB:imagesci:h5create:unrecognizedDatatypeString', options.Datatype));
    end
    
    
    % Set the maxdims parameter to take into account any extendable
    % dimensions.
    maxdims = options.Size;
    if any(options.Extendable)
        unlimited = H5ML.get_constant_value('H5S_UNLIMITED');
        maxdims(options.Extendable) = unlimited;
    end
    
    % Create the dataspace.
    space_id = H5S.create_simple(numel(options.Size), ...
            fliplr(options.Size), fliplr(maxdims));
    cspace_id = onCleanup(@()H5S.close(space_id));
    
    
    % If the dataset is buried a few groups down, then we want to create 
    % all intermediate groups.
    lcpl = H5P.create('H5P_LINK_CREATE');
    clcpl = onCleanup(@()H5P.close(lcpl));
    H5P.set_create_intermediate_group(lcpl,1);
    
    
    dcpl = construct_dataset_creation_property_list(options);
    cdcpl = onCleanup(@()H5P.close(dcpl));
    dapl = 'H5P_DEFAULT';
    
    dset_id = H5D.create(options.file_id,options.Dataset,datatype,space_id,lcpl,dcpl,dapl);
    
catch me
    H5F.close(options.file_id);
%     if file_was_created
%         delete(options.Filename);
%     end
    rethrow(me);
end

H5D.close(dset_id);
return

end

%--------------------------------------------------------------------------
function dcpl = construct_dataset_creation_property_list(options)
% Setup the DCPL - dataset create property list.


dcpl = H5P.create('H5P_DATASET_CREATE');

% Modify the dataset creation property list for the shuffle filter if
% so ordered.
if options.Shuffle
    H5P.set_shuffle(dcpl);
end

% Modify the dataset creation property list for possible chunking and
% deflation.
if ~isempty(options.ChunkSize)
    H5P.set_chunk(dcpl,fliplr(options.ChunkSize));
end
if ~isempty(options.Deflate)
    H5P.set_deflate(dcpl,options.Deflate);
end

% Modify the dataset creation property list for a possible fill value.
if ~isempty(options.FillValue)
    switch(options.Datatype)
        case 'double'
            filltype = 'H5T_NATIVE_DOUBLE';
            fv = double(options.FillValue);
        case 'single'
            filltype = 'H5T_NATIVE_FLOAT';
            fv = single(options.FillValue);
        case 'uint64'
            filltype = 'H5T_NATIVE_UINT64';
            fv = uint64(options.FillValue);
        case 'int64'
            filltype = 'H5T_NATIVE_INT64';
            fv = int64(options.FillValue);
        case 'uint32'
            filltype = 'H5T_NATIVE_UINT';
            fv = uint32(options.FillValue);
        case 'int32'
            filltype = 'H5T_NATIVE_INT';
            fv = int32(options.FillValue);
        case 'uint16'
            filltype = 'H5T_NATIVE_USHORT';
            fv = uint16(options.FillValue);
        case 'int16'
            filltype = 'H5T_NATIVE_SHORT';
            fv = int16(options.FillValue);
        case 'uint8'
            filltype = 'H5T_NATIVE_UCHAR';
            fv = uint8(options.FillValue);
        case 'int8'
            filltype = 'H5T_NATIVE_CHAR';
            fv = int8(options.FillValue);
        otherwise
            H5P.close(dcpl);
            error(message('MATLAB:imagesci:h5create:badFillValueType'));
    end
    H5P.set_alloc_time(dcpl,'H5D_ALLOC_TIME_EARLY');
    H5P.set_fill_value(dcpl,filltype,fv);
end

% Modify the dataset creation property list for the fletcher32 filter if
% so ordered.
if options.Fletcher32
    H5P.set_fletcher32(dcpl);
end

end
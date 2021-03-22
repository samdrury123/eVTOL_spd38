function C = ts_cut_stage(directory,filename,xs,force_read)
% Take important cuts from a TURBOSTREAM solution of a stage

% Default to not force reading of whole file
if exist('force_read','var') == 0
    force_read = 0;
end

% Default to hdf5 name
filename = strrep(filename,'.mat','.hdf5');
mat_name = strrep(filename,'.hdf5','.mat');
cuts_name = strrep(filename,'.hdf5','_cuts.mat');

% Get dates of all filenames
A_hdf5 = dir([directory filename]);
A_mat = dir([directory mat_name]); 
A_cuts = dir([directory cuts_name]);

% Make sure newest modified file is always read
if numel(A_cuts) == 1 
    if (numel(A_mat) == 1 && A_mat.datenum > A_cuts.datenum) || ...
            (numel(A_hdf5) == 1 && A_hdf5.datenum > A_cuts.datenum)
        delete([directory cuts_name]);
    end
end

% Read and cut the whole solution if cut file doesn't exist
if exist([directory cuts_name],'file') == 0 || force_read == 1
    
    % Read in meridional lines if possible
    if exist([directory 'xr_hub.mat'],'file') ~= 0
        load([directory 'xr_hub.mat'])
        load([directory 'xr_cas.mat'])
    else
        xr_hub = []; xr_cas = [];
    end

    % Read in whole TURBOSTREAM file
    if exist([directory filename],'file') == 0 && exist([directory mat_name],'file') == 0
        C = [];
        disp([filename ' not found'])
        return
    else
        g = ts_read_hdf5([directory filename]);
    end
    
    % Take inlet and outlet cuts of the whole machine
    C = ts_cut_rows(g,xr_hub,xr_cas,1,1);
    
    % Take extra cuts at traverse planes
    C{1}.axial{1} = ts_structunstruct_cut(g,[xs(1) 0 ; xs(1) 1],[],2,1);
    C{2}.axial{1} = ts_structunstruct_cut(g,[xs(2) 0 ; xs(2) 1],[],2,1);
    
    % Save cuts
    save([directory cuts_name],'C');
    
else
    
    % Load the cuts from file
    load([directory cuts_name]);

end


end
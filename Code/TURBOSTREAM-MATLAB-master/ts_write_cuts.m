function [g_cut] = ts_write_cuts(job, cuts, use_avg, plot_cuts, remove_files)
% Read in Turbostream .hdf5 file and save unstructured and structured cuts
% to reduce memory usage.
if exist('use_avg','var') == 0
    use_avg = 0;
end

if exist('plot_cuts','var') == 0
    plot_cuts = 0;
end

if exist('remove_files','var') == 0
    remove_files = 0;
end

if use_avg == 0
    if exist([job.directory job.outFile],'file') ~= 0
        if exist(strrep([job.directory job.outFile],'.hdf5','.mat'),'file') == 0
            g = ts_read_hdf5([job.directory job.outFile]);
        else
            load(strrep([job.directory job.outFile],'.hdf5','.mat'))
        end
        outfile = job.outFile;
    else
        disp([job.outFile ' does not exist']);
        g_cut = [];
        return
    end
else
    if exist([job.directory job.avgFile],'file') ~= 0
        if exist(strrep([job.directory job.avgFile],'.hdf5','.mat'),'file') == 0
            g = ts_read_hdf5([job.directory job.avgFile]);
        else
            load(strrep([job.directory job.avgFile],'.hdf5','.mat'))
        end
        outfile = job.avgFile;
    else
        disp([job.avgFile ' does not exist']);
        g_cut = [];
        return
    end
end

% if exist('bids','var') == 0
bids = ts_get_bids(g);
% end
%g = ts_secondary(g);

if plot_cuts == 1
    figure()
    hold on
    axis equal
    colormap bone
    xlabel('X'); ylabel('R'); zlabel('RT')
end

% Get bids of inlet and exit
for bid = bids
    for pid = ts_get_pids(g,bid)
        if g{bid+1}.patch{pid+1}.attribute.kind == 0
            bid_inlet = bid;
        end
        if g{bid+1}.patch{pid+1}.attribute.kind == 1
            bid_outlet = bid;
        end
    end
end

% % Read inlet cut
% if exist('bid_inlet','var') ~= 0
%     g_cut.Inlet = ts_structured_cut(g,bid_inlet,1,1,1,'en',1,'en');
%     if plot_cuts == 1
%         mesh(double(g_cut.Inlet.x),double(g_cut.Inlet.r),double(g_cut.Inlet.rt),double(g_cut.Inlet.rovx));
%     end
% end
% 
% % Read Exit Cut
% if exist('bid_outlet','var') ~= 0
%     g_cut.Outlet = ts_structured_cut(g,bid_outlet,'en','en',1,'en',1,'en');
%     if plot_cuts == 1
%         mesh(double(g_cut.Outlet.x),double(g_cut.Outlet.r),double(g_cut.Outlet.rt),double(g_cut.Outlet.rovx));
%     end
% end    

% Read other structured cuts
if isfield(cuts,'structured') == 1
    for n = 1:length(cuts.structured)
        bid = cuts.structured{n}.bid;
        i1 = cuts.structured{n}.i1;
        i2 = cuts.structured{n}.i2;
        j1 = cuts.structured{n}.j1;
        j2 = cuts.structured{n}.j2;
        k1 = cuts.structured{n}.k1;
        k2 = cuts.structured{n}.k2;
        g_cut.(['S' num2str(n)]) = ts_structured_cut(g,bid,i1,i2,j1,j2,k1,k2);
        if plot_cuts == 1
            mesh(double(g_cut.(['S' num2str(n)]).x),double(g_cut.(['S' num2str(n)]).r),...
                double(g_cut.(['S' num2str(n)]).rt),double(g_cut.(['S' num2str(n)]).rovx))
        end
    end
end

% Read Unstructured Cuts
if isfield(cuts,'unstructured') == 1
    for n = 1:length(cuts.unstructured)
        bids = cuts.unstructured{n}.bids;
        splits = cuts.unstructured{n}.splits;
        xr = cuts.unstructured{n}.xr;
        g_cut.(['C' num2str(n)]) = ts_structunstruct_cut(g,xr,bids,2,1,splits);
        if plot_cuts == 1
            mesh(double(g_cut.(['C' num2str(n)]).x),double(g_cut.(['C' num2str(n)]).r),...
                double(g_cut.(['C' num2str(n)]).rt),double(g_cut.(['C' num2str(n)]).rovx))
        end
    end
end

% Read Blade to Blade Data
if isfield(cuts,'radial') == 1
    for n = 1:length(cuts.radial)
        bid = cuts.radial{n}.bid;
        bids = cuts.radial{n}.bids;
        R = cuts.radial{n}.r;
        r = squeeze(g{bid+1}.r(1,:,1));
        r_nondim = (r - min(r)) / (max(r) - min(r));
        [~, J] = min(abs(r_nondim - R));
        for m = 1:length(bids)
            bid = bids(m);
            g_cut.(['R' num2str(100*R)]){m} = ts_structured_cut(g,bid,1,'en',J,J,1,'en');
            if plot_cuts == 1
                mesh(double(g_cut.(['R' num2str(100*R)]){m}.x),double(g_cut.(['R' num2str(100*R)]){m}.r),...
                    double(g_cut.(['R' num2str(100*R)]){m}.rt),double(g_cut.(['R' num2str(100*R)]){m}.rovx))
            end
        end
        
    end
end

% Save cuts into .mat file
save([job.directory strrep(strrep(outfile,'.hdf5','.mat'),'.mat','_cut.mat')], 'g_cut','-v7.3');

% Plot out cuts if required
if plot_cuts == 1
    disp('%% Output files not deleted %%')
    remove_files = 0;
end
    
% Remove all other files
if remove_files == 1
    files = {job.inFile; job.outFile; job.xdmfFile ; job.avgFile ; job.xavgFile ; job.logFile ;...
        strrep(job.inFile,'.xdmf','.hdf5')};
    if isfield(job,'xinFile')
        files = {files ; job.xinFile};
    end
    for fid = 1:length(files)
        filename = files{fid};
        delete([job.directory filename]);
    end
end

end
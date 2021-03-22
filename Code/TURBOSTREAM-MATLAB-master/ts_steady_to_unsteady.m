function job = ts_steady_to_unsteady(g,frac,job,plot_stuff)
% Write and run a python script to create an unsteady TS input file

% Default to plot
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end


%% Scale and duplication numbers

% Get blade numbers
nblades = zeros(size(g));
for bid = ts_get_bids(g);
    nblades(bid+1) = g{bid+1}.bv.nblade;
end

% Find annulus fraction scale factors and duplication numbers
nb = nblades * frac;
d = round(nb);
s = nb ./ d;


%% Find periodic patches

pps = [];

% Loop over all periodic patches
for bid = ts_get_bids(g)
    for pid = ts_get_pids(g,bid)
        if g{bid+1}.patch{pid+1}.attribute.kind == 5
            
            % Find average rt of one side of patch
            i = g{bid+1}.patch{pid+1}.attribute.ist+1:g{bid+1}.patch{pid+1}.attribute.ien;
            j = g{bid+1}.patch{pid+1}.attribute.jst+1:g{bid+1}.patch{pid+1}.attribute.jen;
            k = g{bid+1}.patch{pid+1}.attribute.kst+1:g{bid+1}.patch{pid+1}.attribute.ken;

            t_1 = g{bid+1}.rt(i,j,k) ./ g{bid+1}.r(i,j,k); t_1 = mean(t_1(:));

            % Find average rt of opposite side of patch
            nxbid = g{bid+1}.patch{pid+1}.attribute.nxbid;
            nxpid = g{bid+1}.patch{pid+1}.attribute.nxpid;
            i = g{nxbid+1}.patch{nxpid+1}.attribute.ist+1:...
                g{nxbid+1}.patch{nxpid+1}.attribute.ien;
            j = g{nxbid+1}.patch{nxpid+1}.attribute.jst+1:...
                g{nxbid+1}.patch{nxpid+1}.attribute.jen;
            k = g{nxbid+1}.patch{nxpid+1}.attribute.kst+1:...
                g{nxbid+1}.patch{nxpid+1}.attribute.ken;

            t_2 = g{nxbid+1}.rt(i,j,k) ./ g{nxbid+1}.r(i,j,k); t_2 = mean(t_2(:));
            
            % Blade pitch
            t_pitch = 2 * pi / double(g{bid+1}.bv.nblade);
            
            % Check distance between patches
            if t_2 - t_1 > 0.5 * t_pitch
                pps = [pps ; bid pid];
            end
        end
    end
end

% Plot the periodics in matching colours accross the passage
g = ts_secondary(g);
if plot_stuff == 1
    
    % Plot the surfaces for reference
    h = figure(); hold on; axis equal;
    ts_plot_surface(g,[0 0 0],h,1,1,ts_get_bids(g));
    
    % Generate colour array
    C = jet(size(pps,1));
    
    % Loop over all periodic pairs
    for n = 1:size(pps,1)
        bid = pps(n,1); pid = pps(n,2);

        % Get patch indices
        i = g{bid+1}.patch{pid+1}.attribute.ist+1:g{bid+1}.patch{pid+1}.attribute.ien;
        j = g{bid+1}.patch{pid+1}.attribute.jst+1:g{bid+1}.patch{pid+1}.attribute.jen;
        k = g{bid+1}.patch{pid+1}.attribute.kst+1:g{bid+1}.patch{pid+1}.attribute.ken;

        % Get coordinates
        x = squeeze(g{bid+1}.x(i,j,k));
        y = squeeze(g{bid+1}.y(i,j,k));
        z = squeeze(g{bid+1}.z(i,j,k));

        % Get colours and plot surfaces
        c = squeeze(repmat(reshape(C(n,:),[1 1 1 3]),size(x)));
        surf(x,z,y,c,'FaceAlpha',0.5);
    end
end


%% Assemble the multi-passage mesh

% % Scale the blades to occupy the right pitch
% for bid = ts_get_bids(g)
%     g{bid+1}.rt = g{bid+1}.rt * s(bid+1);
% end
% 
% % Create matrix of blocks to repeat
% G = cell(max(nb),length(g));
% for bid = ts_get_bids(g)
%     G(1:d(bid+1),bid+1) = g(bid+1);
% end
% 
% % Pitch shift each block
% for m = 1:size(G,2)
%     for n = 1:size(G,1)
%         if isempty(G{n,m}) == 0
%             G{n,m}.rt = G{n,m}.rt + 2 * pi * frac / d(m) * G{n,m}.r;
%         end
%     end
% end




%% Convert to unsteady

% Write hdf5 file
if exist([job.directory job.name '_input.hdf5'],'file') ~= 0
    delete([job.directory job.name '_input.hdf5']);
end
ts_write_hdf5(g,[job.directory job.name '_input.hdf5']);

% Open python conversion script
fid = fopen([job.directory 'convert_unsteady.py'],'w');

% Write headers in python script
fprintf(fid,'%s\n','from ts import ts_tstream_reader, ts_tstream_type, ts_tstream_load_balance');
fprintf(fid,'%s\n\n','from ts import ts_tstream_patch_kind, ts_tstream_steady_to_unsteady');
fprintf(fid,'%s\n\n','tsr = ts_tstream_reader.TstreamReader()');

% Read in hdf5
fprintf(fid,'%s\n\n',['g = tsr.read("' job.hpcdir job.name '_input.hdf5")']);

% Duplicate and scale data
fprintf(fid,'n_duplicate = [%i',d(1));
for n = 2:length(d); fprintf(fid,', %i',d(n)); end;
fprintf(fid,']\nn_scale = [%i',s(1));
for n = 2:length(s); fprintf(fid,', %i',s(n)); end;

% Periodic data
fprintf(fid,']\npid_per = [(%i,%i)',pps(1,1:2));
for n = 2:size(pps,1); fprintf(fid,', (%i,%i)',pps(n,1:2)); end;

% Call conversion function
fprintf(fid,']\n%s\n\n',['g = ts_tstream_steady_to_unsteady.steady_to_unsteady'...
    '(g, n_duplicate, n_scale, pid_per, find_periodic=True)']);

% Load balance and write out
% fprintf(fid,'ts_tstream_load_balance.load_balance(g, %i)\n\n',job.cards * job.ppn);
fprintf(fid,'%s\n',['g.write_xdmf("' job.hpcdir job.name '_input.xdmf", "x", "r", "rt")']);
fprintf(fid,'%s\n',['g.write_hdf5("' job.hpcdir job.name '_input.hdf5")']);
fclose(fid);

% Run python script locally on darwin or remotely
if strcmp(job.host, 'local') == 1
    
    % Run python
    [~,~] = system(['python' job.directory 'convert_unsteady.py']);  
    
else
        
    % Transfer python file
    [~,~] = system(['scp ' job.directory 'convert_unsteady.py ' job.user '@' job.host ':' ...
        job.hpcdir '.']);
    
    % Transfer input file
    [~,~] = system(['scp ' job.directory job.name '_input.hdf5 ' job.user '@' job.host ':'...
        job.hpcdir '.']);    
    
    % Run python
    if strcmp(job.arch,'pascal') == 1
        [~,~] = system(['ssh ' job.user '@' job.host ...
            ' "source /rds/project/hpc/rds-hpc-pullan/ts3/bashrc_module_ts361 && python '...
            job.hpcdir 'convert_unsteady.py"']);
    else
        [~,~] = system(['ssh ' job.user '@' job.host ' "module load turbostream2.5 && python '...
            job.hpcdir 'convert_unsteady.py"']);
    end
    
    % Copy input file back
    [~,~] = system(['scp ' job.user '@' job.host ':' job.hpcdir job.name '_input.hdf5 ' ...
        job.directory '.']);
    [~,~] = system(['scp ' job.user '@' job.host ':' job.hpcdir job.name '_input.xdmf ' ...
        job.directory '.']);
    
end

end

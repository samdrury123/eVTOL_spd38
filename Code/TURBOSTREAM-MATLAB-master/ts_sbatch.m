function job = ts_sbatch(job,g,plot_stuff)
% TS_SBATCH  Submit a TURBOSTREAM job to the Wilkes cluster
%
%   job = TS_SBATCH(job,g)
%
%   job - a data structure with the following components:
%       host - 'local' or Wilkes node name for MATLAB running on Wilkes itself or on a remote machine
%       user - string of username on Wilkes
%       nodes - number of nodes to use
%       walltime - number of hours for the job
%       priority - string of account name to use
%       directory - string of local directory
%       hpcdir - string of Wilkes directory if job is 'remote'
%       name - string part of filename
%   g - optional argument containing TURBOSTREAM cell array write the input file before submission
%   plot_stuff - optional 0 or 1 for showing run name

% Write TS input file
if exist('g','var') ~= 0
    if length(g) > 1; g = ts_load_balance(g,job.cards); end;
    ts_write_hdf5(g,[job.directory job.name '_input.hdf5']);
    ts_write_xdmf(g,[job.directory job.name '_input.xdmf']);
end

% Default to show run name
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

% Process all filenames
job.inname = [job.name '_input.hdf5'];
job.outname = [job.name '.hdf5'];
job.outxdmf = [job.name '.xdmf'];
job.avgname = [job.name '_avg.hdf5'];
job.avgxdmf = [job.name '_avg.xdmf'];
job.logname = ['log_' job.name '.txt'];
job.slurmname = ['slurm_' job.name];

% Preprocess the walltime
hours = floor(job.walltime);
minutes = round((job.walltime-hours)*60);

hours = int2str(hours);
if length(hours) < 2
    hours = ['0' hours];
end

minutes = int2str(minutes);
if length(minutes) < 2
    minutes = ['0' minutes];
end

% Calculate number of nodes required for the job
job.nodes = ceil(job.cards / job.ppn);

% Create the submission script by modifying the template
f = fileread('TS/slurm_submit_template');

% Replace all parameters in the file
f = strrep(f,'-J NAME',['-J ' job.name]);
f = strrep(f,'-p ARCH',['-p ' job.arch]);
f = strrep(f,'-A PRIORITY',['-A ' job.priority]);
f = strrep(f, 'nodes=N', ['nodes=' int2str(job.nodes)]);
f = strrep(f, 'ntasks=T', ['ntasks=' int2str(job.cards)]);
f = strrep(f, 'ppn=P', ['ppn=' int2str(min(job.ppn,job.cards))]);
f = strrep(f, 'mem=M', ['mem=' int2str(job.mem * job.nodes) 'mb']);
f = strrep(f, 'time=W', ['time=' hours ':' minutes ':00']);
f = strrep(f, 'INPUT', job.inname);
f = strrep(f, 'OUTPUT', job.name);
if strcmp(job.host, 'local') == 1
    f = strrep(f, 'workdir="DIRECTORY"', ['workdir="' job.directory '"']);
else
    f = strrep(f, 'workdir="DIRECTORY"', ['workdir="' job.hpcdir '"']);
end
f = strrep(f, 'LOG.txt', job.logname);

% Exclusive node use or not
if strcmp(job.arch,'pascal') == 1
    f = strrep(f, 'gpu:G', ['gpu:' int2str(min(job.ppn,job.cards))]);
else
    f = strrep(f,'#SBATCH --gres','#!SBATCH --gres');
end

% Write path to application based upon which cluster
if strcmp(job.arch,'pascal') == 1
    f = strrep(f, 'application=A', 'application="turbostream"');
elseif strcmp(job.arch,'tesla') == 1
    f = strrep(f, 'application=A',['application="/scratch/turbostream/turbostream_' job.version ...
        '/turbostream/bin/cuda/turbostream"']);
end

% Write execution command based upon which cluster
if strcmp(job.arch,'pascal') == 1
    f = strrep(f, 'CMD=C', 'CMD="mpirun -np $np $application $options"');
elseif strcmp(job.arch,'tesla') == 1
    f = strrep(f, 'CMD=C', ['CMD="mpirun -ppn $mpi_tasks_per_node -np $np gpu_run_task.sh ' ...
        '$application $options"']);
end

% Write module load commands based upon which cluster
if strcmp(job.arch,'pascal') == 1
    f = strrep(f, 'MODULES',['source /rds/project/hpc/rds-hpc-pullan/ts3/ts' ...
        strrep(job.version,'.','') '/bashrc_module_ts' strrep(job.version,'.','') ]);
elseif strcmp(job.arch,'tesla') == 1
    f = strrep(f, 'MODULES',sprintf(['module load default-wilkes\nmodule load turbostream'...
        job.version]));
end

% Write submission file
fid = fopen([job.directory job.slurmname],'w');
fprintf(fid,'%s',f);
fclose(fid);

% Delete log and output files
[~,~] = system(['rm ' job.directory job.logname]);
if strcmp(job.host,'local') == 0
    [~,~] = system(['ssh ' job.user '@' job.host ' rm ' job.hpcdir job.logname]);
    [~,~] = system(['ssh ' job.user '@' job.host ' rm ' job.hpcdir job.outname]);
end

% Submit job to the queue or run on interactive node
nn_try = 10; n_try = 0; status = 1;
while status ~= 0 && n_try < nn_try
    if strcmp(job.host, 'local') == 1

        % Submit job and delete submission file
        [status,result] = system(['sbatch ' job.directory job.slurmname]);

    elseif isfield(job,'node') == 0

        % Transfer input file
        [~,~] = system(['scp ' job.directory job.inname ' ' job.user '@' job.host ':' job.hpcdir '.']);

        % Transfer submission file
        [~,~] = system(['scp ' job.directory job.slurmname ' ' job.user '@' job.host ':' ...
            job.hpcdir '.']);

        % Submit job
        [status,result] = system(['ssh ' job.user '@' job.host ' sbatch ' job.hpcdir job.slurmname]);
        
    else
        
        % Transfer input file
        [~,~] = system(['scp ' job.directory job.inname ' ' job.user '@' job.host ':' job.hpcdir '.']);
        
        % Generate bash script
        job.slurmname = ['submit_' job.name];
        fid = fopen([job.directory job.slurmname],'w');
        fprintf(fid,'%s\n','#!/bin/bash');
        fprintf(fid,'cd %s\n',job.hpcdir);
        fprintf(fid,'%s\n','module purge');
        fprintf(fid,'source /rds/project/hpc/rds-hpc-pullan/ts3/ts%s/bashrc_module_ts%s\n',...
            strrep(job.version,'.',''),strrep(job.version,'.',''));
        fprintf(fid,'mpirun -np %d turbostream %s %s > %s\n',...
            job.cards,job.inname,job.name,job.logname);
        fclose(fid);
        
        % Transfer bash script
        [~,~] = system(['scp ' job.directory job.slurmname ' ' job.user '@' job.host ':' ...
            job.hpcdir '.']);
        
        % Print job started
        if plot_stuff == 1; fprintf(['%% Executing TURBOSTREAM on ' job.name ' %%\n']); end;
        
        % Execute the job 
        [~,~] = system(['ssh -t ' job.user '@' job.host ' ssh ' job.node ' chmod u+x ' ...
            job.hpcdir job.slurmname]);
        [status,result] = system(['ssh -t ' job.user '@' job.host ' ssh ' job.node ' ' ...
            job.hpcdir job.slurmname]);

    end
    
    % Increment try counter
    if status ~= 0; pause(60); n_try = n_try + 1; disp([job.name ' Try ' num2str(n_try)]); end;
end

% Display error message if job submission fails
if status ~= 0
    disp(result)
end

% Set up for job monitoring
if isfield(job,'node') == 0
    
    % Find job id number
    job.id = regexp(result,'job\s+(\d+)','tokens'); job.id = str2double(job.id{1});

    % Initialise job tracking flags
    job.started = 0;
    job.finished = 0;
    
    % Print job started
    if plot_stuff == 1; fprintf(['%% Submitted TURBOSTREAM on ' job.name ' %%\n']); end;
    
else
    
    % Job has finished in interactive mode if command stops hanging
    job.started = 1; job.finished = 1;
    
end

% Delete submission files
[~,~] = system(['rm ' job.directory job.slurmname]);
if strcmp(job.host, 'local') == 0 
    [~,~] = system(['ssh ' job.user '@' job.host ' rm ' job.hpcdir job.slurmname]);
end

end
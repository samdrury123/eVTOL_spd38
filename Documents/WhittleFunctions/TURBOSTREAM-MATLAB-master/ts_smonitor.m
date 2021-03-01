function [job,finished] = ts_smonitor(job,monitor_external,reversed_flow)
% TS_SMONITOR  Monitor a TURBOSTREAM job for completion or failure and
% collect output files
%
%   job = TS_SMONITOR(job,monitor_external,reversed_flow)
%
%   job - structure defining job variables or cell array of multiple jobs

% Default to monitor inside function
if exist('monitor_external','var') == 0
    monitor_external = 0;
end

% Default to cancel jobs with negative flow
if exist('reversed_flow','var') == 0
    reversed_flow = 0;
end

% Assemble cell array if only one job is present
if iscell(job) == 0
    Jobs{1} = job;
else
    Jobs = job;
end

% Loop over all jobs in cell array
finished = zeros(size(Jobs));
while numel(Jobs) > 0    
    for j = length(Jobs):-1:1
        
        % Get current job
        job = Jobs{j};        
        
        % Check if job has been correctly submitted by the right fields in the struct
        if isfield(job,'outname') == 0
            job.started = 0;
            job.finished = 0;
            continue
        end
          
        % Query the queue to see if the job has finished
        if isfield(job,'id') == 1
            if strcmp(job.host, 'local') == 1
                [~,result] = system(['sacct -j' num2str(job.id)]);
            else
                [~,result] = system(['timeout 10s ssh ' job.user '@' job.host ' sacct -j ' ...
                    num2str(job.id)]);
            end
        end
        
        % Check if a job has failed for no good reason and resubmit
%         if (isempty(strfind(result,'FAILED')) == 0 || isempty(strfind(result,'TIMEOUT')) == 0) ...
%                 && isnan(job.finished) == 0
%             job = ts_sbatch(job); Jobs{j} = job;
%             continue
%         end
        
        % Try to copy the logfile across
        if strcmp(job.host, 'local') == 0
            [~,~] = system(['scp ' job.user '@' job.host ':' job.hpcdir job.logname ' ' ...
                job.directory '.']);
        end
        
        % Check if the log file exists
        if exist([job.directory job.logname],'file') == 0
            
            % Continue monitoring if it hasn't started
            job.started = 0;
            job.finished = 0;
            Jobs{j} = job;
            continue
        else
            
            % Raise started flag
            job.started = 1;
            
        end
        
        % Get first and last 20,000 lines of logfile
        [~,log_start] = system(['head -n 20000 ' job.directory job.logname]);
        log_start = strsplit(log_start,'\n').';
        [~,log_end] = system(['tail -n 20000 ' job.directory job.logname]);
        log_end = strsplit(log_end,'\n').';    
        log = [log_start ; log_end];
            
        % Monitor whole logfile for nans
%         nan_1 = regexp(log,'NAN');
%         nan_1 = find(cellfun('isempty',nan_1)==0, 1);
%         nan_2 = regexp(log,'nan');
%         nan_2 = find(cellfun('isempty',nan_2)==0, 1);
%         if ~isempty(nan_1) || ~isempty(nan_2)
%             job.finished = NaN;
%         end
        
        % Only monitor residuals for NaNs
        nan_1 = regexp(log,'TOTAL DAVG   NAN');
        nan_1 = find(cellfun('isempty',nan_1)==0, 1);    
        if ~isempty(nan_1)
            job.finished = NaN;
        end       
        
        % Monitor for negative volumes
        neg = regexp(log,'FOUND NEGATIVE VOLUME!!!');
        neg = find(cellfun('isempty',neg)==0, 1);
        if ~isempty(neg)
            job.finished = NaN;
        end
        
        % Monitor for reversed flow
        if reversed_flow == 0
            m_in = regexp(log,'INLET FLOW =\s+([\w\.\-]+)','tokens');
            m_in(cellfun(@isempty,m_in) == 1) = [];
            if isempty(m_in) == 0
                m_in = cellfun(@str2double,[m_in{:}]);
                if m_in(end) < 0 
                    job.finished = NaN;
                end
            end
        end
        
        % Monitor for termination
        term = regexp(log,'APPLICATION TERMINATED');
        term = find(cellfun('isempty',term)==0, 1);
        if ~isempty(term)
            job.finished = NaN;
        end
        
        % Kill a NaNed job
        if isnan(job.finished) && isfield(job,'id') == 1
            if strcmp(job.host, 'local') == 1
                [~,~] = system(['scancel ' num2str(job.id)]);
            else
                [~,~] = system(['ssh ' job.user '@' job.host ' scancel ' num2str(job.id)]);
            end
            disp(['Job ' num2str(job.id) ' ' job.name ' Cancelled']);
            Jobs(j) = []; finished(j) = nan;
            if isempty(Jobs) == 0
                continue
            else
                return
            end
        end
        
        % Check if the output file exists
        if strcmp(job.host,'local') == 1
            [status,~] = system(['ls ' job.directory job.outname]);
        else
            [status,~] = system(['ssh ' job.user '@' job.host ' ls ' job.hpcdir job.outname]);
        end
        
        % Assign finished status if queue is not running and output file exists
        if job.finished == 0 && isempty(strfind(result,'COMPLETED')) == 0 && status == 0
            job.finished = 1;
        end
        
        % If a job has finished collect output files if necessary
        if job.finished == 1 && strcmp(job.host,'local') == 0
            
            % Collect outfiles
            filenames = {job.logname job.outname job.outxdmf job.avgname job.avgxdmf};
            for f = 1:length(filenames)
                [~,~] = system(['scp ' job.user '@' job.host ':' job.hpcdir filenames{f} ' '...
                    job.directory '.']);
            end
            
            % Delete input file
            [~,~] = system(['rm ' job.directory job.inname]);
            [~,~] = system(['ssh ' job.user '@' job.host ' rm ' job.hpcdir job.inname]);
            [~,~] = system(['rm ' job.directory strrep(job.inname,'.hdf5','.xdmf')]);
            [~,~] = system(['ssh ' job.user '@' job.host ' rm ' job.hpcdir ...
                strrep(job.inname,'.hdf5','.xdmf')]);
            
            % Delete slurm out file
            if strcmp(job.host,'local') == 0 && isfield(job,'node') == 0
                [~,~] = system(['ssh ' job.user '@' job.host ' rm ~/slurm-' num2str(job.id) '.out']);
                [~,~] = system(['ssh ' job.user '@' job.host ' rm ' job.hpcdir...
                    'machine.file.' num2str(job.id)]);
            end
            
            Jobs(j) = [];
            
        elseif job.finished == 1
            
            % Delete input file
            [~,~] = system(['rm ' job.directory job.inname]);
            
            % Delete slurm out file
            [~,~] = system(['rm ./slurm-' num2str(job.id) '.out']);
            [~,~] = system(['rm ' job.directory 'machine.file.' num2str(job.id)]);            
            
            Jobs(j) = [];
            
        end
        
        % Check the output file is valid and resubmit if not
        if job.finished == 1 && (H5F.is_hdf5([job.directory job.outname]) == 0 || ...
                (exist([job.directory job.avgname],'file') ~= 0 && ...
                H5F.is_hdf5([job.directory job.avgname]) == 0))
            disp(['%% No valid output file ' job.name ' %%']);
            job = ts_sbatch(job); Jobs{j} = job;
        end
        
        % Record finished array
        finished(j) = job.finished;

    end
    
    % Return once all jobs have been monitored
    if monitor_external == 1
        job = Jobs;        
        return
    end
    
    % Error if nan found in interactive job
    if isfield(job,'id') == 0 && isnan(job.finished) == 1
        error(['NaN found in interactive job logfile ' job.name]);
    end
    
    % Add in delay
    if isfield(job,'node') == 0; pause(20); end;
    
end


end
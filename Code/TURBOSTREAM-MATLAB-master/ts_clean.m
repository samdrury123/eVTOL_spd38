function [] = ts_clean(directory)
% TS_CLEAN  Reduce disk space by processing TS files in a directory

% Loop over all files
A = dir(directory);
for n = 1:length(A)
    if A(n).isdir == 0 && length(A(n).name) > 4

        % Display progress
        disp(['File ' num2str(n) ' out of ' num2str(length(A)) ' ' A(n).name]);
        
        % Process log files
        if strcmp(A(n).name(1:3),'log') == 1 && strcmp(A(n).name(end-2:end),'txt') == 1
            try; ts_plot_conv([directory A(n).name],0,0,0); catch; delete([directory A(n).name]); end;
        end

        % Delete non-average hdf5 files and process average files
        if strcmp(A(n).name(end-3:end),'hdf5') == 1 
            if isempty(strfind(A(n).name,'avg')) == 1 && isempty(strfind(A(n).name,'unst')) == 1 && ...
                    isempty(strfind(A(n).name,'froz')) == 1
                delete([directory A(n).name])
            elseif isempty(strfind(A(n).name,'probe')) == 1
                ts_read_hdf5([directory A(n).name]);
            end
        end

        % Delete xdmf files
        if strcmp(A(n).name(end-3:end),'xdmf') == 1 
            delete([directory A(n).name])        
        end
    end
    
end


end


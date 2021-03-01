function [] = move_files(dr,search_str,start_str,end_str)
% MOVE_FILES  Search a directory struct for matching files, then rename using strrep

% Define directory varnames to loop over
varnames = fieldnames(dr);

% Loop over all subdirectories
for v = 1:length(varnames)
    
    % Find matching files
    A = dir([dr.(varnames{v}) search_str]);
    
    % Loop over all files, rename and move
    for a = 1:length(A)
        movefile([dr.(varnames{v}) A(a).name],[dr.(varnames{v}) strrep(A(a).name,start_str,end_str)]);
    end
end

end
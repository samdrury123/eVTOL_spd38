function [F,S] = ts_get_files(filestring,use_avg)
% TS_GET_FILES  Search a directory for full TURBOSTREAM solutions
%
%   [F,S] = TS_GET_FILES(filestring,use_avg)
%
%   filestring - input part of filename including directory but not extension
%   use_avg - optional argument, 0 or 1 to select file type
%   F - cell array of filenames
%   S - struct containing fields describing the file contents

% Default to use average if not specified
if exist('use_avg','var') == 0
    use_avg = 1;
end

% Search directory for filenames and allocate in cell array
A = [dir([filestring '*.hdf5']) ; dir([filestring '*.mat'])];
F = cell(length(A),1); for n = 1:length(A); F{n} = A(n).name; end;

% Remove all input files
F(cellfun(@isempty,strfind(F,'input')) == 0) = [];

% Remove all cut files
F(cellfun(@isempty,strfind(F,'cut')) == 0) = [];

% Remove all chic files
F(cellfun(@isempty,strfind(F,'chic')) == 0) = [];

% Remove all paraview files
F(cellfun(@isempty,strfind(F,'highspeed')) == 0) = [];

% Remove all probe files
F(cellfun(@isempty,strfind(F,'probe')) == 0) = [];

% Remove all lost efficiency files
F(cellfun(@isempty,strfind(F,'leta')) == 0) = [];

% Remove all unsteady and frozen rotor files
if isempty(strfind(filestring,'froz')) == 1; F(cellfun(@isempty,strfind(F,'froz')) == 0) = []; end;
if isempty(strfind(filestring,'unst')) == 1; F(cellfun(@isempty,strfind(F,'unst')) == 0) = []; end;

% Convert all .mat files to .hdf5 files
F = strrep(F,'.mat','.hdf5');

% Select only unique filenames
F = unique(F);

% Remove normal names if using averages or remove avg names if not using averages
if use_avg == 1
    F(cellfun(@isempty,strfind(F,'avg')) == 1) = [];
elseif use_avg == 0
    F(cellfun(@isempty,strfind(F,'avg')) == 0) = [];
end

% Split filenames based upon underscores
F_part = regexp(F,'_','split');

% Loop over all filenames and extract data
S = [];
for i = 1:length(F_part)
    
    % Extract numeric and non-numeric subparts from the name
    B = regexp(F_part{i},'[A-Z]+','match');
    C = regexp(F_part{i},'-?\d+\.?\d*|-?\d*\.?\d+','match');
    
    % Check which sub parts have both
    q = cellfun(@isempty,B) == 0 & cellfun(@isempty,C) == 0;
    B = B(q); C = C(q);
    
    % Record in a struct the numeric values
    for j = 1:length(B)
        S.(B{j}{1})(i) = str2double(C{j}{1});
    end
end

end
% Set default figure position
% set(0, 'DefaultFigurePosition', [471 342 960 540])

% Load all folders and subfolders on the default path
Paths = genpath('C:\Users\Sam Drury\Documents\CAMBRIDGE\IIB\eVTOL\Code');
Paths = regexp(Paths,';','split');
Paths(cellfun('isempty',regexp(Paths,'Backup','match')) == 0) = [];
Paths(cellfun('isempty',regexp(Paths,'Figures','match')) == 0) = [];
Paths(cellfun('isempty',regexp(Paths,'Export','match')) == 0) = [];
Path = [];
for n = 1:length(Paths)
    if isempty(strfind(Paths{n},'.git')) == 1
        Path = [Path Paths{n} ';'];
    end
end
addpath(Path);

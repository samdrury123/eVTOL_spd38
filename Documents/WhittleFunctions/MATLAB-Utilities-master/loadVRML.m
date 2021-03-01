function wrl = loadVRML(varargin)
%LOADVRML
%
%Format
%
%   "wrl = loadVRML;"
%
%   "wrl = loadVRML(filename);
%

%program defaults - will be overwritten by the contents of the defaults
%file if it exists.
startPath = getuserdir;

%manages variable input argumetns
fname = '';
if nargin >= 1
    fname = varargin{1};
end

%housekeeping
currentdir = cd;

%gets filename from user if not explicitily specified
if isempty(fname)
    %loads defaults file if in command line mode
    defaultsFile = 'loadVRMLdefaults.mat';
    defaultsFile = fullfile(getuserdir,'MATLAB',defaultsFile);
    if exist(defaultsFile,'file')
        load(defaultsFile);
    end
    
    %gets file from user
    cd(startPath);
    filterspec = {'*.wrl','VRML files (*.wrl)';
        '*.*','All files (*.*)'};
    titstr = 'Open VRML file:';
    [fname,pname] = uigetfile(filterspec,titstr);
    cd(currentdir);
    
    %checks if user pressed canecl
    if fname == 0
        wrl = -1;
        return
    end
    
    %generates full file path
    fname = fullfile(pname,fname);
    
    %saves defaults file ready for next time
    startPath = pname;
    save(defaultsFile,'startPath');
    
else 
    if exist(fname,'file')
       errstr = ['File "',fname,'" not found.'];
       error(errstr);
    end
end

%loads VRML file
[nel,tmp,infoline] = read_vrml(fname);
wrl = tmp(1:nel);

%ends function loadVRML
end
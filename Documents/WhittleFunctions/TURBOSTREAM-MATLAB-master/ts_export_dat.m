function [] = ts_export_dat(g,filename)
% Write traverse data to a text file

% Choose variables to write
varnames = {'x','y','z','Po','Alpha','Beta','P','Vx','Vr','Vt'};

% Open file to write
fid = fopen(filename,'w');

% Write heading of grid size
[nj,nk] = size(g.(varnames{1}));
fprintf(fid,'%g\t%g\n',[nj nk]);

% Loop over all variables
for v = 1:length(varnames)
    
    % Write variable name
    fprintf(fid,'%s\n',varnames{v});
    
    % Write data in a tab delimited array
    fprintf(fid,[repmat('%13.8g\t',[1 nj]) '\n'],g.(varnames{v}));

    % Write a blank line after the array
    fprintf(fid,'\n');
    
end

% Close the file
fclose(fid);

end
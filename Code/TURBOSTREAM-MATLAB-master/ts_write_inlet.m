function A = ts_write_inlet(filename,c)
% TS_WRITE_INLET  Save a cut plane in a text file for use as an inlet boundary condition

% Define variables to write
varnames = {'r' 'To' 'Po' 'Alpha' 'Beta'};

% Assemble inlet boundary condition averages
A = [];
for v = 1:length(varnames)
    A = [A ts_mass_average(c,varnames{v},2)'];
end

% Pad data with zeros
A = [A zeros(size(A,1),2)];

% Write data to file
dlmwrite(filename,A,'roffset',1,'delimiter','\t','precision','%10.8g'); 

end
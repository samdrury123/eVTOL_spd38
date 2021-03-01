function xyz = fr_read_sections(filename)
% FR_READ_SECTIONS  Read a text file of sections taken with faro scanning arm
%
%   xyz = FR_READ_SECTIONS(filename)
%
%   filename - string of text file to read in
%   xyz - output cell array of coordinate points

% Read in all text
s = fileread(filename);

% Find headers of each section
[i1,i2] = regexp(s,'\#[\w\-\s\(\)]+\n','start','end');

% Add end of file onto i1
i1 = [i1 length(s)+1];

% Pre allocate cell array
xyz = cell(length(i2),1);

% Loop over all sections
for n = 1:length(i2)
   
    % Get coordinate points
    xyz_temp = textscan(s(i2(n)+1:i1(n+1)-1),'%f');
    
    % Reshape into array
    xyz{n} = reshape(xyz_temp{1},3,[])';
    
    % Convert to metres
    xyz{n} = xyz{n}/1000;
end

end


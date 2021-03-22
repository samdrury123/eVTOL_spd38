function xyz = bl_read_ibl(filename)
% BL_READ_IBL  Read an ibl file to get back blade definition data

% Read file
f = fileread(filename);

% Read coordinates
curv = regexp(f,'\s+([\w\.\-\+]+)\s+([\w\.\-\+]+)\s+([\w\.\-\+]+)\n','tokens');
xyz_curv = cellfun(@str2double,cat(1,curv{:}));

% Split coordinates into sections
i = find(isnan(xyz_curv(:,2)) == 1); i = [i ; size(xyz_curv,1)+1];

% Pre-allocate section cell array
nj = length(i)-1; xyz = cell(nj,1);

% Loop over each section and save coordinates in metres
for j = 1:nj
    xyz{j} = xyz_curv(i(j)+1:i(j+1)-1,:) / 1000;
end


end


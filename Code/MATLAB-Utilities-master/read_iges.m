function xyz = read_iges(filename)
% Read a GOM scanned IGES file for blade coordinates

fid = fopen(filename,'rt');

marker=0;
file_end=0;
line = fgetl(fid);

while file_end==0
    store = sscanf(line(1:3),'%f');
    
    % If 106 is found then skip line increment marker
    if store==106; 
        line = fgetl(fid);
        marker=marker+1;
        j=0;
    end
    if marker>0
        j=j+1;
        loc=find(line==',' | line==';');
		xyz{marker}(j,:) = [sscanf(line(1:loc(1)-1),'%f') sscanf(line(loc(1)+1:loc(2)-1),'%f')...
            sscanf(line(loc(2)+1:loc(3)-1),'%f')];
    end
    
    % Get the next line before starting the loop again so that it can check for the end of the file
    line = fgetl(fid); 
    
    % File ends with a 'S' so this line looks for it
    if line(1)=='S' 
        file_end = 1;
    end 
end

end
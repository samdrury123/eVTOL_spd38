function [] = write_ibl(filename,xrrt,flip_z,t_shift,nj)
% Write and ibl file with sections for import into pro engineer in mm

if exist('nj','var') == 0
    nj = 31;
end

fid = fopen(filename,'w');

% Maximum number of points
ni = 500;
I = round(linspace(1,size(xrrt,1),ni));
J = round(linspace(1,size(xrrt,2),nj));

% Shift z to specified value
t = xrrt(:,:,3) ./ xrrt(:,:,2);
t = t + t_shift;
x = xrrt(:,:,1); r = xrrt(:,:,2);
[y, z] = pol2cart(t, r);

% Flip if necessary
if flip_z == 1
    z = -z;
end

xyz = cat(3,x,y,z);

fprintf(fid,'%s\n','Closed Index arclength');
for j = J
    fprintf(fid,'%s\n','begin section');
    fprintf(fid,'%s\n','begin curve');
    
    for i = I
        fprintf(fid,'%10.11f %10.11f %10.11f\n',squeeze(xyz(i,j,:)) * 1000);
    end
end

% Leading edge line
fprintf(fid,'%s\n','begin section');
fprintf(fid,'%s\n','begin curve');
for j = 1:size(xyz,2)
    fprintf(fid,'%10.11f %10.11f %10.11f\n',squeeze(xyz(1,j,:)) * 1000);
end

% Trailing edge line
ite = round(size(xyz,1)/2);
fprintf(fid,'%s\n','begin section');
fprintf(fid,'%s\n','begin curve');
for j = 1:size(xyz,2)
    fprintf(fid,'%10.11f %10.11f %10.11f\n',squeeze(xyz(ite,j,:)) * 1000);
end

fclose(fid);

end
function [blade] = Blades(blade, p)
%% Return complete blade structure
% [blade] = Blades(blade)
NUMBEROFSECTIONS = 21;

%% Get blade points from KON_Make_Blade
[blade.s_steps,blade.p_steps,blade.sections,blade.pointspersection,blade.XYZ1,blade.XYZ2,clyoffset] = KON_Make_Blade(blade.sec.chi1,blade.sec.chi2,blade.sec.radius,blade.sec.chord,blade.type,NUMBEROFSECTIONS);

%% Offset blades (3D and zoffset)
blade.zoffset = blade.z - min(min(squeeze(blade.s_steps(:,3,:))));

% Get offset distributions due to sweep and lean
[offset] = SweepLean(blade,p);

blade.zoffset = offset.z;
blade.yoffset = offset.y;
clear offset;

for stepi = 1:blade.sections
    blade.s_steps(:,3,stepi)  = blade.zoffset(stepi) + blade.s_steps(:,3,stepi);
    blade.p_steps(:,3,stepi)  = blade.zoffset(stepi) + blade.p_steps(:,3,stepi);

    blade.s_steps(:,2,stepi)  = blade.s_steps(:,2,stepi) + blade.yoffset(stepi);
    blade.p_steps(:,2,stepi)  = blade.p_steps(:,2,stepi) + blade.yoffset(stepi);

    blade.XYZ1(:,3,stepi)    = blade.zoffset(stepi) + blade.XYZ1(:,3,stepi);
    blade.XYZ2(:,3,stepi)    = blade.zoffset(stepi) + blade.XYZ2(:,3,stepi);

    blade.XYZ1(:,2,stepi)    = blade.XYZ1(:,2,stepi) + blade.yoffset(stepi);
    blade.XYZ2(:,2,stepi)    = blade.XYZ2(:,2,stepi) + blade.yoffset(stepi);
end

%% IBL
% Concatenate PS & SS

for i = 1:size(blade.XYZ1,3)
    blade.XYZ(:,:,i) = [blade.XYZ1(:,:,i);blade.XYZ2(:,:,i)] * 1000; % mm
end
clear blade.XYZ1 clear blade.XYZ2;

% Get centerline values
clzoffset = blade.zoffset;
clradius = linspace(blade.sec.radius(1), blade.sec.radius(end), blade.sections);
blade.centerline = [clradius.',clyoffset.',clzoffset.'] .* 1000; % mm

%% Write .ibl in X,Y,Z (STATOR)
pathname    = [pwd, '/IBL/'];
iblfile = [pathname blade.type '.ibl'];
statorcenterlinefile = fullfile(pathname,'CenterLine.ibl');

%SAM
figure(5);
%subplot(4,6,[5:6,11:12]);
hold on; box on
%SAM
fid = fopen(iblfile, 'w');
fprintf(fid,'Closed Index arclength\n');
for ss = 1:size(blade.XYZ,3)
     fprintf(fid,'begin section\n');
     fprintf(fid,'begin curve\n');
     for pp = 1:size(blade.XYZ,1)
         fprintf(fid,'%f  %f  %f\n', blade.XYZ(pp,1,ss),blade.XYZ(pp,2,ss),blade.XYZ(pp,3,ss));
     end
     %SAM
     plot3(blade.XYZ(:,2,ss),-blade.XYZ(:,3,ss),blade.XYZ(:,1,ss),'-k') ;
     axis equal
     %SAM
end

%SAM
axis equal
axis([-5 15 -60 0 0 60]); view(3);
title('3D Plot - Rotate Me')
%SAM

fclose(fid);

fid = fopen(statorcenterlinefile, 'w');
fprintf(fid,'Closed Index arclength\n');
fprintf(fid,'begin section\n');
fprintf(fid,'begin curve\n');
fprintf(fid,'%f  %f  %f\n', (blade.centerline(1,1)-5),blade.centerline(1,2),blade.centerline(1,3));
for ss = 1:size(blade.centerline,1)
     fprintf(fid,'%f  %f  %f\n', blade.centerline(ss,1),blade.centerline(ss,2),blade.centerline(ss,3));
end
fprintf(fid,'%f  %f  %f\n', (blade.centerline(end,1)+5),blade.centerline(end,2),blade.centerline(end,3));
fclose(fid);

end
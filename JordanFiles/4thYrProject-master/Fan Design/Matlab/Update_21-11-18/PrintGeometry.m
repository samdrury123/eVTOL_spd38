global flagplot

for n1=1:length(file_name)
    load([file_name{n1} '.mat']);
    chi1_rotor     = r.chi1;
    chi2_rotor     = r.chi2;
    R_rotor        = r.radius;
    c_rotor        = r.chord;
    N_rotor        = r.blades;
    chi1_stator    = s.chi2;
    chi2_stator    = s.chi3;
    R_stator       = s.radius;
    
    c_stator       = s.chord;
    N_stator       = s.blades;
    [r_s_steps,r_p_steps,r_sections,r_pointspersection,r_XYZ1,r_XYZ2]= KON_Make_Blade(chi1_rotor,chi2_rotor,R_rotor,c_rotor,'rotor',file_name{n1});    
    [s_s_steps,s_p_steps,s_sections,s_pointspersection,s_XYZ1,s_XYZ2]= KON_Make_Blade(chi1_stator,chi2_stator,R_stator,c_stator,'stator',file_name{n1});
     
    if flagplot == 1
    for i=1:21
     figure(3)
     hold on
%      plot3(r_XYZ1(:,1,i), r_XYZ1(:,2,i), r_XYZ1(:,3,i))
%      plot3(r_XYZ2(:,1,i), r_XYZ2(:,2,i), r_XYZ2(:,3,i))

     hold on
     plot3(s_XYZ1(:,1,i), s_XYZ1(:,2,i), 0.03+s_XYZ1(:,3,i))
     plot3(s_XYZ2(:,1,i), s_XYZ2(:,2,i), 0.03+s_XYZ2(:,3,i))
     
     
     xlabel('x')
     ylabel('y')
     zlabel('z')
     
    end
    end
   
     
    Z_rot          = min(min(squeeze(r_s_steps(:,3,:))));
    Z_stat         = min(min(squeeze(s_s_steps(:,3,:))));
    rotor_offset   = 0.0125 - Z_rot;
    stator_offset  = 0.037318- Z_stat;
    r_s_steps(:,3,:)  = rotor_offset + r_s_steps(:,3,:);
    r_p_steps(:,3,:)  = rotor_offset + r_p_steps(:,3,:);
    s_s_steps(:,3,:)  = stator_offset + s_s_steps(:,3,:);
    s_p_steps(:,3,:)  = stator_offset + s_p_steps(:,3,:);
    
    r_XYZ1(:,3,:)    = rotor_offset + r_XYZ1(:,3,:);
    r_XYZ2(:,3,:)    = rotor_offset + r_XYZ2(:,3,:);
    s_XYZ1(:,3,:)    = stator_offset + s_XYZ1(:,3,:);
    s_XYZ2(:,3,:)    = stator_offset + s_XYZ2(:,3,:) ;
       
    nose=0;
    
%% IBL
% Concatenate PS & SS

for i = 1:size(r_XYZ1,3)
    r_XYZ(:,:,i) = [r_XYZ1(:,:,i);r_XYZ2(:,:,i)] * 1000; % mm
    s_XYZ(:,:,i) = [s_XYZ1(:,:,i);s_XYZ2(:,:,i)] * 1000; % mm
end


% Write .ibl in X,Y,Z (STATOR)
pathname    = [pwd, '/IBL'];
statorfile  = fullfile(pathname, 'Stator.ibl');
rotorfile   = fullfile(pathname,'Rotor.ibl');

fid = fopen(statorfile, 'w');
fprintf(fid,'Closed Index arclength\n');
for ss = 1:size(s_XYZ,3)
     fprintf(fid,'begin section\n');
     fprintf(fid,'begin curve\n');
     for pp = 1:size(s_XYZ,1)
         fprintf(fid,'%f  %f  %f\n', s_XYZ(pp,1,ss),s_XYZ(pp,2,ss),s_XYZ(pp,3,ss));
     end
end

fclose(fid);

% Write .ibl in X,Y,Z (ROTOR)
fid = fopen(rotorfile, 'w');
 
fprintf(fid,'Closed Index arclength\n');
for ss = 1:size(r_XYZ,3)
     fprintf(fid,'begin section\n');
     fprintf(fid,'begin curve\n');
     for pp = 1:size(r_XYZ,1)
         fprintf(fid,'%f  %f  %f\n', r_XYZ(pp,1,ss),r_XYZ(pp,2,ss),r_XYZ(pp,3,ss));
     end
end
fclose(fid);
end
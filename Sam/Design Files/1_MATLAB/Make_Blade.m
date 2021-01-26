%%
function blade = Make_Blade(chi1,chi2,r,c,bladenum,fignum, d)
figure(fignum); 

%%SAM
pathname    = [pwd, './IBL/'];
statorcenterlinefile = fullfile(pathname,'CentreLine.ibl');
if bladenum == 1
    iblfile = [pathname 'rotor.ibl'];
elseif bladenum == 2
    iblfile = [pathname 'stator.ibl'];
end
fid = fopen(iblfile, 'w');
fprintf(fid,'Closed Index arclength\n');
%%SAM

%y needs to be -ve for CW, +ve for CCW
CW = true; %consider setting this as a global
if CW == true; rotation = -1; 
else rotation = 1; end

num_slices = size(chi1,1);
blade.centreline = zeros(num_slices,3);
for rr=1:num_slices
    if bladenum == 1 
        [x,rt,cl] = bladecentroid(chi1(rr),chi2(rr),c(rr),'rotor');
    elseif bladenum == 2
        [x,rt,cl] = bladecentroid(chi1(rr),chi2(rr),c(rr),'stator'); 
    end
    blade.centreline(rr,1) = cl(1); blade.centreline(rr,2) = cl(2);
    y = r(rr) .* sin(rt ./ r(rr));
    z = r(rr) .* cos(rt ./ r(rr));
    
    %This extends the blades to ensure CAD merges at hub and tip
    if rr==1; z=z-0.0001; end
    if rr==num_slices; z=z+0.00015; end
    
    if bladenum == 1; x=x-(num_slices-rr)^3 *1.5/1e7; end %This ensures the rotor TE is more uniform for thick chords
    if bladenum == 2; x=x+0.020; blade.centreline(rr,1)=blade.centreline(rr,1)+0.020; end
    
    % Adding symmetric parabolic lean (y) and parabolic sweep (x) to rotor and stator
    % in degrees, of form kz^2 where z is rr
    r_lean_angle = d.rotor_lean; %These are the angles at the tip and hub
    r_sweep_angle = d.rotor_sweep;
    s_lean_angle = d.stator_lean;
    s_sweep_angle = d.stator_sweep;
    
    % Quadratic variant
    r_lean = tan(r_lean_angle*pi/180)/2/(r(end) - mean(r));
    r_sweep = tan(r_sweep_angle*pi/180)/2/(r(end) - mean(r));
    % Quartic variant
    r_lean4 = tan(r_lean_angle*pi/180)/4/((r(end) - mean(r))^3);
    r_sweep4 = tan(r_sweep_angle*pi/180)/4/((r(end) - mean(r))^3);
    
    s_lean = tan(s_lean_angle*pi/180)/2/(r(end) - mean(r));
    s_sweep = tan(s_sweep_angle*pi/180)/2/(r(end) - mean(r));
    s_lean4 = tan(s_lean_angle*pi/180)/4/((r(end) - mean(r))^3);
    s_sweep4 = tan(s_sweep_angle*pi/180)/4/((r(end) - mean(r))^3);
    
    
    if bladenum == 1
        xshift_fact = 15; %ensures it fits within CAD
        blade.x(rr,:) = x - r_sweep*((z - mean(r)).^2) + r_sweep*((r(end) - mean(r)).^2)/2; %takes away sweep then recentres the blades
        %blade.x(rr,:) = x - r_sweep4*((z - mean(r)).^4) + r_sweep4*((r(end) - mean(r)).^4)/2; %q
        blade.rt(rr,:) = rt;
        blade.r(rr,:) = r(rr);
        blade.t(rr,:) = rt ./ r(rr);
        %blade.y(rr,:) = y - r_lean*((z - mean(r)).^2) + r_lean*((r(end) - mean(r)).^2)/2;
        blade.y(rr,:) = (y - r_lean4*((z - mean(r)).^4) + r_lean4*((r(end) - mean(r)).^4)).*rotation;
        blade.z(rr,:) = z;
        
    elseif bladenum == 2
        xshift_fact = 19; %ensures it fits within CAD
        blade.x(rr,:) = x - s_sweep*((z - mean(r)).^2) + s_sweep*((r(end) - mean(r)).^2)/2;
        blade.centreline(rr,1) = blade.centreline(rr,1) - s_sweep*((r(rr) - mean(r))^2) + s_sweep*((r(end) - mean(r))^2)/2;
        %blade.x(rr,:) = x - s_sweep4*((z - mean(r)).^4) + s_sweep4*((r(end) - mean(r)).^4)/2; %q
        blade.rt(rr,:) = rt;
        blade.r(rr,:) = r(rr);
        blade.t(rr,:) = rt ./ r(rr);
        %blade.y(rr,:) = y - s_lean*((z - mean(r)).^2) + s_lean*((r(end) - mean(r)).^2)/2;
        blade.y(rr,:) = (y + s_lean4*((z - mean(r)).^4) - s_lean4*((r(end) - mean(r)).^4)).*rotation; %note sign change
        blade.centreline(rr,2) = (blade.centreline(rr,2) + s_lean4*((r(rr) - mean(r))^4) - s_lean4*((r(end) - mean(r))^4))*rotation; %note sign change
        blade.z(rr,:) = z;
        blade.centreline(rr,3) = r(rr);
    end
    
    %%SAM
    fprintf(fid,'begin section\n');
    fprintf(fid,'begin curve\n');
    for pp = 1:size(x,2)
        fprintf(fid,'%f  %f  %f\n', blade.z(rr,pp)*1000,blade.y(rr,pp)*1000,blade.x(rr,pp)*1000 + xshift_fact); % +15 ensures the blade fits within CAD 
    end
end
fclose(fid);

%Writing line for hole in stator for wires
if bladenum == 2
    blade.centreline = blade.centreline .* 1000; %mm
    fid2 = fopen(statorcenterlinefile, 'w');
    fprintf(fid2,'Closed Index arclength\n');
    fprintf(fid2,'begin section\n');
    fprintf(fid2,'begin curve\n');
    fprintf(fid2,'%f  %f  %f\n', (blade.centreline(1,3)-4),blade.centreline(1,2),blade.centreline(1,1)+xshift_fact);
    for ss = 1:size(blade.centreline,1)
         fprintf(fid2,'%f  %f  %f\n', blade.centreline(ss,3),blade.centreline(ss,2),blade.centreline(ss,1)+xshift_fact);
    end
    fprintf(fid2,'%f  %f  %f\n', (blade.centreline(end,3)+4),blade.centreline(end,2),blade.centreline(end,1)+xshift_fact);
    fclose(fid2);
end
%%SAM
%if bladenum == 2; blade.x = blade.x + 0.0075; end
subplot(4,6,[5:6,11:12]);hold on; box on
plot3(-blade.y',-blade.x',blade.z','-k') ; axis equal
axis([-.02 .02 -.042 .01 0.01 0.07]); view(3);
%title('3D Plot - Rotate Me')

subplot(4,6,[3:4,9:10]);hold on; box on
plot(-blade.y',-blade.x','-k') ; axis equal
axis([-.02 .02 -.042 .01])
% subplot(4,6,[3:4,9:10]);hold on; box on
% plot(blade.y',blade.x','-k') ; axis equal
% axis([-.02 .02 -.005 .04])
%   
% subplot(4,6,[5:6,11:12]);hold on; box on
% plot3(blade.y',blade.x',blade.z','-k') ; axis equal
% axis([-.02 .02 -.005 .04 0.01 0.07]); view(3);
% title('3D Plot - Rotate Me')


end

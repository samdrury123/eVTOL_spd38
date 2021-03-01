function [offset] = SweepLean(blade, p)
%% Return offset vectors fro sweep and lean profiles
% [offset] = SweepLean(blade, parameter)

%% Sweep and Lean Offsets
% Set distribution and apply in y and z directions

tipangle = blade.sweep;
leangle = blade.lean;

% Tip offset ('winglet', y direction)
dtips = (blade.sec.radius(end)/5)*tan(2*pi*tipangle/360);
pv = polyfit([0 0.15 0.25 0.45 0.65 1],[0 dtips*(0.15/0.65) dtips*(0.25/0.65) dtips*(0.45/0.65) dtips 0],5);
dy_offset = 0.8*polyval(pv,linspace(0,1,blade.sections));

% Leading edge offset (z direction)
dle = 1.25*(blade.sec.radius(end)/5)*tan(2*pi*tipangle/360);
pv = polyfit([0 0.25 0.45 0.6 0.85 1],[0 dle*(0.15/0.65) dle*(0.25/0.65) dle*(0.45/0.65) dle 0],5);
dz_offset = polyval(pv,linspace(0,1,blade.sections));

switch p
    case 'constangle'
        % MID- FREE AND FORCE VORTEX
        dz_offset = polyval(pv,linspace(0,1,blade.sections));
    case 'custom'
        % MID- FREE AND FORCE VORTEX
        dz_offset = polyval(pv,linspace(0,1,blade.sections));
    case 'free'
        if strcmp(blade.type, 'rotor')
            % FREE VORTEX
            curvepoints = linspace(0,1,blade.sections);
            curveoffset = (6 - 22*exp(-curvepoints*2))*10^-3;
            dz_offset = curveoffset + polyval(pv,linspace(0,1,blade.sections));
        end
end

% Set sweep and lean offsets
% offset.z = blade.zoffset.*ones([1,blade.sections]) + dz_offset; 
% offset.y = -dy_offset;

if strcmp(blade.type, 'rotor')
    offset.z = blade.zoffset.*ones([1,blade.sections])% + dz_offset;% + linspace(-1e-3, 2e-3 ,blade.sections) + -1e-3;
    offset.y = zeros([1,blade.sections]);
else
    offset.z = blade.zoffset.*ones([1,blade.sections]);
    offset.y = zeros([1,blade.sections]);
end
end
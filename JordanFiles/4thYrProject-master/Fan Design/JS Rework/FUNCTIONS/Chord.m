function [blade] = Chord(blade, radius, sections, rc, rh, rm)
%% Return chord distribution and number of blades
% 

% Midline chord
blade.m.chord        = (rc - rh)./blade.AR;

% Midline pitch
blade.m.pitch        = blade.m.chord .* blade.m.pitchchord;

% Blade number
blade.Nraw         = 2 * pi * rm ./ blade.m.pitch;

% Spanwise pitch
blade.sec.pitch        = 2*pi .* sections ./ blade.Nraw;
blade.span.pitch        = 2*pi .* radius ./ blade.Nraw;

% Spanwise chord
blade.sec.chord        = blade.sec.pitch ./ blade.sec.pitchchord;
blade.span.chord        = blade.span.pitch./ blade.span.pitchchord;

blade.pitchchord = blade.m.pitchchord;

blade.N = ceil(blade.Nraw);
end
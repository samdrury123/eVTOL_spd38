function [ c_r , c_s ] = ChordLengths( DF_max, angles, r_c, r_h, sc_r, sc_s, AR_rotor, AR_stator )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% 
r_m             = (r_c + r_h)/2;
r_temp          = linspace(r_h,r_c,5)';

% Chord Lengths
r.chord         = (r_c-r_h)/AR_rotor;
s.chord         = (r_c-r_h)/AR_stator;

sc_r            = sc_r * ones(5,1);
sc_s            = sc_s * ones(5,1);

for i = 1:5
% Angles
a1              = angles(i,2);
a2              = angles(i,4);
a3              = angles(i,3);
a4              = angles(i,1); 

% Rotor Pitch to Chord Ratio
DF_r(i)           = 1 - cos(a1)/cos(a2) + (cos(a1) * abs(tan(a2) - tan(a1)) * sc_r(i))/2;

% Stator Pitch to Chord Ratio
DF_s(i)           = 1 - cos(a3)/cos(a4) + (cos(a3) * abs(tan(a4) - tan(a3)) * sc_s(i))/2;

% Set DF to maximum and recalculate the pitch-to-chord ratio
if DF_r(i) >= DF_max
    DF                  = DF_max;
    [sc_r(i)]           = (DF + cos(a1)/cos(a2) - 1)/(0.5 * cos(a1) * abs(tan(a1)-tan((a2))));
end

if DF_s(i) >= DF_max
    DF                  = DF_max;
    [sc_s(i)]           = (DF + cos(a3)/cos(a4) - 1)/(0.5 * cos(a3) * abs(tan(a3)-tan((a4))));
end

end

% Pitch (Blade Spacing)
sc_r_m          = sc_r(3);
sc_s_m          = sc_s(3);

s_r             = sc_r_m * r.chord; % m
s_c             = sc_s_m * s.chord; % m
Circumference   = 2 * pi * r_m; % m Outer circumference

% Number of blades
rblades         = ceil(Circumference / s_r);
sblades         = ceil(Circumference / s_c);

theta_r         = 2 * pi / rblades;
theta_s         = 2 * pi / sblades;

pitch_r         = r_temp * theta_r;
pitch_s         = r_temp * theta_s;

% Chord length at different radii
c_r         = pitch_r ./ sc_r;
c_s         = pitch_s ./ sc_s;

end


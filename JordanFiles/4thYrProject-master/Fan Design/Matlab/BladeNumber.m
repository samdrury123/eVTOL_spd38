function [rblades,sblades,sc_r,sc_s,c_r,c_s] = BladeNumber(AR_rotor,AR_stator,DF_r,DF_s,r_h,r_c,angles,N)
    % Calculates the number of Blades using Lieblein
    
    %% Blade Number
    % Find radial positions
    r_m             = (r_c+r_h)/2;
    r_temp          = linspace(r_h,r_c,N)';
    % Chord Lengths
    r.chord         = (r_c-r_h)/AR_rotor;
    s.chord         = (r_c-r_h)/AR_stator;
    % Initialise pitch to chord ratios
    sc_r            = ones(5,1);
    sc_s            = ones(5,1);
    for i = 1:N
        
        % Angles
        a1              = angles(i,2);
        a2              = angles(i,4);
        a3              = angles(i,3);
        a4              = angles(i,1);
        
        
        % Rotor Pitch to Chord Ratio
        [sc_r(i)]       = (DF_r(i)+cos(a1)/cos(a2)-1) / (0.5 * cos(a1) * ( tan(abs(a1))-tan(abs(a2)) ) );
        
        
        % Stator Pitch to Chord Ratio
        [sc_s(i)]       = (DF_s(i)+cos(a3)/cos(a4)-1) / (0.5*cos(a3)*abs(tan(a3)-tan((a4))));
        
        
    end
    % Pitch (Blade Spacing)
    N_mid           = ceil(N/2);            % Index of mid-blade radius
    sc_r_m          = sc_r(N_mid);          % Rotor pitch to chord ratio
    sc_s_m          = sc_s(N_mid);          % Stator pitch to chord ratio
    s_r             = sc_r_m*r.chord;       % Rotor blade pitch(m)
    s_s             = sc_s_m*s.chord;       % Stator blade pitch(m)
    Circumf         = 2*pi*r_m;             % Circumference at mid-radius (m)
    % Number of blades
    rblades         = ceil(Circumf/s_r);    % Round up the no. of blades
    sblades         = ceil(Circumf/s_s);    % Round up the no. of blades
    % Angle between blade centres
    theta_r         = 2*pi/rblades;
    theta_s         = 2*pi/sblades;
    % Blade pitches
    pitch_r         = r_temp*theta_r;
    pitch_s         = r_temp*theta_s;
    % Chord length at different radii
    c_r             = pitch_r./sc_r;
    c_s             = pitch_s./sc_s;
end
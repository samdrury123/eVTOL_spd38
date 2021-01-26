function [angle, V, R, S, delta, carter] = Carter(angle, V, rc, rh, rm)
% Return angles and velocities based on carters rule
global crdfflag

n = length(angle.b1);

if crdfflag == 1
    DFtip = 0.3;
    DFhub = 0.3;
else
    DFtip = 0.25;
    DFhub = 0.4;
end

DFmid = DFhub - (DFhub-DFtip)*(rm-rh)/(rc-rh);

if n == 1
    DF = DFmid; % Diffusion Factor
%     DF = 0.45;
elseif n == 5
    DF = linspace(DFhub,DFtip,n);
%     DF = 0.45;
else
    DF = linspace(DFhub,DFtip,n);
%     DF = 0.45;
end

if crdfflag == 1
    % Carter's rule
    R.m = 0.23 - angle.b2 ./ 500;
    S.m = 0.23 + (angle.a2 - angle.b3) ./ 500;

    % Pitch to chord
    R.pitchchord = (DF - (1 - (V.rel2 ./ V.abs1))) ./ (V.theta2 ./ (2.*V.rel1));
    S.pitchchord = (DF - (1 - (V.rel3 ./ V.abs2))) ./ (V.theta2 ./ (2.*V.abs2));

    % Metal angles
    angle.chi1 = -angle.b1+angle.i1;
    angle.chi2 = (-angle.b2 - angle.chi1 .* R.m .* R.pitchchord .^ (0.5)) ./ (1 - R.m.*R.pitchchord.^(0.5));
    
    angle.chi3 = angle.a2+angle.i3;
    angle.chi4 = (angle.a3 - abs(angle.chi3) .* S.m .* S.pitchchord .^ (0.5)) ./ (1 - S.m.*S.pitchchord.^(0.5));

    % Deviation
    delta.S = angle.chi4 - angle.a3;
    delta.R = angle.chi2 + angle.b2;

    % Limit deviation
    %delta.S(delta.S > 0) = -10; delta.S(delta.S < -10) = -10;
    %delta.R(delta.R > 0) = -10; delta.R(delta.R < -10) = -10;
    % delta.R = linspace(1,0.75,length(delta.R)) .* delta.R;
    % if max(abs(delta.R)) > 15
    %     delta.R = -5*ones(size(delta.R));
    % end
    % Recalculate exit angle
    angle.chi2 = -angle.b2 + delta.R;
    angle.chi4 = abs(angle.a3) + delta.S;

    carter.S = S.m;
    carter.R = R.m;

    R.delta = delta.R;
    S.delta = delta.S;
else
    % Carter's rule
    R.m = 0.23 - angle.b2 ./ 500;
    S.m = 0.23;

    % Pitch to chord
    R.pitchchord = (DF - (1 - (V.rel2 ./ V.rel1))) ./ (V.theta2 ./ (2.*V.rel1));
    S.pitchchord = (DF - (1 - (V.x ./ V.abs2))) ./ (V.theta2 ./ (2.*V.abs2));

    % Metal angles
    angle.chi1 = abs(angle.b1)+angle.i1;
    angle.chi2 = (-angle.b2 - angle.chi1 .* R.m .* R.pitchchord .^ (0.5)) ./ (1 - R.m.*R.pitchchord.^(0.5));
    angle.chi3 = angle.a2+angle.i3;
    angle.chi4 = (angle.a3 - abs(angle.chi3) .* S.m .* S.pitchchord .^ (0.5)) ./ (1 - S.m.*S.pitchchord.^(0.5));

    % Deviation
    delta.S = angle.chi4 - angle.a3;
    delta.R = angle.chi2 + angle.b2;

    % Limit deviation
    %delta.S(delta.S > 0) = -10; delta.S(delta.S < -10) = -10;
    %delta.R(delta.R > 0) = -10; delta.R(delta.R < -10) = -10;
    % delta.R = linspace(1,0.75,length(delta.R)) .* delta.R;
    % if max(abs(delta.R)) > 15
    %     delta.R = -5*ones(size(delta.R));
    % end
    % Recalculate exit angle
    angle.chi2 = -angle.b2 + delta.R;
    angle.chi4 = abs(angle.a3) + delta.S;

    carter.S = S.m;
    carter.R = R.m;

    R.delta = delta.R;
    S.delta = delta.S;
end
end

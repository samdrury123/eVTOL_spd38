function [R, S] = Assemble(R,S,angle,sections,radius, rc, rh, rm)
%% Assemble variables required for blade generation
% [R, S] = Assemble(R,S,angle,sections,radius);

R.sec.chi1 = angle.sec.chi1;
R.sec.chi2 = angle.sec.chi2;

S.sec.chi1 = angle.sec.chi3;
S.sec.chi2 = angle.sec.chi4;

R = Chord(R, radius, sections, rc, rh, rm);
S = Chord(S, radius, sections, rc, rh, rm);

R.sec.radius = sections;
S.sec.radius = sections;

R.z = 0.0125;
S.z = 0.037318;

R.sweep = 10;
S.sweep = 0;

R.lean = 40;
S.lean = 0;
end
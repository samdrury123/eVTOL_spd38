function P = calc_powereta(P)
%% P = calc_powereta(P)
% Calculate power efficiency = M_f_{measured}/M_f_{ideal}
% P: struct with P.FOM, P.sig

P.powereta = P.FOM ./ sqrt(2*P.sig);
end
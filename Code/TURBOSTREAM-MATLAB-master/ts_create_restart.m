function [g] = ts_create_restart(g, nstep, ilos)
% TS_CREATE_RESTART  Apply finer restart options to a TURBOSTREAM mesh
%
%   g = TS_CREATE_RESTART(g,nstep,poisson_nstep,initial_value)
%
%   g - cell array of TURBOSTREAM solution
%   nstep - number of main loop iterations to use
%   ilos - number to select turbulence model: 0, 1 or 2 for euler, mixing length or S-A

% Set poisson to restart
g{1}.av.poisson_nstep = 0;
g{1}.av.poisson_restart = 1;

% Set turbulence model
g{1}.av.restart = 1;
g{1}.av.ilos = ilos;

% Set mixing length for first part of inlet block
for bid = ts_get_bids(g)
    for pid = ts_get_pids(bid,g)
        kind = g{bid+1}.patch{pid+1}.attribute.kind;
        if ilos == 2 && (kind == 0 || kind == 4)
            g{bid+1}.bv.nimixl = 8;
        else
            g{bid+1}.bv.nimixl = 0;
        end
    end
end

% Lower smoothing for more accurate solution
g{1}.av.cfl = 0.4;
g{1}.av.dampin = 10;
g{1}.av.sfin = 0.5;
g{1}.av.facsecin = 0.005;
g{1}.av.rfmix = 0.2;

% Change step numbers and average over final part of solution
g{1}.av.nstep = nstep;
% g{1}.av.nstep_save_start = round(4*nstep/5);
g{1}.av.nstep_save_start = nstep - 10000;

% Write out entropy generation rates
g{1}.av.write_egen = 1;


end
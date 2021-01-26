function s = exp_dpi_initialise(dpi,P)
% EXP_DPI_INITIALISE  Open a connection to the dpi and set zero
%
%   dpi - input dpi data structure
%   P - desired control pressure in pascals
%   s - ouput serial connection
%
%   Example dpi structure:
%   dpi.port = 'COM1'; dpi.rate = 9600; dpi.term = 'LF';

% Open a serial connection to the dpi
s = open_serial(dpi);

% Put in remote mode
fprintf(s,'R1\n');

% Auto rate
fprintf(s,'J1\n');

% Use pascals
fprintf(s,'S3,U1\n');

% Zero unit
fprintf(s,'O1\n');

% Control at desired pressure in pascals
fprintf(s,'C1\n');
fprintf(s,['P=' num2str(round(P)) '\n']);

end
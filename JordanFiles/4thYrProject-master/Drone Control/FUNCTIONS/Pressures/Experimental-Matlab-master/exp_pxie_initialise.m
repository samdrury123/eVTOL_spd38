function s = exp_pxie_initialise(pxie)
% Open specified channels on national instruments PXIe aquisition system

% Default all channels
if isfield(pxie,'chan') == 0
    pxie.chan = 1:pxie.nchan;
end

% Default maximum rate
if isfield(pxie,'rate') == 0
    pxie.rate = 200000;
end

% Default to 1 second logging
if isfield(pxie,'time') == 0
    pxie.time = 1;
end

% Create session
s = daq.createSession(pxie.sesh);

% Add voltage channels
if strcmp(pxie.type,'V') == 1
    for n = 1:length(pxie.chan)
        ch = s.addAnalogInputChannel(pxie.dev{pxie.chan(n)},pxie.ai(pxie.chan(n)),'Voltage');
        ch.Coupling = pxie.coupling{n};
        ch.TerminalConfig='SingleEnded';
    end
end

% Add RTD channels
if strcmp(pxie.type,'RTD') == 1
    for n = 1:length(pxie.chan)
        ch = s.addAnalogInputChannel(pxie.dev{pxie.chan(n)},pxie.ai(pxie.chan(n)),'RTD');
        ch.RTDType = pxie.class{n}; ch.RTDConfiguration = pxie.wire{n}; 
        ch.R0 = pxie.R0(n); ch.Units = pxie.unit{n};
    end
end

% Set rate & duration
s.Rate = pxie.rate;
s.DurationInSeconds = pxie.time;


end
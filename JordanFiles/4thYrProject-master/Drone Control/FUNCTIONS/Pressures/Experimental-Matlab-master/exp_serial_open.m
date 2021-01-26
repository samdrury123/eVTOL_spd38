function s = exp_serial_open(options)
% Open a connection to a serial port

s = serial(options.port);
s.BaudRate = options.rate;
s.Terminator = options.term;
fopen(s);

end
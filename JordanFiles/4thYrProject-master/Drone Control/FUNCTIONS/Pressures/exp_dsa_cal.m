function dsa = exp_dsa_cal(s,dsa)
% Get a calibration map and constants from a DSA

% Get temperature gain and offsets
flushinput(s); flushoutput(s);
dsa.T_a = zeros(dsa.nchan,1); dsa.T_c = zeros(dsa.nchan,1); 

fprintf(s, 'LIST G\n');
wait_bytes(s)
for i = 1:dsa.nchan
    data = regexp(fscanf(s),'\s+','split');
    dsa.T_a(i) = str2double(data{3});
end

fprintf(s, 'LIST O\n');
wait_bytes(s)
for i = 1:dsa.nchan
    data = regexp(fscanf(s),'\s+','split');
    dsa.T_c(i) = str2double(data{3});
end

% Get master pressure calibration matrix - [T Counts P];
flushinput(s); flushoutput(s);
nT = 4; nP = 9; 
dsa.M = zeros(dsa.nchan,nT,nP,3);
for i = 1:dsa.nchan
    fprintf(s, ['LIST M 10 40 ' num2str(i) '\n']);
    exp_serial_wait(s)
    b_old = 0; b = 1;
    while b_old ~= b
        b_old = b;
        b = s.BytesAvailable;
        pause(0.5);
    end
    for n = 1:nT
        for m = 1:nP
            data = regexp(fscanf(s),'\s+','split');
            dsa.M(i,n,m,1) = str2double(data{2});
            dsa.M(i,n,m,2) = str2double(data{5});
            dsa.M(i,n,m,3) = str2double(data{4});
        end
    end
end

end
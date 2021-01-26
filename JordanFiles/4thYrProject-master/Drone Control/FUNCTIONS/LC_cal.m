importLCcal;

bar = 309.22;
w0 = 50.13;
w1 = 50.01;
w2 = 49.99;
w3 = 50.04;

weights = [bar bar+w0 bar+w0+w1 bar+w0+w1+w2 bar+w0+w1+w2+w3];
readings = [mean(value(24:40)) mean(value(47:56)) mean(value(60:70)) mean(value(74:85)) mean(value(91:102))];

figure; subplot(2,1,1); hold on; plot(value);
subplot(2,1,2); hold on; plot(weights, readings); plot(LC_cal(readings), readings);
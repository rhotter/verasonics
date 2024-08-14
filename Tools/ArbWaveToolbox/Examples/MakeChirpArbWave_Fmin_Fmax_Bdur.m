% Example for creating a chirp waveform

clc

%% Generate the chirp and plot
fstart = 3; % MHz
fstop = 7; % MHz
duration = 1; % us

Fs = 250;
t= 0:1/Fs:duration;

% with hanning window for better performance
w = hann(length(t));
y = chirp(t,fstart,duration,fstop);
y = w'.* y;
ypad = [zeros(1,50), y, zeros(1,100)];

figure,plot(y),hold on, plot(ypad), grid on, hold off

%% save in files
% matfile
ArbWave = ypad';
matfilename = ['Chirp_', num2str(duration),'us_AM_',num2str(fstart),'to',num2str(fstop),'MHz.mat'];
save (matfilename, 'ArbWave')

% % text file
% fileID1 = fopen(['Chirp_', num2str(duration),'us_AM_',num2str(fstart),'to',num2str(fstop),'MHz.txt'],'w');
% fprintf(fileID1,'%0.9f\r\n',y');
% fclose(fileID1);

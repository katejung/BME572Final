clc;
clear all;
close all;

load MAp.mat
load MHp.mat

% training = [1 3 5 7];
training = [3];
test = [2 4 6 8];

f = 15000;
threshold = 0;
start = 10;
m = MAp;

MApt = MAp(training(1),f*start:end);
for i = 1:size(training,2)-1
    MApt = horzcat(MApt,MAp(training(i+1),f*start:end));
end
% Butterworth bandpass filter
d = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',300,'HalfPowerFrequency2',3000, ...
    'SampleRate',15000);
MAplp = filter(d,MApt);
% [b, a] = butter(2,[300/7500 3000/7500],'bandpass');
% MAplp = filter(b,a,MApt);

x1 = 32; % half the window size. Window size = 2x+1

figure;
plot(linspace(1./f,size(MAplp,2)./f,size(MAplp,2)),MAplp);
ylim([-200 200]);
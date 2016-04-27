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
start = 11.5;
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

x1 = 20; % half the window size. Window size = 2x+1
pks = zeros(1,size(MAplp,2));
locs = zeros(1,size(MAplp,2));
% (max(MAplp(i-x1:i+x1)) - min(MAplp(i-x1:i+x1))) >= 20
%%
for i = 1+x1:size(MAplp,2)-x1
    if (max(MAplp(i-x1:i+x1)) == MAplp(i)) && (max(MAplp(i-x1:i+x1)) - min(MAplp(i-x1:i+x1)) > 10)
        pks(i) = MAplp(i);
        locs(i) = i;
    end
end
pks = pks(pks~=0);
locs = locs(locs~=0);
    
thrclass = pks>0;
pkst = pks.*(thrclass);
pkst = pkst(pkst~=0);
locst = locs.*(thrclass);
locst = locst(locst~=0);
%%
figure;
plot(linspace(1./f,size(MAplp,2)./f,size(MAplp,2)),MAplp);
ylim([-200 200]);
hold on;
plot(locst./f,pkst,'.');
hold off;

%% Histogram idea to separate oscillations from others (possible use)
% histogramsize = 1;
% amp = linspace(min(MAplp),max(MAplp),round((max(MAplp)-min(MAplp))/histogramsize));
% count = zeros(1,size(amp,2));
% for i = 1:size(amp,2)-1
%     count(i) = size(MAplp((MAplp>amp(i)) == (MAplp<amp(i+1))),2);
% end
% figure;
% plot(amp,count)

%%
x2 = x1; % Size of window of data
if x2 <= x1
    data = zeros(size(locst,2),2*x2+1);
    for i = 1:size(locst,2)
        data(i,:) = MAplp((locst(i) - x2):(locst(i) + x2));
    end
else
    error('PCA data window bigger. Change the hard code.');
end

[COEFF, SCORE, LATENT] = pca(data);
figure;
plot(LATENT(1:10));

pcadata = data*COEFF;
figure;
scatter(pcadata(:,1),pcadata(:,2));

figure;
scatter3(pcadata(:,1),pcadata(:,2),pcadata(:,3));

%%
figure;
hold on;
% for i = 1:size(data,1)
for i = 1:size(data,1)
plot(data(i,:));
end
hold off;
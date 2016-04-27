% Project_main_JS.m
clear all; clc;
% Load data
load MAp.mat;
load MHp.mat;

%% Plotting settings
linespec = {'linestyle' '-';'linestyle' ':';'linestyle' '-.';'linestyle' '--'}; % Line specs for gray-scale
markerspec = {'linestyle' 'o';'linestyle' '*';'linestyle' '+';'linestyle' '.'}; % Marker specs for gray-scale
color = {'color' 'blue'; 'color' 'red'; 'color' 'green'; ...
    'color' 'cyan'; 'color' 'black'; 'color' 'yellow'};
% fontsize{1,:} for title, fontsize{2,:} for axis, fontsize{3,:} for legend
fontsize = {'Fontsize' 13 'Fontweight' 'bold'; ...
    'Fontsize', 10, 'Fontweight', 'bold'; ...
    'Fontsize', 11, 'Fontweight', 'bold'};
% linespec2{1,:} for line, linespec2{2,:} for marker
linespec2 = {'linewidth' 2; 'MarkerSize' 6};
% Count the number of figures created
fignum = 1;

%% Preprosess the data
% preprocess parameters
samplerate = 15000;
t = 0:(1/samplerate):60; t = t(1:end-1);
baselinec = 0; % Remove first few seconds
intensec = 11.5; % Cutoff for the intense section
poststimulusc = 0; % Cutoff for the post-stimulus section

% Preprocess the data
[MAp_baseline,MAp_intense,MAp_poststimulus] = compartmentize_JS(MAp,baselinec,intensec,poststimulusc);
[MAp_baseline] = filtering_JS(MAp_baseline,0,1); [MAp_intense] = filtering_JS(MAp_intense,1,1); [MAp_poststimulus] = filtering_JS(MAp_poststimulus,0,1);
MAp_filtered = [MAp_baseline MAp_intense MAp_poststimulus];
[MHp_baseline,MHp_intense,MHp_poststimulus] = compartmentize_JS(MHp,baselinec,intensec,poststimulusc);
[MHp_baseline] = filtering_JS(MHp_baseline,0,0); [MHp_intense] = filtering_JS(MHp_intense,1,0); [MHp_poststimulus] = filtering_JS(MHp_poststimulus,0,0);
MHp_filtered = [MHp_baseline MHp_intense MHp_poststimulus];

%% Visualize
fig = figure(fignum); clf(fig)
plot(t(1+samplerate*baselinec:end),MAp_filtered,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MAp filtered data',fontsize{1,:})
fignum = fignum + 1;

fig = figure(fignum); clf(fig)
subplot(3,1,1)
plot(t(1+samplerate*baselinec:samplerate*10),MAp_baseline,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MAp filtered data: Baseline',fontsize{1,:})
subplot(3,1,2)
plot(t(samplerate*10+1:samplerate*intensec),MAp_intense,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MAp filtered data: Intense',fontsize{1,:})
subplot(3,1,3)
plot(t(samplerate*intensec+1:end-poststimulusc),MAp_poststimulus,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MAp filtered data: Post-stimulus',fontsize{1,:})
fignum = fignum + 1;

fig = figure(fignum); clf(fig)
plot(t(1+samplerate*baselinec:end),MHp_filtered,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MHp filtered data',fontsize{1,:})
fignum = fignum + 1;

fig = figure(fignum); clf(fig)
subplot(3,1,1)
plot(t(1+samplerate*baselinec:samplerate*10),MHp_baseline,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MHp filtered data: Baseline',fontsize{1,:})
subplot(3,1,2)
plot(t(samplerate*10+1:samplerate*intensec),MHp_intense,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MHp filtered data: Intense',fontsize{1,:})
subplot(3,1,3)
plot(t(samplerate*intensec+1:end-poststimulusc),MHp_poststimulus,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MHp filtered data: Post-stimulus',fontsize{1,:})
fignum = fignum + 1;

%% spike detection
[pks,locs,peakinterval] = spike_detection_JS(MAp_baseline);
fig = figure(fignum); clf(fig)
plot(t(1+samplerate*baselinec:samplerate*10),MAp_baseline(1,:),linespec{1,:},linespec2{1,:});
hold on
plot(t(locs{1,1}+samplerate*baselinec),MAp_baseline(1,locs{1,1}),'rs','MarkerFaceColor','r');
axis([2.03 2.05 -25 50])
fignum = fignum  +1;

fig = figure(fignum); clf(fig)
hist(peakinterval{1,1},20)
fignum = fignum  +1;

fig = figure(fignum); clf(fig)
plot(t,MAp(1,:),linespec{1,:},linespec2{1,:});
axis([2.03 2.05 -65 15])
fignum = fignum  +1;

[pks,locs,peakinterval] = spike_detection_JS(MAp_poststimulus);
fig = figure(fignum); clf(fig)
plot(t(samplerate*intensec+1:end-poststimulusc),MAp_poststimulus(1,:),linespec{1,:},linespec2{1,:});
hold on
plot(t(locs{1,1}+samplerate*intensec+1),MAp_poststimulus(1,locs{1,1}),'rs','MarkerFaceColor','r');
% axis([14.98 15 -25 50])
fignum = fignum  +1;

fig = figure(fignum); clf(fig)
plot(t,MAp(1,:),linespec{1,:},linespec2{1,:});
% axis([14.98 15 -370 -340])
fignum = fignum  +1;

[pks,locs,peakinterval] = spike_detection_JS(MHp_poststimulus);
fig = figure(fignum); clf(fig)
plot(t(samplerate*intensec+1:end-poststimulusc),MHp_poststimulus(1,:),linespec{1,:},linespec2{1,:});
hold on
plot(t(locs{1,1}+samplerate*intensec+1),MHp_poststimulus(1,locs{1,1}),'rs','MarkerFaceColor','r');
% axis([14.98 15 -25 50])
fignum = fignum  +1;

fig = figure(fignum); clf(fig)
plot(t,MHp(1,:),linespec{1,:},linespec2{1,:});
% axis([14.98 15 -370 -340])
fignum = fignum  +1;


%% Find peak
MAp_train_baseline = MAp_baseline(1:4,:);
MHp_train_baseline = MHp_baseline(1:4,:);
MAp_train_poststimulus = MAp_poststimulus(1:4,samplerate*30+1:end-poststimulusc);
MHp_train_poststimulus = MHp_poststimulus(1:4,samplerate*30+1:end-poststimulusc);
MAp_train_intense = MAp_intense(1:4,:);
MHp_train_intense = MHp_intense(1:4,:);

dataInterest = MAp_train_baseline;
[peakinterval,validPeakCell, validTroughCell, validPeakLocCell, validTroughLocCell] = extractFeatures_height(dataInterest);
dataMedian = 0;
for m = 1:4
    dataMedian = dataMedian + median(abs(dataInterest(m,:))/0.6745);
end
dataMedian = dataMedian/4;

peaks = [];
troughs = [];
peakT = [];
troughT = [];
for i = 1:4
    peaks = [peaks validPeakCell{i,1}];
    troughs = [troughs validTroughCell{i,1}];
    peakT = [peakT validPeakLocCell{i,1}];
    troughT = [troughT validTroughLocCell{i,1}];
end

num = size(peaks,2);
for j = 1:num
    heights(j) = peaks(j)-troughs(j);
    widths(j) = troughT(j)-peakT(j);
end

% overlapping signals are not
figure
for i = 1:4
    windowLength = 100;
    peakTime_trial = validPeakLocCell{i,1};
    peak_trial = validPeakCell{i,1};
    trough_trial = validTroughCell{i,1};
    data_trial = dataInterest(i,:);
    for j  = 1:size(peakTime_trial,2)
        if peakTime_trial(j)+44 < size(peakTime_trial,2)
        if peak_trial(j) > 7*dataMedian 
            plot(data_trial(peakTime_trial(j)-19:peakTime_trial(j)+44),'b');
        else
            plot(data_trial(peakTime_trial(j)-19:peakTime_trial(j)+44),'r');
        end
        end
    end
    hold on  
end
xlabel('time spike at 20 sec')
ylabel('voltage')
title('Spike Alignment')

figure
scatter(troughs,peaks)
xlabel('Spike Minimum (\muV)')
ylabel('Spike Maximum (\muV)')
title('Cluster Maximum vs. Minimum using threshold voltage height')
num = size(peaks,1);
figure
for k = 1:num
    if peaks(k) > 7*dataMedian
        plot(widths(k), heights(k), 'ro')
    else
        plot(widths(k), heights(k), 'bo')
    end
    hold on
end

xlabel('spike width')
ylabel('spike height')
title('Cluster height vs. width using threshold voltage height')

% sample data
figure
plot(dataInterest(1,:));
hold on
plot(validPeakLocCell{1,1}, validPeakCell{1,1},'bo')
hold on
plot(validTroughLocCell{1,1}, validTroughCell{1,1},'mo')
xlabel('time')
ylabel('voltage')
title('sample spike detection')






%% PCA
K =2;% number of principal Comonents

X(:,1) = troughs;
X(:,2) = peaks;
[n,d]=size(X);
C=cov(X);
[V U]=eig(C); %V = eigenvectors U = eigenvalues
L=diag(U); L', U % principal component directions
[sortedU, index] = sort(L, 'descend');
Xproj_PCA = zeros(d,K);
for j = 1:K
    Xproj_PCA(:,j) = V(:,index(j));
end
Y_trained_PCA = X*Xproj_PCA;
Xproj_PCA; % the first column is the 1st PC and the second column is the 2nd PC

% calulate 1st PC line
Xh1_PCA=-35:0.5:5;
Yh1_PCA = Xh1_PCA.*Xproj_PCA(2,1)./Xproj_PCA(1,1);
% calculated 2nd PC line
Xh2_PCA =-35:0.5:5;
Yh2_PCA = Xh2_PCA.*Xproj_PCA(2,2)./Xproj_PCA(1,2);
figure
plot(X(:,1), X(:,2), 'b^', Xh1_PCA,Yh1_PCA,'g-',Xh2_PCA, Yh2_PCA, 'r-')
legend('data', '1st PCA axis', '2nd PCA axis')
figure

scatter(Y_trained_PCA(:,1), Y_trained_PCA(:,2), 'bo')

title('Encoded Training Data')
xlabel(['1st PC [' num2str(Xproj_PCA(1,1)) ',' num2str(Xproj_PCA(2,1)) ']'])
ylabel(['2nd PC [' num2str(Xproj_PCA(1,2)) ',' num2str(Xproj_PCA(2,2)) ']'])


%%

clusterA_peak = [];
clusterA_trough = [];
for i = 1:size(Y_trained_PCA,1)
    I = find(Y_trained_PCA(:,1)>20);
end

for j = 1:4
    clearvars matrix I_Apeak I_Bpeak
    matrix1 = validTroughCell{j,1};
    matrix2 = validPeakCell{j,1};
	matrix(:,1) = matrix1;
    matrix(:,2) = matrix2;
    peakTime = validPeakLocCell{j,1};
    pca_trained = matrix*Xproj_PCA;
    figure
    scatter(pca_trained(:,1), pca_trained(:,2), 'bo')
    I_Apeak = find(pca_trained(:,1)>20);
    I_Bpeak = find(pca_trained(:,1)<=20);
    Apeak_index{j,1} = I_Apeak;
    Bpeak_index{j,1} = I_Bpeak;
end

%% spike alighment according to PCA
for i = 1:4
    windowLength = 100;
    peakTime_trial = validPeakLocCell{i,1};
    peak_trial = validPeakCell{i,1};
    ApeakI = Apeak_index{i,1};
    BpeakI = Bpeak_index{i,1};
    data_trial = dataInterest(i,:);
    for j = 1:size(ApeakI,1)
        peakStart = ApeakI(j);
         plot(data_trial(peakTime_trial(peakStart)-19:peakTime_trial(peakStart)+44),'b');
    end
    for k = 1:size(BpeakI,1)
        peakStart = BpeakI(k);
         plot(data_trial(peakTime_trial(peakStart)-19:peakTime_trial(peakStart)+44),'r');
    end
    
    hold on  
end
xlabel('time spike at 20 sec')
ylabel('voltage')
title('Spike Alignment with PCA classification')









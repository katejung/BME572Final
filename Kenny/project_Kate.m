%%
% Project_main_JS.m
clear all; clc;
% Load data
load MAp.mat;
load MHp.mat;

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

% Visualize
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

%% Feature Analysis
Fs  = 15000;
t = t(1+samplerate*baselinec:samplerate*10);
MAp_train_baseline = MAp_baseline(1:4,:);
MHp_train_baseline = MHp_baseline(1:4,:);
MAp_train_poststimulus = MAp_poststimulus(1:4,samplerate*30+1:end-poststimulusc);
MHp_train_poststimulus = MHp_poststimulus(1:4,samplerate*30+1:end-poststimulusc);
MAp_train_intense = MAp_intense(1:4,:);
MHp_train_intense = MHp_intense(1:4,:);
figure 
subplot(4,1,1)

plot(MAp_poststimulus(1, :))
title('full postStimulus')
subplot(4,1,2)
plot(MAp_train_poststimulus(1,:))
title('post stimulus starting from 30 sec')
subplot(4,1,3)
plot(MHp_poststimulus(1,:))
subplot(4,1,4)
plot(MHp_train_poststimulus(1,:))
%% height
clc
clearvars -except MAp_train_baseline MHp_train_baseline MAp_train_poststimulus MHp_train_poststimulus MAp_train_intense MHp_train_intense
dataInterest = MHp_train_baseline;
[validPeakCell, validTroughCell, validPeakLocCell, validTroughLocCell] = extractFeatures_height(dataInterest);
%%
% 
% [peaks, troughs,  peakT, troughT] = extractFeatures_threshold(dataInterest);
% figure
% peaks = [];
% troughs = [];
% peakT = [];
% troughT = [];
% for i = 1:4
%     peaks = [peaks validPeakCell{i,1}];
%     troughs = [troughs validTroughCell{i,1}];
%     peakT = [peakT validPeakLocCell{i,1}];
%     troughT = [troughT validTroughLocCell{i,1}];
% end
% 
% % for k = 1:4
% %     peak_trial = validPeakCell{i,1};
% %     for l = 1:size(peak_trial, 2)
% %         if peak_trial(l) 
% 
% scatter(troughs,peaks)
% xlabel('Spike Minimum (\muV)')
% ylabel('Spike Maximum (\muV)')
% title('Cluster Maximum vs. Minimum using threshold voltage')

%% height
[peaks, troughs, peakT, troughT] = extractFeatures_height(dataInterest);
num = size(peaks,2);
for j = 1:num
    heights(j) = peaks{j}-troughs{j};
    widths(j) = troughT{j}-peakT{j};
end

dataMedian = 0;
for m = 1:4
    dataMedian = dataMedian + median(abs(dataInterest(m,:))/0.6745);
end
dataMedian = dataMedian/4;



%% spike line up
% overlapping signals are not
figure(1);
% for i = 1:4
    windowLength = 100;
%     peakTime_trial = validPeakLocCell(1,i);
%     peak_trial = validPeakCell(1,i);
%     trough_trial = validTroughCell(1,i);
%     data_trial = dataInterest(i,:);
for i  = 1:size(peakT,2)
    if peakT(i)+44 < size(peakT(i),2)
        if peaks(i) > 20
            set(gca,'Ydata',dataInterest((peakT(i)-19):(peakT(i)+44)));
        else
            set(gca,'Ydata',dataInterest((peakT(i)-19):(peakT(i)+44)));
        end
    end
%     end
    hold on  
end
xlabel('time spike at 20 sec')
ylabel('voltage')
title('Spike Alignment')
hold off
%%
figure
scatter(troughs,peaks)
xlabel('Spike Minimum (\muV)')
ylabel('Spike Maximum (\muV)')
title('Cluster Maximum vs. Minimum using threshold voltage height')
num = size(peaks,1);
figure
% for k = 1:num
%     if peaks(k) > 7*dataMedian
%         plot(widths(k), heights(k), 'ro')
%     else
%         plot(widths(k), heights(k), 'bo')
%     end
%     hold on
% end

xlabel('spike width')
ylabel('spike height')
title('Cluster height vs. width using threshold voltage height')

% sample data
figure
plot(dataInterest(1,:));
hold on
plot(validPeakLocCell(1,1), validPeakCell(1,1),'bo')
hold on
plot(validTroughLocCell(1,1), validTroughCell(1,1),'mo')
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
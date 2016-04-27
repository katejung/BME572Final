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
t_baseline = t(1+samplerate*baselinec:samplerate*10);
MAp_train = MAp_baseline(1:4,:);
plot(MAp_train(1,:));

%% min max feature
arrayLength= size(MAp_train, 2);
numset = size(MAp_train, 1);
troughsArrayA = [];
troughsTimeLocA = [];
windowLength = 500;

for i= 1:numset
    peaksArrayA = [];
    peaksTimeLocA = [];
    dataMedian = median(abs(MAp_train(i,:))/0.6745);
    threshold = 4*dataMedian;
    %     for k = 1:round(arrayLength/windowLength)
    %         startIndex = (k-1)*windowLength;
    %         dataWindow = MAp_train(i,startIndex+1:startIndex+windowLength);
    %         sigma = median(abs(dataWindow)/0.6745);
    %         threshold = 2*sigma;
    %         [peaksInWindow,peaksTimeInWindow] = findpeaks(MAp_train(i,:),'MinPeakHeight',threshold);
    %         peaksArray = [peaksArray peaksInWindow];
    %         peaksTimeLoc = [peaksTimeLoc peaksTimeInWindow+startIndex];
    %     end
    [peaksArrayA,peaksTimeLocA] = findpeaks(MAp_train(i,:),'MinPeakHeight',threshold);
    numPeaks = size(peaksArrayA,2);
    troughsArrayA = [];
    troughsTimeLocA = [];
    for j = 1:numPeaks
        
        peakTime = peaksTimeLocA(j);
        [M,I] = min(MAp_train(i, peakTime:min([peakTime+30, size(MAp_train(i,:),2)])));
        troughsArrayA = [troughsArrayA M];
        troughsTimeLocA = [troughsTimeLocA peakTime+I];
    end
    peakinterval{i,1} = diff(peaksTimeLocA);
    peaksArrayB = [];
    peaksTimeLocB= [];
    for m = 1:numPeaks-1
        peakintervalArray = peakinterval{i,1};
        if peakintervalArray(1,m) > 500
            startIndex = peaksTimeLocA(m)+150;
            endIndex = peaksTimeLocA(m+1)-50;
            sigma = median(abs(MAp_train(i,startIndex:endIndex))/0.6745);
            thresholdB = 3.5*sigma;
            [p,t] = findpeaks(MAp_train(i,startIndex:endIndex),'MinPeakHeight',thresholdB);
            peaksArrayB =[peaksArrayB p];
            peaksTimeLocB = [peaksTimeLocB t+startIndex];
        end
        
    end
    troughsArrayB = [];
    troughsTimeLocB = [];
    for l = 1:size(peaksArrayB,2)
        peakTime = peaksTimeLocB(l);
        [M,I] = min(MAp_train(i, peakTime:min([peakTime+30, size(MAp_train(i,:),2)])));
        troughsArrayB = [troughsArrayB M];
        troughsTimeLocB = [troughsTimeLocB peakTime+I];
    end
    peaksA{i,1} = peaksArrayA;
    peakLocA{i,1} = peaksTimeLocA;
    peaksB{i,1} = peaksArrayB;
    peakLocB{i,1} = peaksTimeLocB;
    troughsB{i,1} = troughsArrayB;
    troughLocB{i,1} = troughsTimeLocB;
    troughsA{i,1} = troughsArrayA;
    troughLocA{i,1} = troughsTimeLocA;
end
figure
plot(MAp_train(1,:));
hold on
plot( peakLocA{1,1},peaksA{1,1}, 'ro')
hold on
plot(troughLocA{1,1}, troughsA{1,1},'bo')
hold on
plot(peakLocB{1,1}, peaksB{1,1}, 'mo')
hold on
plot(troughLocB{1,1}, troughsB{1,1},'go')
hold on

%%
figure
peaks = [];
troughs = [];
for i = 1:numset
    peaks = [peaks peaksA{i,1} peaksB{i,1}];
    troughs = [troughs troughsA{i,1} troughsB{i,1}];
end
scatter(troughs,peaks)
xlabel('Spike Minimum (\muV)')
ylabel('Spike Maximum (\muV)')
title('Cluster Maximum vs. Minimum using threshold voltage')

%% min max feature
clearvars peaks troughs peaksArray
arrayLength= size(MAp_train, 2);
numset = size(MAp_train, 1);
troughsArrayA = [];
troughsTimeLocA = [];
windowLength = 500;

for i= 1:numset
    heightThreshold = 7;
    validPeak = [];
    validPeakLoc = [];
    validTrough =[];
    validTroughLoc = [];
    peaksArray = [];
    peaksTimeLoc = [];
    dataMedian = median(abs(MAp_train(i,:))/0.6745);
    threshold = 1.5*dataMedian;
    
    [peaksArray, peaksTimeLoc] = findpeaks(MAp_train(i,:),'MinPeakHeight',threshold);
    numPeaks = size(peaksArray,2);
    
    numPeaks = size(peaksArray,2);
    troughsArray = [];
    troughsTimeLoc = [];
    for j = 1:numPeaks
        peakTime = peaksTimeLoc(j);
        [M,I] = min(MAp_train(i, peakTime:min([peakTime+40, size(MAp_train(i,:),2)])));
        troughsArray = [troughsArray M];
        troughsTimeLoc = [troughsTimeLoc peakTime+I];
    end
    
    
    validIndices = zeros(1, numPeaks);
    for j = 1:numPeaks
        if peaksArray(j) - troughsArray(j) > heightThreshold
            validIndices(j) = 1;
            validPeak = [validPeak peaksArray(j)];
            validPeakLoc = [validPeakLoc peaksTimeLoc(j)];
            validTrough = [validTrough troughsArray(j)];
            validTroughLoc = [validTroughLoc troughsTimeLoc(j)];
        end
    end
    
    % check if A spike gave false positive
    pTemp = [];
    tTemp = [];
    trTemp = [];
    ttrTemp = [];
    
    
    for a = 1:size(validPeak,2)
        keep = true;
        if (validPeak(a) < 4*dataMedian)
            data = MAp_train(i,:);
            time = validPeakLoc(a);
            for b = time:min([time+50 arrayLength])
                if data(b) > 7*dataMedian
                    keep = false;               
                end
            end
        end
        if keep
            pTemp = [pTemp validPeak(a)];
            tTemp = [tTemp validPeakLoc(a)];
            trTemp = [trTemp validTrough(a)];
            ttrTemp = [ttrTemp validTroughLoc(a)];
        end
    end
    validPeak = pTemp;
    validPeakLoc = tTemp;
    validTrough = trTemp;
    validTroughLoc = ttrTemp;
    
    peakinterval{i,1} = diff(peaksTimeLoc);
    peaks{i,1} = peaksArray;
    peaksLoc{i,1} = peaksTimeLoc;
    troughs{i,1} = troughsArray;
    troughsLoc{i,1} = troughsTimeLoc;
    validIndicesCell{i,1} = validIndices;
    validPeakCell{i,1} = validPeak;
    validPeakLocCell{i,1}= validPeakLoc;
    validTroughCell{i,1} = validTrough;
    validTroughLocCell{i,1} = validTroughLoc;
end
figure
plot(MAp_train(1,:));
hold on
plot( validPeakLocCell{1,1},validPeakCell{1,1}, 'ro')
hold on
plot(validTroughLocCell{1,1}, validTroughCell{1,1}, 'mo')

%%
clearvars peaks troughs
figure
peaks = [];
troughs = [];
peakT = [];
troughT = [];
for i = 1:numset
    peaks = [peaks validPeakCell{i,1}];
    troughs = [troughs validTroughCell{i,1}];
    peakT = [peakT validPeakLocCell{i,1}];
    troughT = [troughT validTroughLocCell{i,1}];
end
scatter(troughs,peaks)
xlabel('Spike Minimum (\muV)')
ylabel('Spike Maximum (\muV)')
title('Cluster Maximum vs. Minimum using threshold voltage height')

num = size(peaks,2);
for j = 1:num
   heights(j) = peaks(j)-troughs(j);
   widths(j) = troughT(j)-peakT(j);
end
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
%%
Fs  = 15000;
t_pc = t(samplerate*30+1:end-poststimulusc);

MAp_train = MAp_poststimulus(1:4,samplerate*30+1:end-poststimulusc);
plot(MAp_train(1,:));
%% min max feature
clearvars peaks troughs peaksArray
arrayLength= size(MAp_train, 2);
numset = size(MAp_train, 1);
troughsArrayA = [];
troughsTimeLocA = [];
windowLength = 500;

for i= 1:numset
    heightThreshold = 7;
    validPeak = [];
    validPeakLoc = [];
    validTrough =[];
    validTroughLoc = [];
    peaksArray = [];
    peaksTimeLoc = [];
    dataMedian = median(abs(MAp_train(i,:))/0.6745);
    threshold = 1.5*dataMedian;
    
    [peaksArray, peaksTimeLoc] = findpeaks(MAp_train(i,:),'MinPeakHeight',threshold);
    numPeaks = size(peaksArray,2);
    
    numPeaks = size(peaksArray,2);
    troughsArray = [];
    troughsTimeLoc = [];
    for j = 1:numPeaks
        peakTime = peaksTimeLoc(j);
        [M,I] = min(MAp_train(i, peakTime:min([peakTime+40, size(MAp_train(i,:),2)])));
        troughsArray = [troughsArray M];
        troughsTimeLoc = [troughsTimeLoc peakTime+I];
    end
    
    
    validIndices = zeros(1, numPeaks);
    for j = 1:numPeaks
        if peaksArray(j) - troughsArray(j) > heightThreshold
            validIndices(j) = 1;
            validPeak = [validPeak peaksArray(j)];
            validPeakLoc = [validPeakLoc peaksTimeLoc(j)];
            validTrough = [validTrough troughsArray(j)];
            validTroughLoc = [validTroughLoc troughsTimeLoc(j)];
        end
    end
    
    % check if A spike gave false positive
    pTemp = [];
    tTemp = [];
    trTemp = [];
    ttrTemp = [];
    
    
    for a = 1:size(validPeak,2)
        keep = true;
        if (validPeak(a) < 4*dataMedian)
            data = MAp_train(i,:);
            time = validPeakLoc(a);
            for b = time:min([time+50 arrayLength])
                if data(b) > 7*dataMedian
                    keep = false;               
                end
            end
        end
        if keep
            pTemp = [pTemp validPeak(a)];
            tTemp = [tTemp validPeakLoc(a)];
            trTemp = [trTemp validTrough(a)];
            ttrTemp = [ttrTemp validTroughLoc(a)];
        end
    end
    validPeak = pTemp;
    validPeakLoc = tTemp;
    validTrough = trTemp;
    validTroughLoc = ttrTemp;
    
    peakinterval{i,1} = diff(peaksTimeLoc);
    peaks{i,1} = peaksArray;
    peaksLoc{i,1} = peaksTimeLoc;
    troughs{i,1} = troughsArray;
    troughsLoc{i,1} = troughsTimeLoc;
    validIndicesCell{i,1} = validIndices;
    validPeakCell{i,1} = validPeak;
    validPeakLocCell{i,1}= validPeakLoc;
    validTroughCell{i,1} = validTrough;
    validTroughLocCell{i,1} = validTroughLoc;
end
%%
figure
plot(MAp_train(1,:));
hold on
plot( validPeakLocCell{1,1},validPeakCell{1,1}, 'ro')
hold on
plot(validTroughLocCell{1,1}, validTroughCell{1,1}, 'mo')

%%
clearvars peaks troughs
figure
peaks = [];
troughs = [];
peakT = [];
troughT = [];
for i = 1:numset
    peaks = [peaks validPeakCell{i,1}];
    troughs = [troughs validTroughCell{i,1}];
    peakT = [peakT validPeakLocCell{i,1}];
    troughT = [troughT validTroughLocCell{i,1}];
end
scatter(troughs,peaks)
xlabel('Spike Minimum (\muV)')
ylabel('Spike Maximum (\muV)')
title('Cluster Maximum vs. Minimum using threshold voltage height')

num = size(peaks,2);
for j = 1:num
   heights(j) = peaks(j)-troughs(j);
   widths(j) = troughT(j)-peakT(j);
end
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
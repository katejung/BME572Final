function [peaks, troughs,  peakT, troughT] = extractFeatures_threshold(data)
arrayLength= size(data, 2);
numset = size(data, 1);
troughsArrayA = [];
troughsTimeLocA = [];
windowLength = 500;

for i= 1:numset
    peaksArrayA = [];
    peaksTimeLocA = [];
    dataMedian = median(abs(data(i,:))/0.6745);
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
    [peaksArrayA,peaksTimeLocA] = findpeaks(data(i,:),'MinPeakHeight',threshold);
    numPeaks = size(peaksArrayA,2);
    troughsArrayA = [];
    troughsTimeLocA = [];
    for j = 1:numPeaks
        
        peakTime = peaksTimeLocA(j);
        [M,I] = min(data(i, peakTime:min([peakTime+30, size(data(i,:),2)])));
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
            sigma = median(abs(data(i,startIndex:endIndex))/0.6745);
            thresholdB = 3.5*sigma;
            [p,t] = findpeaks(data(i,startIndex:endIndex),'MinPeakHeight',thresholdB);
            peaksArrayB =[peaksArrayB p];
            peaksTimeLocB = [peaksTimeLocB t+startIndex];
        end
        
    end
    troughsArrayB = [];
    troughsTimeLocB = [];
    for l = 1:size(peaksArrayB,2)
        peakTime = peaksTimeLocB(l);
        [M,I] = min(data(i, peakTime:min([peakTime+30, size(data(i,:),2)])));
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
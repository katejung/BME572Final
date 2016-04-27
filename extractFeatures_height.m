function [peaks, troughs, peakT, troughT] = extractFeatures_height(data)
arrayLength= size(data, 2);
numset = size(data, 1);
troughsArrayA = [];
troughsTimeLocA = [];

for i= 1:numset
    heightThreshold = 7;
    validPeak = [];
    validPeakLoc = [];
    validTrough =[];
    validTroughLoc = [];
    peaksArray = [];
    peaksTimeLoc = [];
    dataMedian = median(abs(data(i,:))/0.6745);
    threshold = 1.5*dataMedian;
    
    [peaksArray, peaksTimeLoc] = findpeaks(data(i,:),'MinPeakHeight',threshold);
    numPeaks = size(peaksArray,2);
    
    numPeaks = size(peaksArray,2);
    troughsArray = [];
    troughsTimeLoc = [];
    for j = 1:numPeaks
        peakTime = peaksTimeLoc(j);
        [M,I] = min(data(i, peakTime:min([peakTime+40, size(data(i,:),2)])));
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
            tempdata = data(i,:);
            time = validPeakLoc(a);
            for b = time:min([time+50 arrayLength])
                if tempdata(b) > 7*dataMedian
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


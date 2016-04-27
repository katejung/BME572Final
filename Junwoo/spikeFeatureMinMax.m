function [peaks,peakLoc, troughs, troughsLoc, peakinterval] = spikeFeatureMinMax(data,Fs)
arrayLength= size(data, 2);
numset = size(data, 1);
troughsArray;
for i= 1:numset
    sigma = median(abs(data(i,:))/0.6745);
    threshold = 4*sigma;
     [peaks{i,1},peakLoc{i,1}] = findpeaks(data(i,:),Fs,'MinPeakHeight',threshold);
     peakinterval{i,1} = diff(peakLoc{i,1}); 
    numPeaks = size(peaks{i,1});
    for j = 1:numPeaks
        location = peakLoc{i,1};
        [M,I] = min(data(i, location(j):location(j+10000)));
        troughArray = [troughArray M];
        troughsLoc = [troughsLoc I];
    end
end



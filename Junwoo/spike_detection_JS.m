function [output] = spike_detection_JS(input)
% Find spikes in filtered data as suggested in:
% http://www.scholarpedia.org/article/Spike_sorting
% The original article can be found in: reference/"Unsupervised Spike Detection and
% Sorting with Wavelets and Superparamagnetic Clustering".pdf
% Uses the cubic spines interpolatioin algorithm to avoid misalignment due
% to undersampling.

% For this particular project, input is a dataset of column vector (n by m)

datalength = size(input,2);
numset = size(input,1);
output = zeros(numset,datalength);


for i = 1:numset
    % Signma is an estimate of the standard deviation of background noise
     sigma = median(abs(input(i,:))/0.6745);
     % Reference suggests to use 4*sigma as a threshold.
     threshold = 4*sigma;
     
     
     
end




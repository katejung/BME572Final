function [output] = spike_detection_JS(data)
% Find spikes in filtered data as suggested in:
% http://www.scholarpedia.org/article/Spike_sorting
% The original article can be found in: reference/"Unsupervised Spike Detection and
% Sorting with Wavelets and Superparamagnetic Clustering".pdf

% For this particular project, input is a dataset of column vector (n by m)

datalength = size(data,2);
numset = size(data,1);
output = zeros(numset,datalength);

for i = 1:numset
    output(i,:) = filtfilt(d,data(i,:));
end




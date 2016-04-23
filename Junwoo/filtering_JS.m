function [filtered] = filtering_JS(data)
% Filter the data as suggested in:
% http://www.scholarpedia.org/article/Spike_sorting
% REQUIRE Matlab 2014a or later versions

% For this particular project, input is a dataset of column vector (n by m)
% Bandpass filter of 300 Hz to 5000Hz is used to remove noise and trend
d = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',300,'HalfPowerFrequency2',3000, ...
    'SampleRate',15000);

datalength = size(data,2);
numset = size(data,1);
filtered = zeros(numset,datalength);
samplerate = 15000; % Hz

for i = 1:numset
    filtered(i,:) = filtfilt(d,data(i,:));
    % Remove the DC-offset
    filtered(i,:) = filtered(i,:) - mean(filtered(i,:));
end




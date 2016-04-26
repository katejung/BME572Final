function [filtered] = filtering_JS(data)
% Filter the data as suggested in:
% http://www.scholarpedia.org/article/Spike_sorting
% REQUIRE Matlab 2014a or later versions

% For this particular project, input is a dataset of column vector (n by m)
% Bandpass filter of 300 Hz to 5000Hz is used to remove noise and trend

d = designfilt('bandpassfir', 'FilterOrder', 20, ...
    'CutoffFrequency1',50, 'CutoffFrequency2', 3000, ...
    'SampleRate', 15000);

% Bandstop filter around 60 is used to remove 60 Hz electronics noise
d1 = designfilt('bandstopfir', 'FilterOrder', 20, ...
    'CutoffFrequency1',50, 'CutoffFrequency2', 70, ...
    'SampleRate', 15000,'Window', 'hamming');

% Highpass filter of 3Hz is used to remove trend
d2 = designfilt('highpassfir', 'FilterOrder', 20, 'CutoffFrequency', ...
                3, 'SampleRate', 15000, 'Window', 'hamming');

datalength = size(data,2);
numset = size(data,1);
filtered = zeros(numset,datalength);
samplerate = 15000; % Hz

for i = 1:numset
    % Apply a frequency filter
    filtered(i,:) = filtfilt(d,data(i,:));
    filtered(i,:) = filtfilt(d1,filtered(i,:));
    filtered(i,:) = filtfilt(d2,filtered(i,:));
    
    % Remove additional trends with polynomial approximation
    [p,s,mu] = polyfit((1:datalength),filtered(i,:),10);
    trend = polyval(p,(1:datalength),[],mu);
    filtered(i,:) = filtered(i,:) - trend;
    
    % Remove the DC-offset
    filtered(i,:) = filtered(i,:) - mean(filtered(i,:));
    
end






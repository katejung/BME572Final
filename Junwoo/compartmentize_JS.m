function [baseline,intense,poststimulus] = compartmentize_JS(filtered,baselinec,intensec,poststimulusc)
% Filter the data as suggested in:
% http://www.scholarpedia.org/article/Spike_sorting

% For this particular project, input is a dataset of column vector (n by m)

datalength = size(data,2);
numset = size(data,1);
samplerate = 15000; % Hz

for i = 1:numset
    % First 10ms of the data is removed due to a large spike that does not
    % represent the data well.
    baseline(i,:) = filtered(i,1+samplerate*baselinec:samplerate*10); % Before an odor stimulus
    intense(i,:) = filtered(i,samplerate*10+1:samplerate*intensec); % Intense response
    poststimulus(i,:) = filtered(i,samplerate*intensec+1:end-poststimulusc); % After stimulus
end
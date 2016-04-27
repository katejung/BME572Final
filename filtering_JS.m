function [filtered] = filtering_JS(data, dealwithintense,isMAp)
% Filter the data as suggested in:
% http://www.scholarpedia.org/article/Spike_sorting
% REQUIRE Matlab 2014a or later versions

% For this particular project, input is a dataset of column vector (n by m)
% Bandpass filter of 300 Hz to 5000Hz is used to remove noise and trend

d = designfilt('bandpassfir', 'FilterOrder', 20, ...
    'CutoffFrequency1',50, 'CutoffFrequency2', 2000, ...
    'SampleRate', 15000);

% Bandstop filter around 60 is used to remove 60 Hz electronics noise
d1 = designfilt('bandstopfir', 'FilterOrder', 20, ...
    'CutoffFrequency1',50, 'CutoffFrequency2', 70, ...
    'SampleRate', 15000,'Window', 'hamming');

% Highpass filter of 3Hz is used to remove trend
d2 = designfilt('highpassfir', 'FilterOrder', 20, 'CutoffFrequency', ...
    3, 'SampleRate', 15000, 'Window', 'hamming');

% MAp_intense:
% 1: 1218 to 1326
% 2: 1263 to 1570
% 3: 2433 to 2516
% 4: 722 to 871
% 5: 1265 to 1439
% 6: 642 to 852
% 7: 1240 to 1412
% 8: 1279 to 1415
MAp_singularities = {[1153:1353],[1194:1418],[2391:2499], ...
    [689:824],[1219:1415],[656:826], ...
    [1199:1351],[1220:1397]};
MHp_singularities = {[4135:4347],[1356:1532],[1268:1374], ...
    [647:739],[783:924],[892:948], ...
    [1004:1114],[1002:1109]};

datalength = size(data,2);
numset = size(data,1);
filtered = zeros(numset,datalength);
samplerate = 15000; % Hz

for i = 1:numset
    if dealwithintense
        if isMAp
            sing = MAp_singularities{1,i};
        else
            sing = MHp_singularities{1,i};
        end
        % De-trend the datapoints prior to the singularity
        [p,s,mu] = polyfit(1:sing(1)-1,data(i,1:sing(1)-1),6);
        trend = polyval(p,1:sing(1)-1,[],mu);
        filtered(i,1:sing(1)-1) = data(i,1:sing(1)-1) - trend;
        
        % De-trend the datapoints in the singularity
        [p,s,mu] = polyfit(sing,data(i,sing),6);
        trend = polyval(p,sing,[],mu);
        filtered(i,sing) = data(i,sing) - trend;
        
        % De-trend the datapoints after the singularity
        [p,s,mu] = polyfit(sing(end)+1:datalength,data(i,sing(end)+1:end),6);
        trend = polyval(p,sing(end)+1:datalength,[],mu);
        filtered(i,sing(end)+1:end) = data(i,sing(end)+1:end) - trend;
        
    else
        % Remove additional trends with polynomial approximation
        [p,s,mu] = polyfit((1:datalength),data(i,:),10);
        trend = polyval(p,(1:datalength),[],mu);
        filtered(i,:) = data(i,:) - trend;
    end
    
    % Apply a frequency filter
    filtered(i,:) = filtfilt(d,filtered(i,:));
    filtered(i,:) = filtfilt(d1,filtered(i,:));
    filtered(i,:) = filtfilt(d2,filtered(i,:));
    
    % Remove the DC-offset
    filtered(i,:) = filtered(i,:) - mean(filtered(i,:));
    
end






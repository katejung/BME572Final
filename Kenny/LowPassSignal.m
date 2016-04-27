% clc;
% clear all;
% close all;

load MAp.mat
load MHp.mat

% training = [1 3 5 7];
training = [3];
test = [2 4 6 8];

samplerate = 15000;
t = 0:(1/samplerate):60; t = t(1:end-1);
baselinec = 0; % Remove first few seconds
intensec = 11.5; % Cutoff for the intense section
poststimulusc = 0; % Cutoff for the post-stimulus section

f = 15000;
threshold = 0;
start = 11.5;
m = MAp;

MApt = MAp(training(1),f*start+1:end);
for i = 1:size(training,2)-1
    MApt = horzcat(MApt,MAp(training(i+1),f*start+1:end));
end
% Butterworth bandpass filter
[MAp_baseline,MAp_intense,MAp_poststimulus] = compartmentize_JS(MAp,baselinec,intensec,poststimulusc);
% [MAp_baseline] = filtering_JS(MAp_baseline,0,1); [MAp_intense] = filtering_JS(MAp_intense,1,1); [MAp_poststimulus] = filtering_JS(MAp_poststimulus,0,1);
% MAp_filtered = [MAp_baseline MAp_intense MAp_poststimulus];
% [MHp_baseline,MHp_intense,MHp_poststimulus] = compartmentize_JS(MHp,baselinec,intensec,poststimulusc);
% [MHp_baseline] = filtering_JS(MHp_baseline,0,0); [MHp_intense] = filtering_JS(MHp_intense,1,0); [MHp_poststimulus] = filtering_JS(MHp_poststimulus,0,0);
% MHp_filtered = [MHp_baseline MHp_intense MHp_poststimulus];
%%
data = MAp_poststimulus;
dealwithintense =0;
isMAp = 1;
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
MAp_poststimulus1 = filtered;
%%
% MAp_poststimulus2 = nan(8,727500,49);
% for j = 1:49
    data = MAp_poststimulus;
    dealwithintense =0;
    isMAp = 1;
    d = designfilt('bandpassfir', 'FilterOrder', 20, ...
        'CutoffFrequency1',15, 'CutoffFrequency2', 35, ...
        'SampleRate', 15000);
    % d = designfilt('bandpassfir', 'FilterOrder', 20, ...
    %     'CutoffFrequency1',50, 'CutoffFrequency2', 3000, ...
    %     'SampleRate', 15000);
    
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
%     MAp_poststimulus2(:,:,j) = filtered;
MAp_poststimulus3 = filtered;
% end

%%
MAplp = MAp_poststimulus1(1,:);

pks = zeros(1,size(MAplp,2));
locs = zeros(1,size(MAplp,2));
% (max(MAplp(i-x1:i+x1)) - min(MAplp(i-x1:i+x1))) >= 20

%  && (max(MAplp(i-x1:i+x1)) - min(MAplp(i-x1:i+x1)) > 10)
x1 = 16; % half the window size. Window size = 2x+1
for i = 1+x1:size(MAplp,2)-x1
    if (max(MAplp(i-x1:i+x1)) == MAplp(i))
        pks(i) = MAplp(i);
        locs(i) = i;
    end
end
pks = pks(pks~=0);
locs = locs(locs~=0);
thrclass = pks>0;
pkst = pks.*(thrclass);
pkst = pkst(pkst~=0);
locst = locs.*(thrclass);
locst = locst(locst~=0);
figure;
plot(linspace(1./f,size(MAplp,2)./f,size(MAplp,2)),MAplp);
ylim([-200 200]);
hold on;
plot(locst./f,pkst,'.');
hold off;

%% Histogram idea to separate oscillations from others (possible use)
% histogramsize = 1;
% amp = linspace(min(MAplp),max(MAplp),round((max(MAplp)-min(MAplp))/histogramsize));
% count = zeros(1,size(amp,2));
% for i = 1:size(amp,2)-1
%     count(i) = size(MAplp((MAplp>amp(i)) == (MAplp<amp(i+1))),2);
% end
% figure;
% plot(amp,count)

%%
x2 = 16; % Size of window of data
data = zeros(size(locst,2),2*x2+1);
for i = 1:size(locst,2)
    data(i,:) = MAplp((locst(i) - x2):(locst(i) + x2));
end

[COEFF, SCORE, LATENT] = pca(data);
figure;
plot(LATENT(1:x2+1));

pcadata = data*COEFF;
figure;
scatter(pcadata(:,1),pcadata(:,2));
% 
% figure;
% scatter3(pcadata(:,1),pcadata(:,2),pcadata(:,3));

%%
figure;
hold on;
% for i = 1:size(data,1)
for i = 1:size(data,1)
plot(data(i,:),'b');
end
hold off;

%%
plot(C(:,1),'r');
hold on;
plot(MAp_baseline(1,:),'b');
hold off;

%%
plot(C(:,1),C(:,2),'.')
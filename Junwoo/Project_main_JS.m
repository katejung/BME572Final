% Project_main_JS.m
clear all; clc;
% Load data
load MAp.mat;
load MHp.mat;

%% Plotting settings
linespec = {'linestyle' '-';'linestyle' ':';'linestyle' '-.';'linestyle' '--'}; % Line specs for gray-scale
markerspec = {'linestyle' 'o';'linestyle' '*';'linestyle' '+';'linestyle' '.'}; % Marker specs for gray-scale
color = {'color' 'blue'; 'color' 'red'; 'color' 'green'; ...
    'color' 'cyan'; 'color' 'black'; 'color' 'yellow'};
% fontsize{1,:} for title, fontsize{2,:} for axis, fontsize{3,:} for legend
fontsize = {'Fontsize' 13 'Fontweight' 'bold'; ...
    'Fontsize', 10, 'Fontweight', 'bold'; ...
    'Fontsize', 11, 'Fontweight', 'bold'};
% linespec2{1,:} for line, linespec2{2,:} for marker
linespec2 = {'linewidth' 2; 'MarkerSize' 6};
% Count the number of figures created
fignum = 1;

%% Preprosess the data
% preprocess parameters
samplerate = 15000;
t = 0:(1/samplerate):60; t = t(1:end-1);
baselinec = 0; % Remove first few seconds
intensec = 11.5; % Cutoff for the intense section
poststimulusc = 0; % Cutoff for the post-stimulus section

% Preprocess the data
[MAp_baseline,MAp_intense,MAp_poststimulus] = compartmentize_JS(MAp,baselinec,intensec,poststimulusc);
[MAp_baseline] = filtering_JS(MAp_baseline,0,1); [MAp_intense] = filtering_JS(MAp_intense,1,1); [MAp_poststimulus] = filtering_JS(MAp_poststimulus,0,1);
MAp_filtered = [MAp_baseline MAp_intense MAp_poststimulus];
[MHp_baseline,MHp_intense,MHp_poststimulus] = compartmentize_JS(MHp,baselinec,intensec,poststimulusc);
[MHp_baseline] = filtering_JS(MHp_baseline,0,0); [MHp_intense] = filtering_JS(MHp_intense,1,0); [MHp_poststimulus] = filtering_JS(MHp_poststimulus,0,0);
MHp_filtered = [MHp_baseline MHp_intense MHp_poststimulus];

% Visualize
fig = figure(fignum); clf(fig)
plot(t(1+samplerate*baselinec:end),MAp_filtered,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MAp filtered data',fontsize{1,:})
fignum = fignum + 1;

fig = figure(fignum); clf(fig)
subplot(3,1,1)
plot(t(1+samplerate*baselinec:samplerate*10),MAp_baseline,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MAp filtered data: Baseline',fontsize{1,:})
subplot(3,1,2)
plot(t(samplerate*10+1:samplerate*intensec),MAp_intense,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MAp filtered data: Intense',fontsize{1,:})
subplot(3,1,3)
plot(t(samplerate*intensec+1:end-poststimulusc),MAp_poststimulus,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MAp filtered data: Post-stimulus',fontsize{1,:})
fignum = fignum + 1;

fig = figure(fignum); clf(fig)
plot(t(1+samplerate*baselinec:end),MHp_filtered,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MHp filtered data',fontsize{1,:})
fignum = fignum + 1;

fig = figure(fignum); clf(fig)
subplot(3,1,1)
plot(t(1+samplerate*baselinec:samplerate*10),MHp_baseline,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MHp filtered data: Baseline',fontsize{1,:})
subplot(3,1,2)
plot(t(samplerate*10+1:samplerate*intensec),MHp_intense,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MHp filtered data: Intense',fontsize{1,:})
subplot(3,1,3)
plot(t(samplerate*intensec+1:end-poststimulusc),MHp_poststimulus,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MHp filtered data: Post-stimulus',fontsize{1,:})
fignum = fignum + 1;

%% spike detection
[pks,locs,peakinterval] = spike_detection_JS(MAp_baseline);
fig = figure(fignum); clf(fig)
plot(t(1+samplerate*baselinec:samplerate*10),MAp_baseline(1,:),linespec{1,:},linespec2{1,:});
hold on
plot(t(locs{1,1}+samplerate*baselinec),MAp_baseline(1,locs{1,1}),'rs','MarkerFaceColor','r');
axis([2.03 2.05 -25 50])
fignum = fignum  +1;

fig = figure(fignum); clf(fig)
plot(t,MAp(1,:),linespec{1,:},linespec2{1,:});
axis([2.03 2.05 -65 15])
fignum = fignum  +1;

[pks,locs,peakinterval] = spike_detection_JS(MAp_poststimulus);
fig = figure(fignum); clf(fig)
plot(t(samplerate*intensec+1:end-poststimulusc),MAp_poststimulus(1,:),linespec{1,:},linespec2{1,:});
hold on
plot(t(locs{1,1}+samplerate*intensec+1),MAp_poststimulus(1,locs{1,1}),'rs','MarkerFaceColor','r');
% axis([14.98 15 -25 50])
fignum = fignum  +1;

fig = figure(fignum); clf(fig)
plot(t,MAp(1,:),linespec{1,:},linespec2{1,:});
% axis([14.98 15 -370 -340])
fignum = fignum  +1;

[pks,locs,peakinterval] = spike_detection_JS(MHp_poststimulus);
fig = figure(fignum); clf(fig)
plot(t(samplerate*intensec+1:end-poststimulusc),MHp_poststimulus(1,:),linespec{1,:},linespec2{1,:});
hold on
plot(t(locs{1,1}+samplerate*intensec+1),MHp_poststimulus(1,locs{1,1}),'rs','MarkerFaceColor','r');
% axis([14.98 15 -25 50])
fignum = fignum  +1;

fig = figure(fignum); clf(fig)
plot(t,MHp(1,:),linespec{1,:},linespec2{1,:});
% axis([14.98 15 -370 -340])
fignum = fignum  +1;















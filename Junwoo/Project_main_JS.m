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

%% Preprosess the data and find spikes

% preprocess parameters
samplerate = 15000;
t = 0:(1/samplerate):60; t = t(1:end-1);
baselinec = 0.01; % Remove first 10ms to remove the unrealistically large spike
intensec = 10.6; % Cutoff for the intense section
poststimulusc = 0; % Cutoff for the post-stimulus section

% Preprocess the data
[MAp_filtered] = filtering_JS(MAp);
[MAp_baseline,MAp_intense,MAp_poststimulus] = compartmentize_JS(MAp_filtered,baselinec,intensec,poststimulusc);
[MHp_filtered] = filtering_JS(MHp);
[MHp_baseline,MHp_intense,MHp_poststimulus] = compartmentize_JS(MHp_filtered,baselinec,intensec,poststimulusc);

% Visualize
fig = figure(fignum); clf(fig)
plot(t,MAp_filtered,linespec{1,:},linespec2{1,:})
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
plot(t(samplerate*10.6+1:end-poststimulusc),MAp_poststimulus,linespec{1,:},linespec2{1,:})
xlabel('Time (s)', fontsize{2,:}); ylabel('Amplitude (arb. Units)', fontsize{2,:})
title('MAp filtered data: Post-stimulus',fontsize{1,:})
fignum = fignum + 1;

fig = figure(fignum); clf(fig)
plot(t,MHp_filtered,linespec{1,:},linespec2{1,:})
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

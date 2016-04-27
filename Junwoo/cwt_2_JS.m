load MAp.mat;
load MHp.mat;

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
%%
FS = 15000;
F = MAp_baseline(1,1751:1879);

% Build and display the adapted wavelet.
% [Y,X] = pat2cwav(F,'polynomial',5,'continuous');
[Y,X] = pat2cwav(F,'orthconst',0,2);
figure;
plot(X,F/max(abs(F)),'b');
hold on;
plot(X,Y/max(abs(Y)),'r');
title('Form to detect (b) and adapted Wavelet (r)')

% Save the adapted wavelet and add it to the toolbox
locdir = cd;

cd(tempdir);
save adp_FRM2 X Y
wavemngr('add','AdapF2','adpf2',4,'','adp_FRM2.mat',[0 1]);
addpath(tempdir,'-begin');
cd(locdir);

% time = linspace(0,1,length(Z));
% figure;
% plot(time,Z); grid on;

a = figure(5); waveletstuff = VideoWriter('Propagation2','MPEG-4');
open(waveletstuff);
% Define the form to detect.
for jj = 1:10:3000
    
% Construct the signal containing two similar forms.
Z = MAp_baseline(1,jj:FS*0.1+jj);

% We analyze the signal by computing the CWT coefficients of Z using the
% admissible wavelet we constructed to approximate the basic form F.

FS = 15000;

stepSIG = 1/FS;
stepWAV = 1/128;
wname = 'adpf1';
scales  = (1:128)*stepSIG;
WAV = {wname,stepWAV};
SIG = {Z,stepSIG};

c = cwt(SIG,scales,WAV,'scalCNT'); grid
mov = getframe(a);
writeVideo(waveletstuff,mov);

% Detect the maximum of the absolute value of coefficients.
lenSIG = length(Z);
positions = (0:lenSIG-1)*stepSIG;
% fprintf('Instant 1:  %4.2f\n',positions(round(middle)))
% fprintf('Instant 2:  %4.2f\n',positions(round(middle2)))
end
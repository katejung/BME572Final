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

%% Continuous wavelet ananlysis

% Length of the signal and half-length of the pattern.
lenSIG = 897;
L = 45;

% Length, beginning, end and middle of the pattern.
long  = 2*L;
first = 340;
last  = first+long-1;
pat = [-ones(1,L) ones(1,L)];

% Signal with inserted shape.
Z = zeros(1,lenSIG);
Z (first:last) = pat;
scales = 1:256;
wname = 'haar';

% Continuous wavelet transform (CWT).
figure;
c = cwt(Z,scales,wname,'scalCNT'); grid

% Detect the maximum of the absolute value of coefficients.
[~,imax] = max(abs(c(:)));
SIZ = size(c);
ROW = rem(imax,SIZ(1));
if ROW==0 , ROW = SIZ(1); end
[~,COL] = max(abs(c(ROW,:)));
display(sprintf('Detected indices COL and ROW: %3.0f %3.0f',COL,ROW))

% Set the sampling period.
stepSIG = 0.025;

scales  = (1:0.25:10);
positions = (0:lenSIG-1)*stepSIG;

% Sampling rate of the analyzing wavelet.
stepWAV = 1/1024;
wname   = 'haar';
WAV = {wname,stepWAV};

% Caution, SIG is a cell array containing the signal Z and the sampling
% period stepSIG.
SIG = {Z,stepSIG};

% Recall that scales are expressed as multiple of the sampling period.
figure;
c = cwt(SIG,scales,WAV,'scalCNT'); grid
ylabel('Duration')

% Detect the maximum of the absolute value of coefficients.
[~,imax] = max(abs(c(:)));
SIZ = size(c);
ROW = rem(imax,SIZ(1));
if ROW==0 , ROW = SIZ(1); end
[~,COL] = max(abs(c(ROW,:)));
display(sprintf('Detected indices COL and ROW: %3.0f %3.0f',COL,ROW))
display(sprintf('Duration: %4.2f',scales(ROW)))
display(sprintf('Instant:  %4.2f',positions(COL)))

%%
a = figure(4); waveletstuff = VideoWriter('Propagation','MPEG-4');
open(waveletstuff);
% Define the form to detect.
for jj = 1:300
FS = 15000;
F = MAp_poststimulus(1,3053:3181);

% Construct the signal containing two similar forms.
Z = MAp_poststimulus(1,jj:FS*0.1+jj);

% Build and display the adapted wavelet.
% [Y,X] = pat2cwav(F,'polynomial',5,'continuous');
[Y,X] = pat2cwav(F,'orthconst',0,2);
% figure;
% plot(X,F/max(abs(F)),'b');
% hold on;
% plot(X,Y/max(abs(Y)),'r');
% title('Form to detect (b) and adapted Wavelet (r)')

% Save the adapted wavelet and add it to the toolbox
locdir = cd;

cd(tempdir);
save adp_FRM1 X Y
wavemngr('add','AdapF1','adpf1',4,'','adp_FRM1.mat',[0 1]);
addpath(tempdir,'-begin');
cd(locdir);

% time = linspace(0,1,length(Z));
% figure;
% plot(time,Z); grid on;

% We analyze the signal by computing the CWT coefficients of Z using the
% admissible wavelet we constructed to approximate the basic form F.

FS = 15000;

stepSIG = 1/FS;
stepWAV = 1/256;
wname = 'adpf1';
scales  = (1:256)*stepSIG;
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













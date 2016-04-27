% clc;
% clear all;
% close all;

load MAp.mat
load MHp.mat

% training = [1 3 5 7];
training = [3];
test = [2 4 6 8];

clc;
close all;

f = 15000;
threshold = 0;
start = 10;
m = MAp;

MApt = MHp(training(1),1:end);
for i = 1:size(training,2)-1
    MApt = horzcat(MApt,MAp(training(i+1),1:end));
end
% Butterworth bandpass filter
% d = designfilt('bandpassiir','FilterOrder',20, ...
%     'HalfPowerFrequency1',300,'HalfPowerFrequency2',3000, ...
%     'SampleRate',15000);
% MAplp = filter(d,MApt);
[b, a] = butter(2,[15/7500 35/7500],'bandpass');
MAplp1 = filter(b,a,MApt);

x1 = 32; % half the window size. Window size = 2x+1

figure;
plot(linspace(1./f,size(MAplp1,2)./f,size(MAplp1,2)),MAplp1);
% ylim([-200 200]);

%%
figure
plot(linspace(1./f,size(MAplp,2)./f,size(MAplp,2)),MAplp);
hold on
plot(linspace(1./f,size(MAplp1(1:150000),2)./f,size(MAplp1(1:150000),2)),MAplp1(1:150000),'r');
hold off

%%
figure;
plot(linspace(1./f,size(MAplp1,2)./f,size(MAplp1,2)),MAplp1,'r');


%% Find peaks using Kate's code + modification
f = 15000;
dataInterest = MAplp1(f*11.5:end);
% [validPeakCell, validTroughCell, validPeakLocCell, validTroughLocCell] = extractFeatures_height(dataInterest);
[peaks, troughs, peakT, troughT] = extractFeatures_height2(dataInterest);
num = size(peaks,2);
for j = 1:num
    heights(j) = peaks(j)-troughs(j);
    widths(j) = troughT(j)-peakT(j);
end

dataMedian = 0;
for m = 1:1
    dataMedian = dataMedian + median(abs(dataInterest(m,:))/0.6745);
end
dataMedian = dataMedian/4;
window = 40;
data = nan(size(peakT,2),window);
for i = 1:size(data,1)
    data(i,:) = dataInterest((peakT - 19):(peakT + 20));
end
[COEFFICIENT, SCORE, LATENT] = pca(data);
plot(data*COEFFICIENT(:,1),data*COEFFICIENT(:,2));
% K =2;% number of principal Comonents
% 
% X(:,1) = troughs';
% X(:,2) = peaks';
% [n,d]=size(X);
% C=cov(X);
% [V U]=eig(C); %V = eigenvectors U = eigenvalues
% L=diag(U); L', U % principal component directions
% [sortedU, index] = sort(L, 'descend');
% Xproj_PCA = zeros(d,K);
% for j = 1:K
%     Xproj_PCA(:,j) = V(:,index(j));
% end
% Y_trained_PCA = X*Xproj_PCA;
% Xproj_PCA; % the first column is the 1st PC and the second column is the 2nd PC
% 
% % calulate 1st PC line
% Xh1_PCA=-35:0.5:5;
% Yh1_PCA = Xh1_PCA.*Xproj_PCA(2,1)./Xproj_PCA(1,1);
% % calculated 2nd PC line
% Xh2_PCA =-35:0.5:5;
% Yh2_PCA = Xh2_PCA.*Xproj_PCA(2,2)./Xproj_PCA(1,2);
% figure
% plot(X(:,1), X(:,2), 'b^', Xh1_PCA,Yh1_PCA,'g-',Xh2_PCA, Yh2_PCA, 'r-')
% legend('data', '1st PCA axis', '2nd PCA axis')
% figure
% 
% scatter(Y_trained_PCA(:,1), Y_trained_PCA(:,2), 'bo')
% 
% title('Encoded Training Data')
% xlabel(['1st PC [' num2str(Xproj_PCA(1,1)) ',' num2str(Xproj_PCA(2,1)) ']'])
% ylabel(['2nd PC [' num2str(Xproj_PCA(1,2)) ',' num2str(Xproj_PCA(2,2)) ']'])
clear
clc
figNum = 1;

MAp = load('MAp.mat');
MHp = load('MHp.mat');

for i = 1:4
    train(i,:) = MAp.MAp(i,:);
end

figure (figNum)
figNum = figNum+1;

for i = 1:4
    subplot(2,2,i)
    plot(train(i,:))
     title('raw data')
end
%%
data = train(i,:);
samples = size(data,2);
dt = 60/samples;
t = dt:dt:60;

plot(t,data)
xlim([0,10])

baselineT = dt:dt:10;
baselineData = data(1:size(baselineT,2));



Fs = 15000;

fc = 2000;
Wn = (2/Fs)*fc;
b = fir1(20,Wn,'low',kaiser(21,3));

fvtool(b,1,'Fs',Fs)
y = filter(b,1,baselineData);
figure(figNum)
figNum= figNum+1;
plot(baselineData)
hold on
plot(y,'r')
xlim([0,20000])
xlabel('Time (s)')
ylabel('Amplitude')
legend('Original Signal','Filtered Data')
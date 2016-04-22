clear
clc
figNum = 1;

MAp = load('MAp.mat');
MHp = load('MHp.mat');

for i = 1:4
    train_MAp(i,:) = MAp.MAp(i,:);
    train_MHp(i,:) = MHp.MHp(i,:);
end

figure (figNum)
figNum = figNum+1;

for i = 1:4
    subplot(2,2,i)
    plot(train_MAp(i,:))
     title('raw data')
end
%%
data_MAp = train_MAp(i,:);
data_MHp = train_MHp(i,:);
samples = size(data_MAp,2);
dt = 60/samples;
t = dt:dt:60;

figure (figNum)
figNum = figNum+1;
plot(t,data_MAp)
xlim([0,10])

baselineT = dt:dt:10;
baselineData_MAp = data_MAp(1:size(baselineT,2));
baselineData_MHp = data_MHp(1:size(baselineT,2));


Fs = 15000;

fc = 2000;
Wn = (2/Fs)*fc;
b = fir1(20,Wn,'low',kaiser(21,3));

fvtool(b,1,'Fs',Fs)
y_MAp = filter(b,1,baselineData_MAp);
y_MHp = filter(b,1,baselineData_MHp);
figure(figNum)
figNum= figNum+1;
plot(baselineData_MAp)
plot(baselineData_MHp)
hold on
plot(y_MAp,'r')
xlim([0,20000])
xlabel('Time (s)')
ylabel('Amplitude')
legend('Original Signal','Filtered Data')

%%
Y_MAp = fft(y_MAp);
Y_MHp = fft(y_MHp);
L = length(Y_MAp);

P2_MAp = abs(Y_MAp/L);
P2_MHp = abs(Y_MHp/L);

P1_yMAp = P2_MAp(1:L/2+1);
P1_yMAp(2:end-1) = 2*P1_yMAp(2:end-1);
P1_yMHp = P2_MHp(1:L/2+1);
P1_yMHp(2:end-1) = 2*P1_yMHp(2:end-1);


f = Fs*(0:(L/2))/L;
plot(f, P1_yMHp, 'r') 
hold on
plot(f,P1_yMAp)

title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0,2000])
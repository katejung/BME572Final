% test.m

% Wavelet decomposion
% C = zeros(size(MAp_baseline,2)-63,64);
% for i = 1:size(MAp_baseline,2)-1
%     [C(i,:) L] = wavedec(MAp_baseline(1,i:i+63),4,'haar');
% end
% 
% figure;
% plot(MAp_baseline(1,:))
% hold on
% plot(peakT,peaks,'.')

figure(1);clf(1);plot(t(samplerate*intensec+1:end),MAp_poststimulus(1,:))
[upper, lower] = envelope(MAp_poststimulus(1,:),1500,'peak');
hold on; figure(1);plot(t(samplerate*intensec+1:end),upper(1,:));

figure(2);plot(t(samplerate*intensec+1:samplerate*29),MAp_poststimulus(1,1:samplerate*(29-intensec)));
hist(MAp_poststimulus(1,1:samplerate*29),100)

X = 1:samplerate*(29-intensec); Y = sin(X/1500); 
figure(3); hist(Y);

[wn,yn] = ica_self([MAp_poststimulus(1,:)',MHp_poststimulus(1,:)']);

subplot(2,2,1)
plot(MAp_baseline(1,:))
subplot(2,2,2)
plot(MHp_baseline(1,:))
subplot(2,2,3)
plot(yn(1,:))
subplot(2,2,4)
plot(yn(2,:))

% figure(4); plot(t(samplerate*intensec+1:end),MAp_poststimulus(1,:));
% hold on; plot(peakT,peaks);
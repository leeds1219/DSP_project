load lead2.mat
% Notch filter design
f0 = 60;
fs = 500;
r = 0.95;
w = 2*pi*f0/fs;
j = sqrt(-1);
bb = [1 -(exp(j*w)+exp(-j*w)) 1];
aa = [1 -r*(exp(j*w)+exp(-j*w)) r*r];

% Apply notch filter
yy = filter(bb,aa, lead2);
N = length(lead2);
t = (0:N-1)/fs;

% Plot in time domain
figure;
subplot(2,1,1);
plot(t,lead2);
title('Input Signal');
xlabel('Time(s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t,yy);
title('Filtered Signal')
xlabel('Time(s)');
ylabel('Amplitude');

% Plot in frequency domain
HH = fft(yy);
hh = fft(lead2);
f = (0:N-1)*(fs/N);

figure;
subplot(2,1,1);
plot(f,abs(hh));
title('Original');
xlabel('Hz');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(f, abs(HH));
title('Filtered')
xlabel('Hz')
ylabel('Magnitude')
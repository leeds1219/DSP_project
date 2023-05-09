[xx,fs] = audioread('SunshineSquare.wav');
xx = xx';
n = 0:length(xx)-1;


% Plot the real and imaginary parts separately
subplot(2,1,1)
plot(n, real(xx))
title('Real part of x[n]')
xlabel('n')
ylabel('Amplitude')
subplot(2,1,2)
plot(n, imag(xx))
title('Imaginary part of x[n]')
xlabel('n')
ylabel('Amplitude')

% Plot the magnitude and phase
subplot(2,1,1)
plot(n, abs(xx))
title('Magnitude of x[n]')
xlabel('n')
ylabel('Amplitude')
subplot(2,1,2)
plot(n, angle(xx))
title('Phase of x[n]')
xlabel('n')
ylabel('Radians')
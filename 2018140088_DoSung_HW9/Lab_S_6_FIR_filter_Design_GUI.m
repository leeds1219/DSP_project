% filterdesign
% file export를 통해서 꺼낼 수 있다.
% rectangle
% hamming
%ww = -pi:2*pi/10000:pi;
%b = freqz(num,1,ww);
%f1=figure;
%plot(ww,abs(b));
% find the stop/pass band freq
% rectangle
%b_rec=freqz(num_rectangle,1,ww);
%f2=figure;
%plot(ww,abs(b_rec));

% Define passband and stopband frequencies and ripple/attenuation
%w_p = 0.68*pi;
%w_s = 0.72*pi;
%d_p = 0.01;
%d_s = 0.01;

% Define ideal filter response
%n = 1000; % number of frequency samples
%w = linspace(0, pi, n); % frequency axis (radians)
%h_ideal = zeros(1, n); % ideal magnitude response
%pb_range = (w >= 0) & (w <= w_p); % passband range
%sb_range = (w >= w_s) & (w <= pi); % stopband range
%h_ideal(pb_range) = 1; % passband magnitude
%h_ideal(sb_range) = 0; % stopband magnitude

% Plot ideal filter response
%figure;
%plot(w, h_ideal);
%xlabel('Frequency (radians)');
%ylabel('Magnitude');
%title('Ideal Lowpass Filter Response');
%axis([0 pi 0 1.1]); % set axis limits
%grid on;


% 3
w_p = 0.68*pi; % in rad/sample
w_s = 0.72*pi; % in rad/sample
d_p = 0.05;
d_s = 0.01;
fs = 2; % no specified fs... 
f = [w_p/(2*pi), w_s/(2*pi)];
a = [1, 0];
dev = [d_p, d_s];

% (a) filter design
% Calculate filter order
[n,fo,ao,w] = firpmord(f,a,dev,fs);
% n = n + 1;
n = n + 2;
% Design filter coefficients
b = firpm(n,fo,ao,w);

% (b) plot impulse response 
f1 = figure;
impz(b);
% Recall that the filter coefficients of the FIR filter are
% the values of the impulse response. Also, the DTFT of the impulse response is the frequency response
% of the FIR filter.

% (c) Make a plot of the frequency response magnitude versus omega hat
% plot frequency response magnitude
[H,w] = freqz(b,1);

% Plot frequency response magnitude vs omega hat
f2=figure;
plot(w/pi, abs(H));
xlabel('Omega hat');
ylabel('Magnitude');

% Calculate slope of linear portion of plot
L = n+3; % filter order + 1 + 2 = 149
width = w_s-w_p; % 0.04*pi rad/sample
C = L*width; % cutoff frequency
fc=C*fs/(2*pi);
slope = -1/(fc*2); % -0.0828dB/Hz
%lin_range = w/pi >= 0.68 & w/pi <= 0.72;
%p = polyfit(w(lin_range)/pi, abs(H(lin_range)), 1);
%slope = -p(1); % -0.0709dB/Hz


load 'ecg_withnoise.mat'
% Define sampling frequency and duration
fs = 500; % Sampling frequency (Hz)
duration = 10; % Signal duration (seconds)

% Generate time vector
t = 0:1/fs:duration-1/fs; % Time vector

% Compute the FFT of the ECG signal
fft_ecg = fft(ecg);

% Create frequency vector corresponding to the FFT
L = length(ecg); % Length of the signal
f = fs*(0:(L/2))/L; % Frequency vector

% Define the desired frequency range for the BPF
f_low = 0.5; % Lower cutoff frequency (Hz)
f_high = 50; % Upper cutoff frequency (Hz)

% Create logical mask to zero out unwanted frequencies
mask = (f >= f_low) & (f <= f_high);

% Apply the mask to the FFT
filtered_fft_ecg = fft_ecg;
filtered_fft_ecg(~mask) = 0;

% Reconstruct the filtered ECG signal using inverse FFT
filtered_ecg = ifft(filtered_fft_ecg);

% Plot original and filtered ECG signals in time domain
figure;
subplot(2,1,1);
plot(t, ecg);
title('Original ECG Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, real(filtered_ecg)); % Use real() to remove small imaginary components
title('Filtered ECG Signal');
xlabel('Time (s)');
ylabel('Amplitude')

% Plot single-sided magnitude spectrum of ECG data in subplots
figure;

for n = 0:1:4
    t_start=2*n;
    t_end=2*(n+1);
    n_start = round(t_start * fs) + 1; % Start sample index (+1 to account for 1-based indexing)
    n_end = round(t_end * fs); % End sample index
    ecg_segment = filtered_ecg(n_start:n_end);

    N = numel(ecg_segment); % Length of the ECG segment
    ecg_spectrum = fft(ecg_segment);

    ecg_spectrum_mag = abs(ecg_spectrum(1:N/2+1)); % Take only the non-negative frequencies
    ecg_spectrum_mag(2:end-1) = 2 * ecg_spectrum_mag(2:end-1); % Double the amplitudes (except for DC and Nyquist frequencies)

    f = (0:N/2) * fs / N; % Frequency axis in Hz

    subplot(5,1,n+1);
    plot(f, ecg_spectrum_mag);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title('Single-Sided Magnitude Spectrum of ECG Data');
end

% Define sampling frequency and duration
fs = 500; % Sampling frequency (Hz)
duration = 10; % Signal duration (seconds)

% Generate time vector
t = 0:1/fs:duration-1/fs; % Time vector

% Design parameters for the band-reject filter
f_stop_low = 6.8; % Lower stopband frequency (Hz)
f_stop_high = 12; % Upper stopband frequency (Hz)

% Normalize frequencies
Ws = [f_stop_low, f_stop_high] / (fs/2);

% Design the band-reject filter using FIR filter design
order = 100; % Filter order
b = fir1(order, [Ws(1) Ws(2)], 'stop');

% Apply the band-reject filter to the input signal
filtered_ecg_filtered = filtfilt(b, 1, filtered_ecg);

% Plot the original and filtered ECG signals
figure;
subplot(2,1,1);
plot(t, filtered_ecg);
title('Original Filtered ECG Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, filtered_ecg_filtered);
title('Filtered ECG Signal');
xlabel('Time (s)');
ylabel('Amplitude');

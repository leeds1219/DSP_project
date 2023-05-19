load 'ecg 12ch.txt' 
fs = 500; % Sampling frequency (Hz)
t_start = 1; % Start time (seconds)
t_end = 2; % End time (seconds)
n_start = round(t_start * fs); % Start sample index
n_end = round(t_end * fs); % End sample index
ecg_data = ecg_12ch(:,2);
ecg_segment = ecg_data(n_start:n_end);

filter_length = 16;
averaging_filter = ones(filter_length, 1) / filter_length;

filtered_ecg_no_padding = conv(ecg_segment, averaging_filter, 'same');

padding_length = filter_length - 1;
padded_ecg_segment = [zeros(floor(padding_length/2), 1); ecg_segment; zeros(ceil(padding_length/2), 1)];
filtered_ecg_with_padding = conv(padded_ecg_segment, averaging_filter, 'same');

f1 = figure;
plot(filtered_ecg_no_padding);
hold on
plot(filtered_ecg_with_padding);

xlabel('Sample');
ylabel('Amplitude');
title('Comparison of Filtered ECG Results');
legend('No Zero-Padding', 'With Zero-Padding');

hold off;

% 6
t_start=0;
t_end=10;
n_start = round(t_start * fs) + 1; % Start sample index (+1 to account for 1-based indexing)
n_end = round(t_end * fs); % End sample index
ecg_segment = ecg_12ch(n_start:n_end, 2);

N = numel(ecg_segment); % Length of the ECG segment
ecg_spectrum = fft(ecg_segment);

ecg_spectrum_mag = abs(ecg_spectrum(1:N/2+1)); % Take only the non-negative frequencies
ecg_spectrum_mag(2:end-1) = 2 * ecg_spectrum_mag(2:end-1); % Double the amplitudes (except for DC and Nyquist frequencies)

f = (0:N/2) * fs / N; % Frequency axis in Hz

f2 = figure;
plot(f, ecg_spectrum_mag);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Single-Sided Magnitude Spectrum of ECG Data');

% 7 LPF
cutoff_freq = 25; % Cutoff frequency of the LPF in Hz
N = numel(ecg_segment); % Length of the ECG segment

frequency_response = zeros(1, N);
frequency_response(1:round(cutoff_freq*N/fs)+1) = 1;
frequency_response(end-round(cutoff_freq*N/fs)+1:end) = 1;

filtered_ecg_spectrum = ecg_spectrum .* frequency_response.';

filtered_ecg_segment = ifft(filtered_ecg_spectrum,'symmetric');

t = (0:N-1) / fs; % Time axis in seconds
f3 = figure;
plot(t, real(filtered_ecg_segment));
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered ECG Segment (Ideal LPF)');

% 8 HPF
cutoff_freq = 0.5; % Cutoff frequency of the HPF in Hz
N = numel(ecg_segment); % Length of the ECG segment

frequency_response = ones(1, N);
frequency_response(1:round(cutoff_freq*N/fs)+1) = 0;
frequency_response(end-round(cutoff_freq*N/fs)+1:end) = 0;

filtered_ecg_spectrum = ecg_spectrum .* frequency_response.';

filtered_ecg_segment = ifft(filtered_ecg_spectrum,'symmetric');

f4 =figure;
t = (0:N-1) / fs; % Time axis in seconds
plot(t, real(filtered_ecg_segment));
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered ECG Segment (Ideal HPF)');


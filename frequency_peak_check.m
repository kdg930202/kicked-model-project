clearvars
% Parameters
A = 1;                       % Amplitude
f = 10;                      % Frequency in Hz (for fixed frequency case)
fs = 1000;                   % Sampling frequency in Hz
T = 1;                       % Duration in seconds
t = 0:1/fs:T-1/fs;           % Time vector

% Cosine signal
y = A * cos(2*pi*f*t);

% FFT
n = length(y);
Y = fft(y);
f_axis = (0:n-1)*(fs/n);     % Frequency axis
magY = abs(Y)/n;             % Magnitude of FFT

% One-sided spectrum (due to symmetry)
half_n = floor(n/2);
f_plot = f_axis(1:half_n);
magY_plot = 2*magY(1:half_n);   % Multiply by 2 to conserve energy

% Plot time-domain signal
figure;
subplot(2,1,1);
plot(t, y, 'b');
xlabel('Time (s)');
ylabel('Amplitude');
title('Time-Domain Signal');
grid on;

% Plot frequency-domain signal (spectrum)
subplot(2,1,2);
plot(f_plot, magY_plot, 'r');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum via FFT');
grid on;

% Mark peak frequency
[~, idx_peak] = max(magY_plot);
peak_freq = f_plot(idx_peak);
hold on;
plot(peak_freq, magY_plot(idx_peak), 'ko', 'MarkerFaceColor', 'g');
text(peak_freq, magY_plot(idx_peak), sprintf('  Peak: %.1f Hz', peak_freq));

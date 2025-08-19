kick_polariton_omega_0_3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What does 0 frequency mean?
%A 0 Hz (zero frequency) component in the Fourier Transform corresponds to a constant, 
%non-oscillating part of the signal — in other words, its average (mean) value over time.
%It's also called the DC component 
%(short for "direct current", from electrical engineering terminology).

%Why does a peak appear at 0 frequency?
%Let's say you have a signal f(t)
%The Fourier Transform of that signal, F(ω), represents 
% how much of each frequency exists in f(t).
%If the average value of f(t) is not zero, then there is a constant (flat) part in your signal.
%This constant part corresponds to 0 Hz, and so the transform shows a large value (a peak) at zero frequency.

%Case 1: Constant Signal
%Let's say : f(t) = 5
%This is just a flat line --- a constant signal
%The Fourier Transform of this is a delta function at 0 frequency with amplitude 5.
%Why? Because the signal doesn't oscillate; it's purely "0 Hz".

%Case 2: Sinusoidal Signal + Offset
%f(t) = sin(2*pi*10*t) + 3
%Here : 
%The sine wave is at 10Hz
%The 3 is a constant offset.

%So the Fourier transform will show:
%A peak at 10 Hz
%A big spike at 0 Hz due to the 3 offset

%How to remove it?
%If you're not interested in the DC component
%Subtract the mean of the signal before applying the Fouier Transform
%f_{zero-mean}(t) = f(t) - mean(f)
%Then do the FFT. The peak at 0 will disappear.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


inten_0 = abs(psiC_x0_V_0).^2;
inten_2 = abs(psiC_x0_V_2).^2;
% inten_0 = inten_0 - mean(inten_0);
% inten_2 = inten_2 - mean(inten_2);

[f_V0, mag_V0] = spect_from_wavepac(inten_0_new,tt);
[f_V2, mag_V2] = spect_from_wavepac(inten_2_new,tt);



figure()
subplot(4,1,1)
plot(tt, abs(psiC_x0_V_0).^2)
title('|\Psi_c(x=0,t)|^2, Vx=0',FontSize=10)
xlabel('time')

subplot(4,1,2)
plot(2*pi*f_V0, mag_V0)
xlabel('frequency')
xlim([0,15])

subplot(4,1,3)
plot(tt, abs(psiC_x0_V_2).^2)
title('|\Psi_c(x=0,t)|^2, Vx=2',FontSize=10)
xlabel('time')

subplot(4,1,4)
plot(2*pi*f_V2, mag_V2)
xlim([0,15])
xlabel('frequency')

function [freq, spec] = spect_from_wavepac(wave, t)


Nt = numel(t);
dt = t(2) - t(1);
fs = 1/dt;

% FFT (shift both spectrum and frequency grid)
X = fftshift(fft(wave));
spec = abs(X).^2;
spec = spec / max(spec);

% Frequency grid centered at zero
if mod(Nt,2)==0
    k = -Nt/2 : Nt/2-1;
else
    k = -(Nt-1)/2 : (Nt-1)/2;
end
freq = k * (fs/Nt);

%freq = 2pifreq;
end

% function [f, magnitude] = getFrequencySpectrum(t, signal)
% % getFrequencySpectrum - Computes FFT of a time-domain signal
% % 
% % Syntax:
% %   [f, magnitude] = getFrequencySpectrum(t, signal)
% %
% % Inputs:
% %   t      - Time vector (1D array)
% %   signal - Signal values corresponding to time vector (1D array)
% %
% % Outputs:
% %   f         - Frequency vector (Hz)
% %   magnitude - Magnitude spectrum (amplitude at each frequency)
% 
%     % Ensure column vectors
%     t = t(:);
%     signal = signal(:);
% 
%     % Time step and sampling frequency
%     dt = mean(diff(t));
%     fs = 1 / dt;
%     n = length(signal);
% 
%     % Compute FFT
%     Y = fft(signal);
%     magY = abs(Y) / n;  % Normalize magnitude
%     magY = magY(1:floor(n/2));  % One-sided spectrum
%     magY = 2 * magY;  % Multiply by 2 to preserve energy (except DC)
% 
%     % Frequency vector
%     f = (0:floor(n/2)-1) * (fs / n);
% 
%     % Output
%     magnitude = magY;
% end

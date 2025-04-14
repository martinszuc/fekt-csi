function [psd, f] = myperiodogram(x, Fs)
%% Odhad PSD pomocí peridogramu
% [psd, f] = myperiodogram(x, Fs)
% x   - vektor vzorků vstupního signálu
% Fs  - vzorkovací kmitočet
%
% psd - vektor odhadnuté výkonové spektrální hustoty
% f   - kmitočtová osa

% Number of samples
N = length(x);

% Calculate FFT
NFFT = 2^nextpow2(N);
X = fft(x, NFFT);

% Calculate periodogram
P = (1/N) * abs(X).^2;

% Frequency axis
f = (-NFFT/2:NFFT/2-1)*Fs/NFFT;

% Shift to center zero frequency
psd = fftshift(P);
end
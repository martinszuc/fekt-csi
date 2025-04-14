function [psd, f] = myblatuk(x, Fs)
%% Odhad PSD pomocí Blackmanovy-Tukeyovy metody
% [psd, f] = myblatuk(x, Fs)
% x   - vektor vzorků vstupního signálu
% Fs  - vzorkovací kmitočet
%
% psd - vektor odhadnuté výkonové spektrální hustoty
% f   - kmitočtová osa

% Calculate autocorrelation function
[r, lags] = xcorr(x, 'biased');  % 'biased' means r is divided by length(x)

% Number of samples
N = length(x);

% Calculate FFT of autocorrelation function
NFFT = 2^nextpow2(length(r));
R = fft(r, NFFT);

% Calculate double-sided spectrum
P = abs(R);

% Frequency axis
f = (-NFFT/2:NFFT/2-1)*Fs/NFFT;

% Shift to center zero frequency
psd = fftshift(P);
end
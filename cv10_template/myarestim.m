function [psd, f] = myarestim(x, Fs, P)
%% Odhad PSD pomocí AR modelu
% [psd, f] = myarestim(x, Fs, P)
% x   - vektor vzorků vstupního signálu
% Fs  - vzorkovací kmitočet
% P   - řád modelu
%
% psd - vektor odhadnuté výkonové spektrální hustoty
% f   - kmitočtová osa

% Default model order
if(3 > nargin || isempty(P))
    P = 10;
end

% Estimate AR model coefficients using Yule-Walker method
[a, variance] = aryule(x, P);

% Number of frequency points for spectrum calculation
NFFT = 512;

% Calculate frequency response
[h, f_normalized] = freqz(1, a, NFFT, 'whole');

% Convert normalized frequency to Hz
f = f_normalized * Fs / (2*pi);
f = f - Fs/2; % Shift to center at zero

% Calculate PSD
psd = sqrt(variance) * abs(h);
end
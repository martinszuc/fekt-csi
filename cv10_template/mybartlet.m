function [psd, f] = mybartlet(x, Fs, K)
%% Odhad PSD pomocí Bartlettovy metody
% [psd, f] = mybartlet(x, Fs, K)
% x   - vektor vzorků vstupního signálu
% Fs  - vzorkovací kmitočet
% K   - počet segmentů
%
% psd - vektor odhadnuté výkonové spektrální hustoty
% f   - kmitočtová osa

% Default number of segments
if(3 > nargin || isempty(K))
    K = 5;
end

% Total length of signal
N = length(x);

% Length of each segment
L = floor(N/K);

% Initialize matrix for storing periodograms
P = zeros(K, L);

% Calculate periodogram for each segment
for i = 0:K-1
    % Extract segment
    segment = x(i*L+1:(i+1)*L);
    
    % Calculate FFT
    X = fft(segment, L);
    
    % Calculate periodogram and store in row i+1
    P(i+1,:) = (1/L) * abs(X).^2;
end

% Average periodograms
psd_unsorted = mean(P, 1);

% Frequency axis
f = (-L/2:L/2-1)*Fs/L;

% Shift to center zero frequency
psd = fftshift(psd_unsorted);
end
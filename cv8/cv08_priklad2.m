%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cv08_zadani_2.m
% Full MATLAB code for "Příklad 2" from the docs:
%  1) Loads the hychirp signal
%  2) Plots and plays it
%  3) Computes & plots its FT (via mydft)
%  4) Computes & plots its STFT (spectrogram)
%  5) Computes & plots its CWT (cwt)
% Compare the time-frequency representations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% Load the hychirp signal
% The file "hychirp.mat" should be on your MATLAB path or in the current folder.
load('hychirp.mat');  % this should give a variable named hychirp

% If you do not know the sampling frequency of hychirp, you can guess or
% check documentation for an official value. Here we'll assume:
Fs = 2000;  % sample rate in Hz (adjust if your data says otherwise)

% Put hychirp in a column vector
x = hychirp(:);

%% Time-domain plot
N = length(x);
t = (0:N-1)/Fs;  % time axis in seconds

figure
plot(t, x)
title('Time-domain Plot of hychirp')
xlabel('Time (s)')
ylabel('Amplitude')
grid on

%% Listen to the signal (uncomment if you want to hear it)
%soundsc(x, Fs);

%% 1) Compute & plot the Fourier Transform (FT)
[X, faxis] = mydft(x, Fs);

figure
subplot(2,1,1)
plot(faxis, abs(X))
grid on
title('Magnitude of FT (hychirp)')
xlabel('Frequency (Hz)')
ylabel('|X(f)|')

subplot(2,1,2)
plot(faxis, angle(X)*180/pi)
grid on
title('Phase of FT (hychirp)')
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')

%% 2) Short-time Fourier transform (STFT) via spectrogram
% Try different window lengths, overlaps, etc. to see the effect on resolution.
winLength  = round(0.05 * Fs);   % 50 ms window
overlap    = round(0.4 * winLength);
nfft       = 2^nextpow2(winLength);  % typical NFFT choice

figure
spectrogram(x, hamming(winLength), overlap, nfft, Fs, 'yaxis');
title('STFT via spectrogram() (Window = 50 ms)')
colorbar

%% 3) Continuous Wavelet Transform (CWT)
% MATLAB’s built-in function "cwt" automatically displays a scalogram.
figure
cwt(x, Fs)
title('Continuous Wavelet Transform of hychirp')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local helper function: MYDFT
% Simple DFT routine similar to "myfft" or "mydft" from the examples.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, f] = mydft(x, Fs, f)
    if nargin < 1
        error('Must provide signal vector x at least.');
    end
    x = x(:);
    N = length(x);

    % Default Fs if not given
    if nargin < 2
        Fs = 1;
    end

    % If frequencies are not given, use a default from -Fs/2..Fs/2
    if nargin < 3
        f = ((0:N-1) - floor((N-1)/2))/N * Fs;
    end

    % If scalar "f" is given, interpret it as max frequency: -f..f
    if isscalar(f)
        f = -f : Fs/N : f;
    end

    % Allocate result
    X = zeros(length(f), 1);
    n = (0:N-1).';

    % Compute the DFT sum
    for k = 1:length(f)
        X(k) = (1/N) * (x.' * exp(-1i*2*pi*f(k)*n/Fs));
    end
end

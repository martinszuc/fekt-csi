%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cv08_priklad2.m
% Příklad 2 – Fourier, STFT, and CWT analysis of the 'hychirp' signal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% Load the signal
load('hychirp.mat'); % should load variable named 'hychirp'
Fs = 200;            % hychirp is sampled at 200 Hz (as per MATLAB examples)
t = (0:length(hychirp)-1)/Fs;

%% Plot and play the signal
figure;
plot(t, hychirp)
title('HyChirp signal in time domain')
xlabel('Time [s]')
ylabel('Amplitude')
grid on

% soundsc(hychirp, Fs); % Uncomment to listen

%% Compute DFT using custom function
disp('Computing DFT using mydft...');
figure;
mydft(hychirp, Fs);

%% STFT using spectrogram (window length 50 ms and 200 ms)
winLens = [0.05, 0.2];  % in seconds
for i = 1:length(winLens)
    winLenSamples = round(winLens(i)*Fs);
    window = hann(winLenSamples);
    overlap = round(0.5 * winLenSamples); % 50% overlap
    [S, F, T_stft] = spectrogram(hychirp, window, overlap, [], Fs);

    figure;
    imagesc(T_stft, F, 10*log10(abs(S).^2))
    axis xy
    colormap jet
    colorbar
    title(['STFT with window length = ' num2str(winLens(i)*1000) ' ms'])
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
end

%% Continuous Wavelet Transform (CWT)
figure;
cwt(hychirp, Fs);
title('Continuous Wavelet Transform of HyChirp')

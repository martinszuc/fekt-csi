%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cv08_zadani_1.m
% Full MATLAB code for Příklad 1:
% Demonstrates generating a multi-component signal, computing its DFT,
% and performing STFT with varying window lengths.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% Basic parameters
Fs = 4000;              % sampling frequency
Ts = 1/Fs;              % sampling period
T  = 2;                 % total signal duration in seconds
N  = T * Fs;            % total number of samples
t  = 0 : Ts : T - Ts;   % time vector

%% Generate signals
f0 = 1000; f1 = 1250; f2 = 1500; f3 = 1750; f4 = 500;
sin0 = sin(2*pi*f0*t);
sin1 = sin(2*pi*f1*t);
sin2 = sin(2*pi*f2*t);
sin3 = sin(2*pi*f3*t);
sin4 = sin(2*pi*f4*t);

%% Construct a signal with interruptions and frequency segments
% Short "blip" from sin4 repeated, plus segments of sin0..sin3
y2 = [zeros(N/100, 1)'  sin4(1 : N/100)]';    % short piece
y2 = repmat(y2, 50, 1);                      % repeat that piece
y  = y2 + [ ...
    sin0(1 : floor(N/4)-1),          ...
    sin1(floor(N/4) : floor(N/2)-1), ...
    sin2(floor(N/2) : floor(3*N/4)-1), ...
    sin3(floor(3*N/4) : end)        ...
    ]';

%% Optional playback
% soundsc(y, Fs);

%% Compute DFT for the entire 2 s signal
% (Calls our local mydft function below.)
mydft(y, Fs);

%% Plot STFT with different window lengths
figure

% Window length = 0.1 s
subplot(2,2,1)
winLength = 0.1;  % seconds
[fAxis, SpecCoeffs] = MySTFT(y, Fs, winLength);
imagesc(linspace(0, T, size(SpecCoeffs,2)), fAxis', log10(abs(SpecCoeffs).^2));
title(['Délka okna: ' num2str(winLength) ' s'])
set(gca,'YDir','normal')
xlabel('t (s)')
ylabel('f (Hz)')

% Window length = 0.05 s
subplot(2,2,2)
winLength = 0.05; % seconds
[fAxis, SpecCoeffs] = MySTFT(y, Fs, winLength);
imagesc(linspace(0, T, size(SpecCoeffs,2)), fAxis', log10(abs(SpecCoeffs).^2));
title(['Délka okna: ' num2str(winLength) ' s'])
set(gca,'YDir','normal')
xlabel('t (s)')
ylabel('f (Hz)')

% Window length = 0.02 s
subplot(2,2,3)
winLength = 0.02; % seconds
[fAxis, SpecCoeffs] = MySTFT(y, Fs, winLength);
imagesc(linspace(0, T, size(SpecCoeffs,2)), fAxis', log10(abs(SpecCoeffs).^2));
title(['Délka okna: ' num2str(winLength) ' s'])
set(gca,'YDir','normal')
xlabel('t (s)')
ylabel('f (Hz)')

% Window length = 0.01 s
subplot(2,2,4)
winLength = 0.01; % seconds
[fAxis, SpecCoeffs] = MySTFT(y, Fs, winLength);
imagesc(linspace(0, T, size(SpecCoeffs,2)), fAxis', log10(abs(SpecCoeffs).^2));
title(['Délka okna: ' num2str(winLength) ' s'])
set(gca,'YDir','normal')
xlabel('t (s)')
ylabel('f (Hz)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local function: MYSTFT
% Computes the Short-Time Fourier Transform for a given signal x,
% sampling rate Fs, and window length WinLength (in seconds).
% Returns the vector of positive frequencies 'f' and
% the STFT matrix 'SpecCoeffs'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, SpecCoeffs] = MySTFT(x, Fs, WinLength)
    % Convert window length from seconds to samples
    L = round(WinLength * Fs);
    
    % Hann window
    w = hann(L,'periodic');
    
    % Ensure x is a column
    x = x(:);
    
    % Zero-pad x if needed so that length(x) is multiple of L
    Nx = length(x);
    remainder = mod(Nx, L);
    if remainder ~= 0
        x = [x; zeros(L - remainder, 1)];
        Nx = length(x);
    end
    
    % Number of frames (no overlap used here)
    numFrames = Nx / L;
    
    % Allocate space for STFT matrix
    SpecCoeffs = zeros(floor(L/2)+1, numFrames);
    
    % Frequency vector for entire FFT
    freqAxis = (0 : L-1)*(Fs/L);
    
    % Loop over frames
    for idx = 1 : numFrames
        startIdx = (idx-1)*L + 1;
        endIdx   = startIdx + L - 1;
        
        segment = x(startIdx:endIdx);
        segmentWindowed = segment .* w;
        X = fft(segmentWindowed);
        
        % Keep only positive frequencies
        SpecCoeffs(:, idx) = X(1 : floor(L/2) + 1);
    end
    
    % Return frequency axis for positive frequencies as a column
    f = freqAxis(1 : floor(L/2) + 1).';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local function: MYDFT
% Example discrete Fourier transform. 
% Usage: mydft(signal, Fs)
% or:    [X,f] = mydft(signal, Fs, fVector)
% Where fVector can be either a vector of frequencies or a max frequency.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, f] = mydft(x, Fs, f)
    if nargin < 1
        error('Must provide signal vector x at least.');
    end
    
    % Ensure x is column vector
    x = x(:);
    N = length(x);
    
    % Default Fs if not given
    if nargin < 2
        Fs = 1;
    end
    
    % If f is not given, compute DFT for 0..(N-1)/N * Fs 
    if nargin < 3
        f = ((0:N-1) - floor((N-1)/2))/N * Fs;
    end
    
    % If f is a single scalar, interpret it as the max frequency
    if isscalar(f)
        f = -f : Fs/N : f;
    end
    
    % Allocate result
    X = zeros(length(f), 1);
    n = (0 : N-1).';
    
    % Compute the DFT
    for k = 1 : length(f)
        X(k) = (1/N) * (x.' * exp(-1i*2*pi*f(k)*n / Fs));
    end
    
    % If no output requested, plot
    if nargout < 1
        figure
        subplot(2,1,1)
        plot(f, abs(X))
        grid on
        title('Magnitude of Discrete Fourier Transform')
        xlabel('Frequency [Hz]')
        ylabel('|X(f)|')
        
        subplot(2,1,2)
        plot(f, angle(X)*180/pi)
        grid on
        title('Phase of Discrete Fourier Transform')
        xlabel('Frequency [Hz]')
        ylabel('Phase [degrees]')
    end
end

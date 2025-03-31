%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cv08_zadani_1.m
% Řešení Příkladu 1 – Generování signálu a výpočet STFT s různou délkou okna
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% Základní proměnné
Fs = 4000;              % vzorkovací frekvence
Ts = 1/Fs;              % vzorkovací perioda
T  = 2;                 % délka signálu v sekundách
N  = T * Fs;            % celkový počet vzorků
t  = 0 : Ts : T - Ts;   % časový vektor

%% Generování harmonických signálů
f0 = 1000; f1 = 1250; f2 = 1500; f3 = 1750; f4 = 500;
sin0 = sin(2*pi*f0*t);
sin1 = sin(2*pi*f1*t);
sin2 = sin(2*pi*f2*t);
sin3 = sin(2*pi*f3*t);
sin4 = sin(2*pi*f4*t);

%% Přerušovaný signál
y2 = [zeros(N/100, 1)'  sin4(1 : N/100)]';
y2 = repmat(y2, 50, 1);

y  = y2 + [ ...
    sin0(1 : floor(N/4)-1),          ...
    sin1(floor(N/4) : floor(N/2)-1), ...
    sin2(floor(N/2) : floor(3*N/4)-1), ...
    sin3(floor(3*N/4) : end)        ...
    ]';

%% Přehrání signálu
%soundsc(y, Fs);

%% Výpočet DFT
mydft(y, Fs);

%% Výpočet STFT a zobrazení pro různé délky oken
figure

% Délka okna = 0.1 s
subplot(2,2,1)
winLength = 0.1;
[fAxis, SpecCoeffs] = MySTFT(y, Fs, winLength);
imagesc(linspace(0, T, size(SpecCoeffs,2)), fAxis', log10(abs(SpecCoeffs).^2));
title(['Délka okna: ' num2str(winLength) ' s'])
set(gca,'YDir','normal')
xlabel('t (s)')
ylabel('f (Hz)')

% Délka okna = 0.05 s
subplot(2,2,2)
winLength = 0.05;
[fAxis, SpecCoeffs] = MySTFT(y, Fs, winLength);
imagesc(linspace(0, T, size(SpecCoeffs,2)), fAxis', log10(abs(SpecCoeffs).^2));
title(['Délka okna: ' num2str(winLength) ' s'])
set(gca,'YDir','normal')
xlabel('t (s)')
ylabel('f (Hz)')

% Délka okna = 0.02 s
subplot(2,2,3)
winLength = 0.02;
[fAxis, SpecCoeffs] = MySTFT(y, Fs, winLength);
imagesc(linspace(0, T, size(SpecCoeffs,2)), fAxis', log10(abs(SpecCoeffs).^2));
title(['Délka okna: ' num2str(winLength) ' s'])
set(gca,'YDir','normal')
xlabel('t (s)')
ylabel('f (Hz)')

% Délka okna = 0.01 s
subplot(2,2,4)
winLength = 0.01;
[fAxis, SpecCoeffs] = MySTFT(y, Fs, winLength);
imagesc(linspace(0, T, size(SpecCoeffs,2)), fAxis', log10(abs(SpecCoeffs).^2));
title(['Délka okna: ' num2str(winLength) ' s'])
set(gca,'YDir','normal')
xlabel('t (s)')
ylabel('f (Hz)')

%% Funkce: MySTFT – výpočet krátkodobé Fourierovy transformace
function [f, SpecCoeffs] = MySTFT(x, Fs, WinLength)
    % Převod délky okna na počet vzorků
    L = round(WinLength * Fs);
    
    % Okno – Hannovo
    w = hann(L,'periodic');
    
    % Vstupní signál jako sloupcový vektor
    x = x(:);
    
    % Doplnění nulami, pokud délka není násobkem okna
    Nx = length(x);
    remainder = mod(Nx, L);
    if remainder ~= 0
        x = [x; zeros(L - remainder, 1)];
        Nx = length(x);
    end
    
    % Počet segmentů
    numFrames = Nx / L;
    
    % Předalokace STFT matice
    SpecCoeffs = zeros(floor(L/2)+1, numFrames);
    freqAxis = (0 : L-1)*(Fs/L);

    % Smyčka přes segmenty
    for idx = 1 : numFrames
        startIdx = (idx-1)*L + 1;
        endIdx   = startIdx + L - 1;
        
        segment = x(startIdx:endIdx) .* w;
        X = fft(segment);
        SpecCoeffs(:, idx) = X(1 : floor(L/2) + 1);
    end
    
    % Vrácení kladných frekvencí
    f = freqAxis(1 : floor(L/2) + 1).';
end

%% Funkce: mydft – výpočet diskrétní Fourierovy transformace
function [X, f] = mydft(x, Fs, f)
    if nargin < 1, error('Musí být zadán signál x'); end
    x = x(:);
    N = length(x);
    if nargin < 2, Fs = 1; end
    if nargin < 3, f = ((0:N-1) - floor((N-1)/2))/N * Fs; end
    if isscalar(f), f = -f : Fs/N : f; end

    X = zeros(length(f), 1);
    n = (0:N-1).';

    for k = 1:length(f)
        X(k) = (1/N) * (x.' * exp(-1i*2*pi*f(k)*n / Fs));
    end

    if nargout < 1
        figure
        subplot(2,1,1)
        plot(f, abs(X)), grid on, title('Modul DFT'), xlabel('f (Hz)'), ylabel('|X(f)|')
        subplot(2,1,2)
        plot(f, angle(X)*180/pi), grid on, title('Fáze DFT'), xlabel('f (Hz)'), ylabel('Fáze (°)')
    end
end

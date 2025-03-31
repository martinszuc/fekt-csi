%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cv08_zadani_2.m
% Řešení Příkladu 2 – Analýza signálu hychirp pomocí FT, STFT a CWT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% Načtení signálu hychirp
load('hychirp.mat');  % musí být v aktuální složce
x = hychirp(:);

% Vzorkovací frekvence (upravte dle potřeby)
Fs = 2000;

%% Časová osa a zobrazení signálu
N = length(x);
t = (0:N-1)/Fs;

figure
plot(t, x)
title('Signál hychirp v čase')
xlabel('Čas (s)')
ylabel('Amplituda')
grid on

%% Přehrání signálu
%soundsc(x, Fs);

%% Výpočet Fourierovy transformace
[X, faxis] = mydft(x, Fs);

figure
subplot(2,1,1)
plot(faxis, abs(X)), grid on
title('Modul Fourierovy transformace')
xlabel('f (Hz)')
ylabel('|X(f)|')

subplot(2,1,2)
plot(faxis, angle(X)*180/pi), grid on
title('Fáze Fourierovy transformace')
xlabel('f (Hz)')
ylabel('Fáze (°)')

%% Krátkodobá Fourierova transformace (STFT)
winLength = round(0.05 * Fs);       % délka okna: 50 ms
overlap   = round(0.4 * winLength); % překryv
nfft      = 2^nextpow2(winLength);

figure
spectrogram(x, hamming(winLength), overlap, nfft, Fs, 'yaxis');
title('STFT pomocí spectrogram() (okno 50 ms)')
colorbar

%% Kontinuální vlnková transformace (CWT)
figure
cwt(x, Fs)
title('Kontinuální vlnková transformace (CWT)')

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
end

clear
clc
close all

%% Základní proměnné
Fs = 4000;
Ts = 1/Fs;
T = 2;
N = T*Fs;
t = 0:Ts:T-Ts;

%% Generování signálů
f0 = 1000;
f1 = 1250;
f2 = 1500;
f3 = 1750;
f4 = 500;

sin0 = sin(f0*2*pi*t);
sin1 = sin(f1*2*pi*t);
sin2 = sin(f2*2*pi*t);
sin3 = sin(f3*2*pi*t);
sin4 = sin(f4*2*pi*t);

%% Přerušovaný signál
y2 = [zeros(N/100, 1)' sin4(1:N/100)]';
y2 = repmat(y2, 50, 1);

y = y2 + [sin0(1:floor((N/4)-1)) sin1(floor(N/4):floor((N/2)-1)) sin2(N/2:floor((3*N/4)-1)) sin3(3*N/4:end)]';

%% Přehrání
%soundsc(y, Fs);

%% DFT
mydft(y, Fs);

%% Vykreslení STFT

% délka okna = 0.1 s
figure
subplot(2, 2, 1)
winLength = 0.1;
[f, SpecCoeffs] = MySTFT(y, Fs, winLength);
hp = imagesc(linspace(0, T, size(SpecCoeffs, 2)), f', log10(abs(SpecCoeffs).^2));
title(['Délka okna: ' num2str(winLength) ' s'])
set(gca, 'YDir', 'normal')
xlabel('t (s)')
ylabel('f (Hz)')

% délka okna = 0.05 s
subplot(2, 2, 2)
xlabel('t (s)')
ylabel('f (Hz)')
winLength = 0.05;
[f, SpecCoeffs] = MySTFT(y, Fs, winLength);
imagesc(linspace(0, T, size(SpecCoeffs, 2)), f', log10(abs(SpecCoeffs).^2))
title(['Délka okna: ' num2str(winLength) ' s'])
set(gca, 'YDir', 'normal')
xlabel('t (s)')
ylabel('f (Hz)')

% délka okna = 0.02 s
subplot(2, 2, 3)
xlabel('t (s)')
ylabel('f (Hz)')
winLength = 0.02;
[f, SpecCoeffs] = MySTFT(y, Fs, winLength);
imagesc(linspace(0, T, size(SpecCoeffs, 2)), f', log10(abs(SpecCoeffs).^2))
title(['Délka okna: ' num2str(winLength) ' s'])
set(gca, 'YDir', 'normal')
xlabel('t (s)')
ylabel('f (Hz)')

% délka okna = 0.01 s
subplot(2, 2, 4)
xlabel('t (s)')
ylabel('f (Hz)')
winLength = 0.01;
[f, SpecCoeffs] = MySTFT(y, Fs, winLength);
imagesc(linspace(0, T, size(SpecCoeffs, 2)), f', log10(abs(SpecCoeffs).^2))
title(['Délka okna: ' num2str(winLength) ' s'])
set(gca, 'YDir', 'normal')
xlabel('t (s)')
ylabel('f (Hz)')

function [f, positiveFreqs] = MySTFT(x, Fs, WinLength)
    %% Návrh okna
    
    
    %% Segmentace signálu
    

    %% Doplnění signálu nulami pro nesoudělné délky okna
    
    
    %% For smyčka pro výpočet FFT pro každý segment
    

    %% Nastavení frekvencí (chceme pouze kladné)
    

end
%% Příklad 3.2: Rekonstrukce signálu pomocí 3 nejvýznamnějších FFT komponent
% Tento skript načte soubor data.mat, který obsahuje proměnné:
%   s  - sloupcový vektor s průběhem signálu
%   tn - vektor časové osy
%   N  - celková délka signálu s
%
% Skript pomocí FFT zjistí 3 nejvýznamnější frekvenční komponenty, 
% zrekonstruuje signál jako součet 3 harmonických funkcí a vykreslí 
% původní signál a rekonstruovaný signál v jednom grafu.

clear; close all; clc;

%% Načtení dat
load('data.mat');  % Ujistěte se, že data.mat je ve vaší pracovní složce

if ~exist('s', 'var') || ~exist('tn', 'var') || ~exist('N', 'var')
    error('V souboru data.mat chybí jedna nebo více požadovaných proměnných (s, tn, N).');
end

%% Parametry a FFT
% Vzorkovací frekvence vypočtená z časového vektoru
Fs = 1/mean(diff(tn));

% Výpočet FFT signálu
S = fft(s);
Nfft = length(s);  % Mělo by odpovídat hodnotě N

% Frekvenční osa: kompletní spektrum
f = (0:Nfft-1) * (Fs/Nfft);

% Pro rekonstrukci a analýzu používáme pouze jednosměrné (pozitivní) frekvence
half_N = floor(Nfft/2) + 1;
S_half = S(1:half_N);
f_half = f(1:half_N);

%% Výpočet amplitudového spektra
% Normalizace FFT (pro jednosměrné spektrum)
amp = abs(S_half)/Nfft;
% Komponenty kromě DC mají dvojnásobnou amplitudu (protože jsou symetrické)
amp(2:end) = 2 * amp(2:end);

%% Nalezení 3 nejvýznamnějších komponent
[sortedAmp, sortIdx] = sort(amp, 'descend');
top3Idx = sortIdx(1:3);
% Pro přehlednost se komponenty seřadí podle frekvence vzestupně
top3Idx = sort(top3Idx);

fprintf('3 nejvýznamnější komponenty:\n');
for i = 1:length(top3Idx)
    idx = top3Idx(i);
    freq = f_half(idx);
    amplitude = amp(idx);
    phase = angle(S_half(idx));
    fprintf('Komponenta %d: Frekvence = %.2f Hz, Amplituda = %.2f, Fáze = %.2f rad\n',...
        i, freq, amplitude, phase);
end

%% Rekonstrukce signálu pomocí 3 harmonických funkcí
s_reconstructed = zeros(size(s));
for i = 1:length(top3Idx)
    idx = top3Idx(i);
    freq = f_half(idx);
    % Pro DC (idx==1) není nutné zdvojovat amplitudu, u ostatních již byla amplituda upravena
    amplitude = amp(idx);
    phase = angle(S_half(idx));
    % Přidání příslušné harmonické složky
    s_reconstructed = s_reconstructed + amplitude * cos(2*pi*freq*tn + phase);
end

%% Vykreslení původního a rekonstruovaného signálu
figure;
plot(tn, s, 'b', 'LineWidth', 1.5);
hold on;
plot(tn, s_reconstructed, 'r--', 'LineWidth', 1.5);
hold off;
xlabel('Čas (s)');
ylabel('Amplituda');
title('Původní signál a rekonstruovaný signál (3 nejvýznamnější FFT komponenty)');
legend('Původní signál', 'Rekonstruovaný signál');
grid on;

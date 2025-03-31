function [X, f] = mydft(x, Fs, f)
%% MPC-CSI Cislicove zpracovani signalu
%
% Priklad vypoctu diskretni Fourierovy transformace.
%
% [X, f] = mydft(x, Fs, f)
%
%  x  - vektor vzorku vstupniho signalu
%  Fs - vzorkovaci kmitocet
%  f  - vektor kmitoctu, pro ktere se ma pocitat DFT
%       nebo nejvyssi kmitocet, do ktereho chci pocitat
%
%  X  - vzorky diskretni fourierovy transformace pro zadane kmitocty
%  f  - vektor kmitoctu, uzitecny, jestlize zadane f byl skalar oznacujici
%       nejvyssi kmitocet

%% Kontrola parametru
if nargin < 1
    error('Musi byt zadan alespon vektor vzorku signalu')
end
% Delka signalu - predpoklada se sloupcovy vektor ci matice, kde jsou
% signaly ulozene jako sloupce
N = size(x, 1);
% Predpokladame vzorkovaci kmitocet 1 Hz, neni-li zadany
if nargin < 2
    Fs = 1;
end
% Predpokladame vypocet pro N kmitoctu od 0 do (N-1)/N*Fs
if nargin < 3
    f = ((0:N-1) - floor((N-1)/2))/N*Fs;
end
if isscalar(f)
    f = -f:Fs/N:f;
end

%% Vygenerovani jadra transformace
n = 0:N-1;
n = n(:);
X = zeros(length(f), 1);

%% Vypocet Fourierovy transformace
for k = 0:length(f)-1
    X(k+1) = 1/N*x.'*exp(-1i*2*pi*f(k+1)*n/Fs);
end

%% Pokud nejsou zadany vystupni parametry, proved zobrazeni
if nargout < 1
    figure
    subplot(2, 1, 1)
    plot(f, abs(X))
    grid('on')
    title('Modul diskretni Fourierovy transformace')
    xlabel('\rightarrow {\it f} [Hz]')
    ylabel('\rightarrow |{\it X}(e^{i\omega})| [-]')
    subplot(2, 1, 2)
    plot(f, 180/pi*angle(X))
    grid('on')
    title('Faze diskretni Fourierovy transformace')
    xlabel('\rightarrow {\it f} [Hz]')
    ylabel('\rightarrow arg({\it X}(e^{i\omega})) [-]')
end
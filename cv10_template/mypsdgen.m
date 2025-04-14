function [x, Fs, a, d] = mypsdgen()
%% Generování náhodného signálu
% x  - vzorky vygenerovaného signálu
% Fs - vzorkovací kmitočet
% a  - koeficienty inovačního filtru
% r  - rozptyl vstupního bílého šumu

% parametry signálu
Fs = 8e3;
N = 2*Fs;

% rozptyl vstupního bílého šumu
d = 0.001;
% generování bílého šumu
w = sqrt(d) * randn(N,1);

% parametry inovačního procesu
a = [1];
% příklad modelu
%a = [1, -1.64, 1.68, -1.07, 0.28];
% parametry pro hlásku A
%w = 1e-2*w; a = [ 1.0000   -1.5204    1.1468   -1.2432    1.4715   -0.6261   -0.0690   -0.1181    0.0828    0.3242   -0.2634];
% parametry pro hlásku O
%w = 1e-2*w; a = [ 1.0000   -1.6583    1.1632   -1.7167    2.4952   -1.5786    0.9489   -1.2130    0.6905    0.1074   -0.1251];

% filtrace
x = filter(1, a, w);
end
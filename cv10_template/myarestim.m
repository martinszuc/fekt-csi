function [ psd, f] = myarestim( x, Fs, P)
%% Odhad PSD pomocí AR modelu
% [psd, f] = myarestim( x, Fs, P)
% x   - vektor vzorků vstupního signálu
% Fs  - vzorkovací kmitočet
% P   - řád modelu
%
% psd - vektor odhadnuté výkonové spektrální hustoty
% f   - kmitočtová osa

% výchozí řád modelu
if( 3 > nargin || isempty( P))
    P = 10;
end

% výpočet kmitočtové osy pro oboustranné spektrum
f = linspace( -1, 1);

psd = inf( size( f));


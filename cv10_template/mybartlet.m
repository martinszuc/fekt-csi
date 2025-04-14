function [ psd, f] = mybartlet( x, Fs, K)
%% Odhad PSD pomocí Bartlettovy metody
% [psd, f] = mybartlet( x, Fs, K)
% x   - vektor vzorků vstupního signálu
% Fs  - vzorkovací kmitočet
% K   - počet segmentů
%
% psd - vektor odhadnuté výkonové spektrální hustoty
% f   - kmitočtová osa

% výchozí počet segmentů
if( 3 > nargin || isempty( K))
    K = 5;
end

% výpočet kmitočtové osy pro oboustranné spektrum
f = linspace( -1, 1);

psd = inf( size( f));

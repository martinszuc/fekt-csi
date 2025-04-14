function [ psd, f] = myperiodogram( x, Fs)
%% Odhad PSD pomocí peridogramu
% [psd, f] = myperiodogram( x, Fs)
% x   - vektor vzorků vstupního signálu
% Fs  - vzorkovací kmitočet
%
% psd - vektor odhadnuté výkonové spektrální hustoty
% f   - kmitočtová osa

% výpočet kmitočtové osy pro oboustranné spektrum
f = linspace( -1, 1);

psd = inf( size( f));

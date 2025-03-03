close all
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Zde doplnit vstupní parametry %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Zde zavolat funkci dp a hp %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vstupní parametry funkcí dp a hp:
% 1) fm
% 2) fvz

% uvnitø funkce: vıpoèet zesílení g
% vyuijte funkci parallel (paralelní spojení);
% pøímá cesta (nebo dopøedná vazba) bez zesílení/zeslabení 

% vıstupní parametry funkce:
% 1) vektor koeficientù èitatele (Hc) 
% 2) vektor koeficientù jmenovatele (Hj)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Zde vypoèítat frekvenèní charakteristiku %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%
%%% Vykreslení %%%%
%%%%%%%%%%%%%%%%%%%
% prostudujte, jak je zaobrazení zadáno a pøizpùsobte tomu volání funkcí
% v pøedchozím kroku
% u modulového spektra pøidejte v grafu "tick", kterı bude roven -3 dB,
% a je rozlišeno propustné a nepropustné pásmo 
% (není tím myšlena úprava v grafickém oknì, ale ji samotné volání grafu)

figure

subplot(2,1,1)
semilogx(f, 20*log10(abs(H1)))
hold on
semilogx(f, 20*log10(abs(H2)))


subplot(2,1,2)
semilogx(f, unwrap(angle(H1))/pi) % angle normálnì odeèítá sudé násobky pi, proto unwrap (rozbalení fáze)
hold on
semilogx(f, unwrap(angle(H2))/pi)


subplot(2,1,1)
set(gca, 'xlim', [20 fvz/2])
set(gca, 'ylim', [-15 0])
title('\bfModulová kmitoètová charakteristika')
xlabel('{{\itf}} (Hz) \rightarrow')
ylabel('Modul (dB) \rightarrow')
legend('{\it\bfH}_D_P ({\itf})', '{\it\bfH}_H_P ({\itf})')
grid on

subplot(2,1,2)
set(gca, 'xlim', [20 fvz/2])
title('\bfFázová kmitoètová charakteristika')
xlabel('{{\itf}} (Hz) \rightarrow')
ylabel('Fáze (\pirad) \rightarrow')
legend('{\it\bfH}_D_P ({\itf})', '{\it\bfH}_H_P ({\itf})')
grid on
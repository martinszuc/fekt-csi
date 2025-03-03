close all
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Zde doplnit vstupn� parametry %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Zde zavolat funkci dp a hp %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vstupn� parametry funkc� dp a hp:
% 1) fm
% 2) fvz

% uvnit� funkce: v�po�et zes�len� g
% vyu�ijte funkci parallel (paraleln� spojen�);
% p��m� cesta (nebo dop�edn� vazba) bez zes�len�/zeslaben� 

% v�stupn� parametry funkce:
% 1) vektor koeficient� �itatele (Hc) 
% 2) vektor koeficient� jmenovatele (Hj)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Zde vypo��tat frekven�n� charakteristiku %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%
%%% Vykreslen� %%%%
%%%%%%%%%%%%%%%%%%%
% prostudujte, jak je zaobrazen� zad�no a p�izp�sobte tomu vol�n� funkc�
% v p�edchoz�m kroku
% u modulov�ho spektra p�idejte v grafu "tick", kter� bude roven -3 dB,
% a� je rozli�eno propustn� a nepropustn� p�smo 
% (nen� t�m my�lena �prava v grafick�m okn�, ale ji� samotn� vol�n� grafu)

figure

subplot(2,1,1)
semilogx(f, 20*log10(abs(H1)))
hold on
semilogx(f, 20*log10(abs(H2)))


subplot(2,1,2)
semilogx(f, unwrap(angle(H1))/pi) % angle norm�ln� ode��t� sud� n�sobky pi, proto unwrap (rozbalen� f�ze)
hold on
semilogx(f, unwrap(angle(H2))/pi)


subplot(2,1,1)
set(gca, 'xlim', [20 fvz/2])
set(gca, 'ylim', [-15 0])
title('\bfModulov� kmito�tov� charakteristika')
xlabel('{{\itf}} (Hz) \rightarrow')
ylabel('Modul (dB) \rightarrow')
legend('{\it\bfH}_D_P ({\itf})', '{\it\bfH}_H_P ({\itf})')
grid on

subplot(2,1,2)
set(gca, 'xlim', [20 fvz/2])
title('\bfF�zov� kmito�tov� charakteristika')
xlabel('{{\itf}} (Hz) \rightarrow')
ylabel('F�ze (\pirad) \rightarrow')
legend('{\it\bfH}_D_P ({\itf})', '{\it\bfH}_H_P ({\itf})')
grid on
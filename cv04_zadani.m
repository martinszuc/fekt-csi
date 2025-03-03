close all
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Zde doplnit vstupn� parametry %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fm = 500;
fvz = 22050;
fft_N = 1024;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Zde zavolat funkci dp a hp %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Hc_dp, Hj_dp] = dp(fm, fvz);
[Hc_hp, Hj_hp] = hp(fm, fvz);

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
f = linspace(20, fvz/2, fft_N/2);
H1 = freqz(Hc_dp, Hj_dp, 2*pi*f/fvz);
H2 = freqz(Hc_hp, Hj_hp, 2*pi*f/fvz);

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
yline(-3, 'r--', '-3 dB')

subplot(2,1,2)
semilogx(f, unwrap(angle(H1))/pi) % angle norm�ln� ode��t� sud� n�sobky pi, proto unwrap (rozbalen� f�ze)
hold on
semilogx(f, unwrap(angle(H2))/pi)


subplot(2,1,1)
set(gca, 'xlim', [20 fvz/2])
set(gca, 'ylim', [-15 0])
title('\bfModulovo kmitoctova charakteristika')
xlabel('{{\itf}} (Hz) \rightarrow')
ylabel('Modul (dB) \rightarrow')
legend('{\it\bfH}_D_P ({\itf})', '{\it\bfH}_H_P ({\itf})')
grid on

subplot(2,1,2)
set(gca, 'xlim', [20 fvz/2])
title('\bfFazovo kmitoctova charakteristika')
xlabel('{{\itf}} (Hz) \rightarrow')
ylabel('Faze (\pirad) \rightarrow')
legend('{\it\bfH}_D_P ({\itf})', '{\it\bfH}_H_P ({\itf})')
grid on

function [Hc, Hj] = dp(fm, fvz)
    g = (tan(pi*fm/fvz) - 1) / (tan(pi*fm/fvz) + 1);

    Hc_phase = [g, 1];
    Hj_phase = [1, g];

    Hc_direct = 0.5;
    Hj_direct = 1;

    [Hc, Hj] = parallel(Hc_direct, Hj_direct, 0.5*Hc_phase, Hj_phase);
end

function [Hc, Hj] = hp(fm, fvz)
    g = (tan(pi*fm/fvz) - 1) / (tan(pi*fm/fvz) + 1);

    Hc_phase = [g, 1];
    Hj_phase = [1, g];

    Hc_direct = 0.5;
    Hj_direct = 1;

    Hc_phase_neg = -0.5*Hc_phase;

    [Hc, Hj] = parallel(Hc_direct, Hj_direct, Hc_phase_neg, Hj_phase);
end

function [num, den] = parallel(b1, a1, b2, a2)
    b1 = b1(:).';
    a1 = a1(:).';
    b2 = b2(:).';
    a2 = a2(:).';

    num = conv(b1, a2) + conv(b2, a1);
    den = conv(a1, a2);

    if den(1) ~= 0
        num = num / den(1);
        den = den / den(1);
    end
end
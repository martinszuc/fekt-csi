% Priklad 8.4
% Mame-li nahodny proces s nulovou strednou hodnotou, akej inej statistickej 
% velicine odpoveda hodnota autokorelacie na vzdialenost 0, R(t,t)?

% Definicia problemu
fprintf('Priklad 8.4: Comu odpoveda hodnota autokorelacie R(t,t) pre nahodny proces s nulovou strednou hodnotou?\n\n');

% Teoreticke vysvetlenie
fprintf('Autokorelacna funkcia nahodneho procesu je definovana ako:\n');
fprintf('R(t1, t2) = E[X(t1) * X(t2)]\n\n');

fprintf('Pre t1 = t2 = t dostaneme:\n');
fprintf('R(t, t) = E[X(t) * X(t)] = E[X(t)^2]\n\n');

fprintf('Pre proces s nulovou strednou hodnotou, E[X(t)] = 0, plati:\n');
fprintf('Var[X(t)] = E[(X(t) - E[X(t)])^2] = E[X(t)^2 - 2*X(t)*E[X(t)] + E[X(t)]^2]\n');
fprintf('Var[X(t)] = E[X(t)^2] - 2*E[X(t)]*E[X(t)] + E[X(t)]^2\n');
fprintf('Var[X(t)] = E[X(t)^2] - 0 + 0\n');
fprintf('Var[X(t)] = E[X(t)^2] = R(t,t)\n\n');

fprintf('Zaver: Pre nahodny proces s nulovou strednou hodnotou, hodnota autokorelacie R(t,t)\n');
fprintf('       odpoveda rozptylu procesu Var[X(t)] v case t.\n\n');

% Numericka demonstracia pomocou simulacie
fprintf('Numericka demonstracia:\n');

% Generovanie nahodneho procesu s nulovou strednou hodnotou
n = 10000;  % dlzka procesu
X = randn(n, 1);  % standardny normalny proces s E[X] = 0, Var[X] = 1
X = X - mean(X);  % zabezpecuje presne nulovu strednu hodnotu

% Vypocet statistickych velicin
mean_X = mean(X);
var_X = var(X);
autocorr_0 = xcorr(X, 0, 'biased');  % autokorelacia pri oneskoreni 0

fprintf('Stredna hodnota: %.6f (teoreticky 0)\n', mean_X);
fprintf('Rozptyl: %.6f (teoreticky 1)\n', var_X);
fprintf('Autokorelacia pri oneskoreni 0: %.6f\n\n', autocorr_0);

% Vizualizacia nahodneho procesu
figure;
subplot(2,1,1);
plot(1:100, X(1:100));
title('Ukazka nahodneho procesu s nulovou strednou hodnotou');
xlabel('t');
ylabel('X(t)');
grid on;

subplot(2,1,2);
[acf, lags] = xcorr(X(1:100), 20, 'biased');
stem(lags, acf);
title('Autokorelacna funkcia');
xlabel('Oneskorenie');
ylabel('R(tau)');
grid on;

% Histogram pre vizualizaciu distribucii
figure;
histogram(X, 50, 'Normalization', 'pdf');
hold on;
x_range = linspace(min(X), max(X), 1000);
plot(x_range, normpdf(x_range, 0, sqrt(var_X)), 'r', 'LineWidth', 2);
title('Histogram nahodneho procesu');
xlabel('Hodnota');
ylabel('Hustota pravdepodobnosti');
legend('Histogram', 'Teoreticka distribucina funkcia');
grid on;

% Zaver
fprintf('Zaver: Hodnota autokorelacie na vzdialenost 0, R(t,t), pre nahodny proces\n');
fprintf('       s nulovou strednou hodnotou je rovnaka ako rozptyl tohto procesu.\n');
fprintf('       To vidime aj z numerickych vysledkov, kde hodnota R(0) je blizka hodnote\n');
fprintf('       rozptylu procesu.\n');
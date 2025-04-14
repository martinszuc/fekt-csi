%% Priklad 3: Diskretna nahodna premenna s nulovou strednou hodnotou a nulovym rozptylom
% Ako vyzera pravdepodobnostna funkcia diskretnej nahodnej veliciny,
% ktora ma nulovu strednu hodnotu a nulovy rozptyl?

fprintf('\nPriklad 3: Diskretna nahodna premenna s nulovou strednou hodnotou a nulovym rozptylom\n');

% Teoreticke riesenie:
% Nech X je diskretna nahodna premenna s hodnotami x_i a pravdepodobnostami p(x_i)
% Stredna hodnota: E(X) = sum(x_i * p(x_i)) = 0
% Rozptyl: Var(X) = E[(X - E(X))^2] = E[X^2] = sum(x_i^2 * p(x_i)) = 0
%
% Aby bol rozptyl nulovy, musi platit sum(x_i^2 * p(x_i)) = 0
% Pretoze x_i^2 >= 0 a p(x_i) > 0, jedina moznost aby suma bola 0 je, 
% ze pravdepodobnost p(x_i) = 0 pre vsetky x_i rozne od 0, alebo x_i = 0
% pre vsetky i s nenulovou pravdepodobnostou.
%
% Jedina moznost je teda, ze X = 0 s pravdepodobnostou 1 (degenerovania distribucuia)

fprintf('Pravdepodobnostna funkcia diskretnej nahodnej veliciny s nulovou strednou hodnotou\n');
fprintf('a nulovym rozptylom vyzera nasledovne:\n');
fprintf('p(0) = 1\n');
fprintf('p(x) = 0 pre vsetky x != 0\n\n');

% Demonstracia pomocou simulacie
x_values = [0];
p_values = [1];

% Vypocet strednej hodnoty a rozptylu na overenie
mean_value = sum(x_values .* p_values);
variance = sum((x_values - mean_value).^2 .* p_values);

fprintf('Overenie pomocou vypoctu:\n');
fprintf('Stredna hodnota: %.2f\n', mean_value);
fprintf('Rozptyl: %.2f\n', variance);

% Graficke znazornenie
figure;
stem(x_values, p_values, 'LineWidth', 2);
grid on;
title('Pravdepodobnostna funkcia s nulovou strednou hodnotou a nulovym rozptylom');
xlabel('x');
ylabel('p(x)');
xlim([-1, 1]);
ylim([0, 1.1]);

%% Priklad 4: Analyza platov
% Vypocet priemerneho a medianoveho platu pred a po zmene

fprintf('\nPriklad 4: Analyza platov\n');

% Platy pred zmenou
platy_pred = [20, 10, 20, 30, 5, 15, 20, 40, 10, 25];

% Platy po zmene (posledna hodnota zmenena z 25 na 50)
platy_po = [20, 10, 20, 30, 5, 15, 20, 40, 10, 50];

% Vypocet priemerneho platu
priemerny_plat_pred = mean(platy_pred);
priemerny_plat_po = mean(platy_po);

fprintf('Priemerny plat pred zmenou: %.1f\n', priemerny_plat_pred);
fprintf('Priemerny plat po zmene: %.1f\n', priemerny_plat_po);
fprintf('Rozdiel v priemernych platoch: %.1f\n', priemerny_plat_po - priemerny_plat_pred);

% Vypocet medianoveho platu
medianovy_plat_pred = median(platy_pred);
medianovy_plat_po = median(platy_po);

fprintf('Medianovy plat pred zmenou: %.1f\n', medianovy_plat_pred);
fprintf('Medianovy plat po zmene: %.1f\n', medianovy_plat_po);
fprintf('Rozdiel v medianovych platoch: %.1f\n', medianovy_plat_po - medianovy_plat_pred);

% Graficke porovnanie platov pred a po zmene
figure;

% Zobrazenie platov pred zmenou
subplot(2, 1, 1);
bar(platy_pred, 'FaceColor', [0.3 0.6 0.9]);
hold on;
line([0, 11], [priemerny_plat_pred, priemerny_plat_pred], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
line([0, 11], [medianovy_plat_pred, medianovy_plat_pred], 'Color', 'g', 'LineWidth', 2, 'LineStyle', '-.');
title('Platy pred zmenou');
xlabel('Zamestnanec');
ylabel('Plat');
legend('Plat', 'Priemer', 'Median', 'Location', 'best');
ylim([0, 55]);
grid on;

% Zobrazenie platov po zmene
subplot(2, 1, 2);
bar(platy_po, 'FaceColor', [0.8 0.4 0.2]);
hold on;
line([0, 11], [priemerny_plat_po, priemerny_plat_po], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
line([0, 11], [medianovy_plat_po, medianovy_plat_po], 'Color', 'g', 'LineWidth', 2, 'LineStyle', '-.');
title('Platy po zmene');
xlabel('Zamestnanec');
ylabel('Plat');
legend('Plat', 'Priemer', 'Median', 'Location', 'best');
ylim([0, 55]);
grid on;

% Porovnanie distribucii platov
figure;
subplot(1, 2, 1);
boxplot([platy_pred', platy_po'], 'Labels', {'Pred zmenou', 'Po zmene'});
title('Porovnanie distribucii platov');
ylabel('Plat');
grid on;

subplot(1, 2, 2);
histogram(platy_pred, 0:5:55, 'FaceColor', [0.3 0.6 0.9], 'FaceAlpha', 0.5);
hold on;
histogram(platy_po, 0:5:55, 'FaceColor', [0.8 0.4 0.2], 'FaceAlpha', 0.5);
title('Histogram platov');
xlabel('Plat');
ylabel('Pocet zamestnancov');
legend('Pred zmenou', 'Po zmene', 'Location', 'best');
grid on;

% Analyza vplyvu zmeny na statistiky
fprintf('\nAnalyza vplyvu zmeny:\n');
fprintf('Zmena jednej hodnoty z 25 na 50 zvysila priemer o %.1f jednotiek.\n', priemerny_plat_po - priemerny_plat_pred);
fprintf('To je %.1f%% narust priemerneho platu.\n', 100*(priemerny_plat_po/priemerny_plat_pred - 1));

if medianovy_plat_po ~= medianovy_plat_pred
    fprintf('Median sa zmenil o %.1f jednotiek.\n', medianovy_plat_po - medianovy_plat_pred);
    fprintf('To je %.1f%% narust medianoveho platu.\n', 100*(medianovy_plat_po/medianovy_plat_pred - 1));
else
    fprintf('Median sa nezmenil, co ukazuje, ze je robustnejsi voci odlahlejsim hodnotam.\n');
end

% Zaver
fprintf('\nZaver:\n');
fprintf('Tento priklad ilustruje, preco median moze byt vhodnejsou statistikou ako priemer\n');
fprintf('pri analyze platov. Median nie je tak citlivy na extremne hodnoty, zatial co\n');
fprintf('priemer moze byt vyrazne ovplyvneny niekolkymi vysokymi alebo nizkymi hodnotami.\n');
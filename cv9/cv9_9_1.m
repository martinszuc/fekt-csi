% Priklad 9.1 - minimalna verzia pre MATLAB R2024b
% Tato verzia je velmi zjednodusena a pouziva len najzakladnejsie funkcie, aby fungovala
% na novych verziach MATLABu

% Nacitanie dat
fprintf('Nacitavam data z data.mat...\n');
load('data.mat');
fprintf('Data uspesne nacitane: %d hodnot\n', length(data));

% Rozdelenie dat
train_length = 72;
test_length = length(data) - train_length;
train_data = data(1:train_length);
test_data = data(train_length+1:end);

% Zobrazenie dat
figure;
plot(1:length(data), data);
hold on;
plot([train_length, train_length], [min(data), max(data)], 'r--');
title('Data s rozdelenim na trenovacie a testovacie casti');
legend('Data', 'Hranica');
grid on;

% Odhad AR modelu
p = 24;  % rad AR modelu
ar_coeffs = aryule(train_data, p);
fprintf('AR model radu %d uspesne odhadnuty\n', p);

% Predikcia - bez pouzitia filter funkcie a pociatocnych podmienok
% Pouzijeme jednoduchy algoritmus predikcie
pred = zeros(test_length, 1);
history = train_data(end-p+1:end);  % poslednych p vzoriek trenovacich dat

for i = 1:test_length
    % Predikcia buducej hodnoty pomocou AR modelu
    pred_value = -ar_coeffs(2:end) * flipud(history);
    
    % Ulozenie predikcie
    pred(i) = pred_value;
    
    % Aktualizacia historie
    history = [history(2:end); pred_value];
end

% Alternativny sposob, pouzitie celej postupnosti s postupnym doplnenim
pred_alt = zeros(test_length, 1);
full_sequence = [train_data; zeros(test_length, 1)];

for i = 1:test_length
    idx = train_length + i;
    % Vypocet novej hodnoty
    full_sequence(idx) = -ar_coeffs(2:end) * flipud(full_sequence(idx-p:idx-1));
    % Ulozenie predikcie
    pred_alt(i) = full_sequence(idx);
end

% Kontrola, ze oba sposoby daju rovnaky vysledok
fprintf('Rozdiel medzi metodami predikcie: %e\n', norm(pred - pred_alt));

% Vypocet chyby predikcie
mse = mean((pred - test_data).^2);
rmse = sqrt(mse);
fprintf('Stredna kvadraticka chyba (MSE): %.6f\n', mse);
fprintf('Odmocnina strednej kvadratickej chyby (RMSE): %.6f\n', rmse);

% Zobrazenie predikcie
figure;
plot(1:length(data), data);
hold on;
plot(train_length+1:length(data), pred, 'r--', 'LineWidth', 2);
plot([train_length, train_length], [min(data), max(data)], 'k--');
title(sprintf('AR(%d) predikcia', p));
legend('Skutocne data', 'Predikcia', 'Hranica');
grid on;

% Detailnejsie porovnanie
figure;
plot(train_length+1:length(data), test_data, 'b-', 'LineWidth', 1.5);
hold on;
plot(train_length+1:length(data), pred, 'r--', 'LineWidth', 2);
title('Porovnanie predikcie so skutocnymi hodnotami');
legend('Skutocne data', 'Predikcia');
xlabel('Index');
ylabel('Hodnota');
grid on;

% Zaver
fprintf('\nZaver:\n');
fprintf('- AR model radu %d bol pouzity na predikciu %d hodnot.\n', p, test_length);
fprintf('- Celkova chyba predikcie (RMSE) je %.6f.\n', rmse);
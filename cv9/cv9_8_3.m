% Priklad 8.3
% Model nahodneho procesu je definovany pomocou prenosovej funkcie inovacnej reprezentacie
% H(z) = (1 + 2z^(-2)) / (1 + z^(-1) - 3z^(-3))

% Definicia prenosovej funkcie H(z)
num = [1, 0, 2];    % citatej (1 + 2z^(-2))
den = [1, 1, 0, -3]; % menovatej (1 + z^(-1) - 3z^(-3))

% Klasifikacia modelu
fprintf('Klasifikacia modelu na zaklade prenosovej funkcie:\n');
fprintf('H(z) = (1 + 2z^(-2)) / (1 + z^(-1) - 3z^(-3))\n\n');

% Overime, ci model obsahuje AR a/alebo MA zlozku
has_AR = any(den(2:end) ~= 0);
has_MA = any(num(2:end) ~= 0);

if has_AR && has_MA
    fprintf('Model je typu ARMA, pretoze obsahuje aj AR, aj MA zlozku.\n\n');
elseif has_AR
    fprintf('Model je typu AR, pretoze obsahuje len AR zlozku.\n\n');
elseif has_MA
    fprintf('Model je typu MA, pretoze obsahuje len MA zlozku.\n\n');
else
    fprintf('Model je biely sum (nema AR ani MA zlozku).\n\n');
end

% Urcime rady zloziek
p = length(den) - 1;  % rad AR zlozky
q = length(num) - 1;  % rad MA zlozky

fprintf('Rad AR zlozky (p): %d\n', p);
fprintf('Koeficienty AR zlozky:\n');
for i = 2:length(den)
    fprintf('a%d = %g\n', i-1, den(i));
end
fprintf('\n');

fprintf('Rad MA zlozky (q): %d\n', q);
fprintf('Koeficienty MA zlozky:\n');
for i = 2:length(num)
    fprintf('b%d = %g\n', i-1, num(i));
end
fprintf('\n');

% Zjednoduseny zapis prenosovej funkcie
fprintf('Zjednoduseny zapis prenosovej funkcie:\n');
fprintf('H(z) = (1 + 2z^(-2)) / (1 + z^(-1) - 3z^(-3))\n\n');

% Diferencna rovnica
fprintf('Diferencna rovnica modelu:\n');
fprintf('x[n] + %g*x[n-1] + %g*x[n-2] + %g*x[n-3] = w[n] + %g*w[n-1] + %g*w[n-2]\n\n', ...
    den(2), den(3), den(4), num(2), num(3));

% Zobrazime polozero mapu
figure;
zplane(num, den);
title('Polozero mapa modelu H(z)');
grid on;

% Zobrazime frekvencnu charakteristiku
figure;
[h, w] = freqz(num, den, 1024);
subplot(2,1,1);
plot(w/pi, abs(h));
title('Amplitudova charakteristika H(z)');
xlabel('Normalizovana frekvencia (x\pi rad/vzorka)');
ylabel('Amplituda');
grid on;

subplot(2,1,2);
plot(w/pi, angle(h));
title('Fazova charakteristika H(z)');
xlabel('Normalizovana frekvencia (x\pi rad/vzorka)');
ylabel('Faza (rad)');
grid on;

% Generujeme kratku vzorku signalu z modelu pre vizualizaciu
n = 200;
noise = randn(n, 1);
signal = filter(num, den, noise);

figure;
subplot(2,1,1);
plot(1:n, noise);
title('Vstupny biely sum w[n]');
xlabel('n');
ylabel('Amplituda');
grid on;

subplot(2,1,2);
plot(1:n, signal);
title('Vystupny signal modelu x[n]');
xlabel('n');
ylabel('Amplituda');
grid on;

% Zaver
fprintf('Zaver: Model je typu ARMA(%d,%d).\n', p, q);
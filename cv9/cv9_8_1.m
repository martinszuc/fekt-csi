% Priklad 8.1
% AR model radu p = 5 s koeficientmi a = [1, 0, 0, 0, 2]

% Koeficienty AR modelu
a = [1, 0, 0, 0, 2];
p = length(a) - 1;  % rad modelu je 5

% Vypiseme diferencnu rovnicu
fprintf('Diferencna rovnica AR modelu radu p = %d:\n', p);
fprintf('x[n] + 0*x[n-1] + 0*x[n-2] + 0*x[n-3] + 2*x[n-5] = w[n]\n');
fprintf('alebo ekvivalentne:\n');
fprintf('x[n] = -2*x[n-5] + w[n]\n\n');

% Vypocitame a vykreslime inovacny filter H(z)
fprintf('Inovacny filter H(z):\n');
fprintf('H(z) = 1 / (1 + 0*z^(-1) + 0*z^(-2) + 0*z^(-3) + 2*z^(-5))\n');
fprintf('H(z) = 1 / (1 + 2*z^(-5))\n\n');

% Vypocitame a vykreslime belici filter G(z)
fprintf('Belici filter G(z):\n');
fprintf('G(z) = H^(-1)(z) = 1 + 2*z^(-5)\n\n');

% Zobrazime polozero mapu inovacneho filtra
figure;
zplane([1], a);
title('Polozero mapa inovacneho filtra H(z)');
grid on;

% Zobrazime frekvencnu charakteristiku inovacneho filtra
figure;
[h, w] = freqz(1, a, 1024);
subplot(2,1,1);
plot(w/pi, abs(h));
title('Amplitudova charakteristika inovacneho filtra H(z)');
xlabel('Normalizovana frekvencia (x\pi rad/vzorka)');
ylabel('Amplituda');
grid on;

subplot(2,1,2);
plot(w/pi, angle(h));
title('Fazova charakteristika inovacneho filtra H(z)');
xlabel('Normalizovana frekvencia (x\pi rad/vzorka)');
ylabel('Faza (rad)');
grid on;

% Generujeme kratku vzorku signalu z AR modelu pre vizualizaciu
n = 100;
noise = randn(n, 1);
signal = filter(1, a, noise);

figure;
subplot(2,1,1);
plot(1:n, noise);
title('Vstupny biely sum w[n]');
xlabel('n');
ylabel('Amplituda');
grid on;

subplot(2,1,2);
plot(1:n, signal);
title('Vystupny signal AR modelu x[n]');
xlabel('n');
ylabel('Amplituda');
grid on;
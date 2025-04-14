% Priklad 8.2
% MA model radu q = 4 s koeficientmi b = [1, -1, 0, 2]

% Koeficienty MA modelu
b = [1, -1, 0, 2];
q = length(b) - 1;  % rad modelu je 3

% Vypiseme diferencnu rovnicu
fprintf('Diferencna rovnica MA modelu radu q = %d:\n', q);
fprintf('x[n] = w[n] - w[n-1] + 0*w[n-2] + 2*w[n-3]\n');
fprintf('alebo ekvivalentne:\n');
fprintf('x[n] = w[n] - w[n-1] + 2*w[n-3]\n\n');

% Vypocitame a vypiseme inovacny filter H(z)
fprintf('Inovacny filter H(z):\n');
fprintf('H(z) = 1 - z^(-1) + 0*z^(-2) + 2*z^(-3)\n');
fprintf('H(z) = 1 - z^(-1) + 2*z^(-3)\n\n');

% Vypocitame a vypiseme belici filter G(z)
fprintf('Belici filter G(z):\n');
fprintf('G(z) = H^(-1)(z) = 1 / (1 - z^(-1) + 2*z^(-3))\n\n');

% Zobrazime polozero mapu inovacneho filtra
figure;
zplane(b, [1]);
title('Polozero mapa inovacneho filtra H(z)');
grid on;

% Zobrazime frekvencnu charakteristiku inovacneho filtra
figure;
[h, w] = freqz(b, 1, 1024);
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

% Generujeme kratku vzorku signalu z MA modelu pre vizualizaciu
n = 100;
noise = randn(n, 1);
signal = filter(b, 1, noise);

figure;
subplot(2,1,1);
plot(1:n, noise);
title('Vstupny biely sum w[n]');
xlabel('n');
ylabel('Amplituda');
grid on;

subplot(2,1,2);
plot(1:n, signal);
title('Vystupny signal MA modelu x[n]');
xlabel('n');
ylabel('Amplituda');
grid on;
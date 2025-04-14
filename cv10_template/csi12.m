function csi12()
% Výkonová spektrální hustota

close all
echo('off', 'all');

volba = 1;

while( volba)
    clear('all');

    volba = vyber();
    switch( volba)
        case 0
            volba = 0;
        case 1
            disp('Bílý šum');
            echo(mfilename, 'on');
            % počet průběhů
            M = 4;
            % délka signálu
            N = 1000;
            % vzorkovací kmitočet
            Fs = 8e3;
            % jednotlivé osy
            t = (0:N-1).'/Fs;
            n = (0:N-1).' - floor( N/2);
            fp = n/N*Fs;
            % modul
            A = zeros( size(fp));
            A(fp~=0) = 1;
            % fáze - musí být lichá
            phi = 2*pi*randn( ceil(N/2)*2-1, M);
            phi = phi - flipud( phi);
            phi = [zeros( N-length(phi), M); phi];
            % ideální spektrum
            Xi = A.*exp( j*phi);
            % výpočet časového průběhu
            x = real( ifft( N*fftshift( Xi)));
            % výpočet skutečného spektra
            X = fftshift( fft( x)/N);
            echo('off', 'all');

            figure( volba);
            set( gcf, 'Name', 'Bílý šum');
            plotsum( fp, [X, Xi], t, x);
        case 2
            disp('Růžový šum');
            echo(mfilename, 'on');
            % počet průběhů
            M = 4;
            % délka signálu
            N = 1000;
            % vzorkovací kmitočet
            Fs = 48e3;
            % jednotlivé osy
            t = (0:N-1).'/Fs;
            n = (0:N-1).' - floor( N/2);
            fp = n/N*Fs;
            % modul
            A = zeros( size(fp));
            A(fp~=0) = 1./sqrt(abs( fp(fp~=0)))*sqrt(Fs/N);
            % fáze - musí být lichá
            phi = 2*pi*randn( ceil(N/2)*2-1, M);
            phi = phi - flipud( phi);
            phi = [zeros( N-length(phi), M); phi];
            % ideální spektrum
            Xi = A.*exp( j*phi);
            % výpočet časového průběhu
            x = real( ifft( N*fftshift( Xi)));
            % výpočet skutečného spektra
            X = fftshift( fft( x)/N);
            echo('off', 'all');

            figure( volba);
            set( gcf, 'Name', 'Růžový šum');
            plotsum( fp, [X, Xi], t, x);
        case 3
            disp('Hnědý/červený šum');
            echo(mfilename, 'on');
            % počet průběhů
            M = 4;
            % délka signálu
            N = 1000;
            % vzorkovací kmitočet
            Fs = 8e3;
            % jednotlivé osy
            t = (0:N-1).'/Fs;
            n = (0:N-1).' - floor( N/2);
            fp = n/N*Fs;
            % modul
            A = zeros( size(fp));
            A(fp~=0) = 1./(abs( fp(fp~=0)))*(Fs/N);
            % fáze - musí být lichá
            phi = 2*pi*randn( ceil(N/2)*2-1, M);
            phi = phi - flipud( phi);
            phi = [zeros( N-length(phi), M); phi];
            % ideální spektrum
            Xi = A.*exp( j*phi);
            % výpočet časového průběhu
            x = real( ifft( N*fftshift( Xi)));
            % výpočet skutečného spektra
            X = fftshift( fft( x)/N);
            echo('off', 'all');

            figure( volba);
            set(gcf, 'Name', 'Hnědý šum');
            plotsum( fp, [X, Xi], t, x);
        case 4
            disp('Neparametrické metody');
            echo(mfilename, 'on');
            % generace signálu
            [ x, Fs, a, d] = mypsdgen();

            %% neprima metoda vypoctu
            [ psd1, f1] = myblatuk( x, Fs);
            %% prima metoda vypoctu
            [psd2, f2] = myperiodogram( x, Fs);
            %% Bartlettova metoda
            [psd3, f3] = mybartlet( x, Fs);
            %%
            echo('off', 'all');

            N = length( x);
            t = (0:N-1)/Fs;
            figure( volba);
            set( gcf, 'Name', 'Neparametrické metody');
            subplot( 2, 1, 1);
            plot( t, x);
            grid('on');
            title( 'Časový průběh');
            xlabel('\rightarrow t [s]');
            ylabel('\rightarrow x[nT] [-]');
            subplot( 2, 1, 2);
            % Kmitočtová charakteristika inovačního filtru
            H = freqz( 1, a, f1, Fs);
            plot( f1, 10*log10( d * abs( H).^2), 'LineWidth', 2);
            hold( 'on');
            plot( f1, 10*log10( psd1));
            plot( f2, 10*log10( psd2));
            plot( f3, 10*log10( psd3));
            hold( 'off');
            grid('on');
            title('Výkonová spektrální hustota');
            xlabel( '\rightarrow {\it f} [Hz]');
            ylabel( '\uparrow {\it P}_{\it  xx} [dB]');
            legend( ['Inovační proces'], ...
                ['Nepřímá metoda \Delta f = ' num2str( Fs/length( psd1))], ...
                ['Přímá metoda \Delta f = ' num2str( Fs/length( psd2))], ...
                ['Bartletova metoda \Delta f = ' num2str( Fs/length( psd3))]);

        case 5
            disp('Parametrické metody');
            echo(mfilename, 'on');
            % odhadovaný řád AR procesu
            P = 4;
            % generování signálu
            [x, Fs, a, d] = mypsdgen();
            %% odhad parametru
            [psd, f1] = myarestim( x, Fs);
            %%
            echo('off', 'all');

            %
            N = length( x);
            t = (0:N-1)/Fs;
            figure( volba);
            set( gcf, 'Name', 'Parametrické metody');
            subplot( 2, 1, 1);
            plot( t, x);
            grid('on');
            title( 'Časový průběh');
            xlabel('\rightarrow t [s]');
            ylabel('\rightarrow x[nT] [-]');
            subplot( 2, 1, 2);
            % Odhad periodogramem
            ff = ((0:N-1) - floor( N/2))*Fs/N;
            plot( ff, 10*log10( 1/N * fftshift( abs( fft( x)).^2)));
            hold( 'on');
            % Kmitočtová charakteristika inovačního filtru
            H = freqz( 1, a, f1, Fs);
            plot( f1, 10*log10( d * abs( H).^2), 'LineWidth', 2);
            % odhad výkonové spektrální hustoty
            plot( f1, 20*log10( psd));
            hold( 'off');
            grid('on');
            title('Výkonová spektrální hustota');
            xlabel( '\rightarrow {\it f} [Hz]');
            ylabel( '\uparrow {\it P}_{\it  xx} [dB]');
            legend( ['Periodogram'],...
                ['Inovační proces'], ...
                ['AR model']);

        otherwise
            disp('Neznámý příklad');
    end
end
end

function [c] = vyber()

drawnow();

disp( ' ');
disp( 'Výkonová spektrální hustota - příklady'); 
disp( '1 - Bílý šum');
disp( '2 - Růžový šum');
disp( '3 - Hnědý/červený šum');
disp( '4 - Neparametrické metody');
disp( '5 - Parametrické metody');
disp( '-------------------------------------');
disp( '0 - Konec');
disp( '-------------------------------------');
c = input( 'Vyber příklad [0]: ');
if( isempty( c))
    c = 0;
end

clc
disp(' ');

end

function plotsum(f, X, t, x)

subplot( 2, 2, 1);
semilogx( f( f > 0), 20*log10( abs( X(f>0,:))));
grid('on');
title( 'Modulová kmitočtová charakteristika hustota');
xlabel( '\rightarrow f [Hz]');
ylabel( '\rightarrow P_x(f) [dB]');

subplot( 2, 2, 3);
semilogx( f( f > 0), 180/pi*angle( X(f>0,:)));
grid('on');
title( 'Fázová kmitočtová charakteristika');
xlabel( '\rightarrow f [Hz]');
ylabel( '\rightarrow arg\{X(f)\} [°]');

subplot( 2, 2, 2);
plot( t, x);
grid('on');
title( 'Časový průběh');
xlabel( '\rightarrow n [-]');
ylabel( '\rightarrow x[nT] [s]');

subplot( 2, 2, 4);
r = zeros( 2*size( x, 1) - 1, size( x, 2));
[r(:,1), l] = xcorr( x(:,1),'biased');
for( k = 2:size( x, 2))
    r(:,k) = xcorr( x(:,k), 'biased');
end

plot( l, r);
grid('on');
title( 'Autokorelační funkce');
xlabel( '\rightarrow m [-]');
ylabel( '\rightarrow r_{xx}[m] [-]');
end
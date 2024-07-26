clc
clear 
close all

index_fig = 0;

%% Experimento 1

% Parametros do problema
M = 100;
f= 0.15;
k= 7;

t = linspace(0,25,1000); % Vetor de tempo
u = 1+cos(pi*(t/3)); % Sinal de entrada

lambda= 2; % Parametro do filtro

G= tf(1, [M f k]); % Resposta original da planta
x= lsim(G, u, t);

%% Experimento 1 - Item (a) 

s= tf("s");
filtro = tf(1 , [1 2*lambda lambda^2]); % Declaracao do filtro
z= lsim(filtro, u, t); % Sinal z do modelo parametrico

index_fig= index_fig + 1;
figure (index_fig) 
hold on 
grid on 
plot (t, u)
plot (t, z)
title ("Experimento 1(a) - Comparação entre u e z")
legend ('u', 'z')
xlabel ('Tempo [s]')
hold off

% Aplicacoes de filtro para o vetor phi
xdd_f= lsim(s^2*filtro, x, t);     
xd_f= lsim(s*filtro, x, t);
x_f= lsim(filtro, x, t); 

phi= [xdd_f xd_f x_f]'; %Vetor phi com filtro 

index_fig= index_fig + 1;
figure (index_fig)
hold on
grid on 
plot (t, phi)
title ("Experimento 1(a) - Sinal \phi")
hLeg = legend ('$\ddot{x}$', "$\dot{x}$", '$x$');
set(hLeg,'Interpreter','latex');
xlabel ('Tempo [s]')

index_fig= index_fig + 1;
figure (index_fig)
bode(filtro)
title("Experimento 1(a) - Resposta em frequÃªncia do filtro")
grid on 

%% Experimento 1 - Item (b)

% Novo modelo parametrico
z_b = z - M*xdd_f; % Novo sinal de z com M conhecido
phi_b = [xd_f x_f]'; % Vetor phi com filtro

% Derivadas do sinal de saída sem filtro 
xd= diff(x)./(diff(t)'); 
xdd= diff(xd)./(diff(t(1,1:end-1))');

z_sf= (u(1,1:end-2))' - M.*xdd;

index_fig= index_fig + 1;
figure (index_fig)
hold on 
grid on
plot(t, z_b)
plot(t(1,1:end-2), z_sf)
title ("Experimento 1(b) - Comparação z com filtro e z")
legend ('z com filtro', 'z sem filtro')
xlabel ('Tempo [s]')

index_fig= index_fig + 1;
figure (index_fig)
hold on
grid on
plot(t, phi_b)
title ("Experimento 1(b) -Sinal \phi")
hLeg = legend ("$\dot{x}$", '$x$');
set(hLeg,'Interpreter','latex');
xlabel ('Tempo [s]')
clc
clear 
close all

index_fig = 0;

%% Experiência 2

% Parametros do problema
K1= 1;
K2= 12;
K3= 100;

t= linspace(1,5e7,10e5); % Vetor de tempo
r= 2.7*cos(4e-7*t); % Sinal de entrada

lambda= 1; % Parametro do filtro

G2= tf(1, [K1 K2 K3]); % FT real
y= lsim(G2, r, t); % Simulação original do problema

%% Experiência 2 - Item (a)

s= tf("s");

filtro = tf(1 , [1 2*lambda lambda^2]); % Implementacao do filtro
z= lsim(filtro, r, t); % Sinal z do modelo parametrico

index_fig= index_fig + 1;
figure (index_fig) 
hold on 
grid on 
plot (t, r)
plot (t, z, '--', 'LineWidth', 2)
title ("Experiência 2(a) - Comparação entre r e z")
legend ('r', 'z')
xlabel ('Tempo [s]')
hold off

% Sinais obtidos a partir do filtro
ydd_f= lsim(s^2*filtro, y, t);     
yd_f= lsim(s*filtro, y, t);
y_f= lsim(filtro, y, t); 

phi= [ydd_f yd_f y_f]';  % Vetor phi com filtro

index_fig= index_fig + 1;
figure (index_fig)
hold on
grid on 
plot (t, phi)
title (" Experiência 2(a) - Sinal \phi")
hLeg = legend ('$\ddot{y}$', "$\dot{y}$", '$y$');
set(hLeg,'Interpreter','latex');
xlabel ('tempo (s)')

index_fig= index_fig + 1;
figure (index_fig)
bode(filtro*s^2);                                                                                                                                                   
title("Experiência 2(a) - Resposta do filtro s^2")
grid on 

%% Experiencia 2 - Item (b)
z_b= lsim(filtro, y, t);

ydd_f= lsim(s^2*filtro, y, t);     
yd_f= lsim(s*filtro, y, t);

r_f = lsim(filtro, r, t);

phi_b= [ydd_f yd_f r_f]';  

index_fig= index_fig + 1;
figure (index_fig)
hold on 
grid on
plot(t, z_b) 
title ("Experiência 2(b) - Sinal z")
legend ('z')
xlabel ('tempo (s)')

index_fig= index_fig + 1;
figure (index_fig)
hold on
grid on
plot(t, phi_b)
title ("Experiência 2(b) - Sinal \phi")
hLeg = legend ('$\ddot{y}$', "$\dot{y}$", '$r$');
set(hLeg,'Interpreter','latex');
xlabel ('Tempo [s]')

%% Experiencia 2 - Item (c)
z_c= z - K2*yd_f;  % Sinal z com K2 conhecido
phi_c= [ydd_f y_f]; % Sinal phi com K2 conhecido

index_fig= index_fig + 1;
figure (index_fig)
hold on 
grid on
plot(t, z_c) 
title ("Experiência 2(c) - Sinal z")
legend ('z')
xlabel ('Tempo [s]')

index_fig= index_fig + 1;
figure (index_fig)
hold on
grid on
plot(t, phi_c)
title ("Experiência 2(c) - Sinal \phi")
hLeg = legend ("$\ddot{y}$", '$y$');
set(hLeg,'Interpreter','latex');
xlabel ('Tempo [s]')
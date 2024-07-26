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

lambda= 1; % Parametro do filtro

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

%% Experiencia 4

% Sistema da questão 1

% Parametros do problema
M = 100;
f= 0.15;
k= 7;

t = linspace(0,25,1000); % Vetor de tempo
u = 1+cos(pi*(t/3)); % Sinal de entrada
z= lsim(filtro, u, t); % Sinal z do modelo parametrico

lambda= 1; % Parametro do filtro

G= tf(1, [M f k]); % Resposta original da planta
x= lsim(G, u, t);

s= tf("s");
filtro = tf(1 , [1 2*lambda lambda^2]); % Declaracao do filtro

% Aplicacoes de filtro para o vetor phi
xdd_f= lsim(s^2*filtro, x,t);     
xd_f= lsim(s*filtro, x, t);
x_f= lsim(filtro, x, t); 

%Vetor phi com filtro 
phi= [xdd_f xd_f x_f]'; 

%% Experiência 4 - MMQ

% Matriz inicial P
P0 = 10000*[200 0 0;
      0 1 0;
      0 0 100];
  
% Matriz inicial theta
theta_0 = [1;
           0;
           1;];
       
% Passo de simulacao
dt = t(2)-t(1);

% Progressão de theta
theta_v = [theta_0];
P_m= [];
% Algoritmo MMQ
for i=1:length(t)
    phi_i = phi(:,i); % Valor de phi pro instante i
    
    ms = 1+phi_i'*phi_i; % Valor de ms para o instante i
    
    P_d = -P0*phi_i*phi_i'*P0/(ms^2); % Calculo da derivada de P
    P = P0 + P_d*dt; % Atualização de P
    
    z_i = z(i);  % Valor de z para o instante i
    z_hat = theta_0'*phi_i; % Estimativa de z com os parametros atuais
    
    e = (z_i - z_hat)/ms^2; % Calculo de erro
    
    theta_d = P*e*phi_i;  % Calculo da derivada de theta
    theta = theta_0 + theta_d*dt; % Atualização de theta
    
    theta_v(:, end+1) = theta; % Armazenamento da progressao de theta_v
    P_sin= svds(P);
    P_m(:, end+1)= P_sin;
    
    % Valores para a proxima iteração
    P0 = P;
    theta_0 = theta;
end

%% Experimento 4 - Plots

% Plot M

linha_referencia = 100*ones(1,length(t)); % Valor de referência para M

index_fig = index_fig + 1;
figure(index_fig)
grid on
hold on
plot(t, theta_v(1,1:length(t)))
plot(t, linha_referencia, '--')
ylim([0,110])
title("Experiencia 4 - Parametro M");
hLeg = legend ("$\hat{M}$", 'M');
set(hLeg,'Interpreter','latex');
xlabel ('Tempo [s]');
ylabel("Parametro");

%Plot f

linha_referencia = 0.15*ones(1,length(t)); % Valor de referência para f

index_fig = index_fig + 1;
figure(index_fig)
grid on
hold on
plot(t, theta_v(2,1:length(t)))
plot(t, linha_referencia, '--')
title("Experiencia 4 - Parametro f");
hLeg = legend ("$\hat{f}$", 'f');
set(hLeg,'Interpreter','latex');
xlabel ('Tempo [s]');
ylabel("Parametro");

%Plot k

linha_referencia = 7*ones(1,length(t)); % Valor de referência para k

index_fig = index_fig + 1;
figure(index_fig)
grid on
hold on
plot(t, theta_v(3,1:length(t)))
plot(t, linha_referencia, '--')
title("Experiencia 4 - Parametro k");
hLeg = legend ("$\hat{k}$", 'k');
set(hLeg,'Interpreter','latex');
xlabel ('Tempo [s]');
ylabel("Parametro");

%% Experimento 5 - modelo parametrico - não funcionou

Ra= 1.36;
La= 0.0036;
Km= 0.838;
Kb=Km;
J= 0.001;
b= 0.268;
t= linspace(0,10,10000);
Va= sin(pi*t/3)+sin(pi*t/2)+sin(5*pi*t/3)+sin(100*pi*t/3)+sin(0.1*pi*t/3)+sin(400*pi*t/3)+sin(38*pi*t/3)+sin(0.001*pi*t/3)+sin(42*pi*t/3);
lambda= 1;
s= tf("s");

K1= Km/(J*La);
K2= (J*Ra + b*La)/(J*La);
K3= (Kb*Km)/(J*La); 

G_5= tf(K1, [1 K2 K3]);
w= lsim(G_5, Va, t);

filtro= 1/(s+lambda)^2;

z= lsim(filtro*s^2,w,t);
phi_1= lsim(filtro,Va,t);
phi_2= lsim(-filtro*s,w,t);
phi_3= lsim(-filtro,w,t);

phi= [phi_1 phi_2 phi_3]';


%% Experimento 5 - modelo parametrico - Tentativa 2

% Parâmetros do problema
Ra= 1.36;
La= 0.0036;
Km= 0.838;
Kb=Km;
J= 0.001;
b= 0.268;
t= linspace(0,10,10000);

lambda= 1;

% Constantes simplificadas do problema
K1= Km/(J*La);
K2= (J*Ra + b*La)/(J*La);
K3= (Kb*Km)/(J*La); 

% Vetor de tempo
t = linspace(0,25,1000); 
% Sinal de entrada (Nesse caso, um sinal rico)
Va= sin(pi*t/3)+sin(pi*t/2)+sin(5*pi*t/3);%+sin(100*pi*t/3)+sin(0.1*pi*t/3)+sin(400*pi*t/3)+sin(38*pi*t/3)+sin(0.001*pi*t/3)+sin(42*pi*t/3);

% Função de transferência
G_5= tf(K1, [1 K2 K3]); 
% Obtenção do sinal de entrada w
w= lsim(G_5, Va, t); 

% Declaração do filtro
s= tf("s");
filtro= 1/(s+lambda)^2;

% Sinais do modelo paramétrico
z= lsim(filtro,Va,t);

phi_1= lsim(filtro*s^2,w,t);
phi_2= lsim(filtro*s,w,t);
phi_3= lsim(filtro,w,t);

phi= [phi_1 phi_2 phi_3]';

%% Experimento 5 - MMQ

% Matriz inicial P
P0= 10*[1000 0 0; 0 10 0; 0 0 10];

% Matriz inicial theta
theta_0= [0;0;1];

dt = t(2)-t(1);
theta_m= [];

for i=1:length(t)
    phi_i = phi(:,i); % Valor de phi pro instante i
    
    ms = 1+phi_i'*phi_i; % Valor de ms para o instante i
    
    P_d = -P0*phi_i*phi_i'*P0/(ms^2); % Calculo da derivada de P
    P = P0 + P_d*dt; % Atualização de P
    
    z_i = z(i);  % Valor de z para o instante i
    z_hat = theta_0'*phi_i; % Estimativa de z com os parametros atuais
    
    e = (z_i - z_hat)/(ms^2); % Calculo de erro
    
    theta_d = P*e*phi_i;  % Calculo da derivada de theta
    theta = theta_0 + theta_d*dt; % Atualização de theta
    
    theta_m(:, end+1) = theta; % Armazenamento da progressao de theta_v

    % Valores para a proxima iteração
    P0 = P;
    theta_0 = theta;
end

%% Experimento 5 - plots

one= ones(1,length(t));
linha_referencia1 = one*(1/K1);
linha_referencia2 = one*(K2/K1);
linha_referencia3 = one*(K3/K1);

index_fig = index_fig + 1;
figure(index_fig)
grid on
hold on
plot(t, theta_m')
plot(t, linha_referencia1, '--')
plot(t, linha_referencia2, '--')
plot(t, linha_referencia3, '--')
title("Experiencia 5 - Parametro M");
hLeg = legend ("1/K1", "K2/K1", "K3/K1","1/K1", "K2/K1", "K3/K1");
set(hLeg,'Interpreter','latex');
xlabel ('Tempo [s]');
ylabel("Parametro");

%% Experiência 7 - (a)

clear;
clc;
close all;
index_fig = 0;

% Parâmetros da planta
a = -1;
b = 5;
% Parâmetros do modelo de referência
am = 10;
bm = 1;

% Valores ideais de controle
k_star = (a+am)/b;
l_star = bm/b;

% Vetor de tempo
dt = 0.01;
t = 0:dt:10;

% Sinal de entrada
r = 5*cos(t) + 2*sin(pi*t);

% Valores de x para o controlador ideal
x0 = 0;
x = zeros(1,length(t));
x(1) = x0;

% Função de transferência modelo de referência
s = tf('s');
G = bm/(s+am);
xm = lsim(G, r, t);

% Valores iniciais de adaptação
k0 = 1;
l0 = 1;
gamma1 = 50;
gamma2 = 0.5;

% Valores de x para o controlador adaptativo
x_hat = zeros(1,length(t));
x0_hat = 0;
x_hat(1) = x0_hat;

% Vetores para progressão dos ganhos e do erro
k = zeros(1,length(t));
k(1) = k0;
l = zeros(1,length(t));
l(1) = l0;
e = zeros(1,length(t));
e(1) = x_hat(1)-xm(1);

for i=2:length(t)
    % X com controlador adaptativo
    e(i-1) = x_hat(i-1) - xm(i-1);
    
    k_dot = gamma1*e(i-1)*x_hat(i-1)*sign(b);
    k(i) = k_dot*dt + k(i-1);
    
    l_dot = -gamma2*e(i-1)*r(i-1)*sign(b);
    l(i) = l_dot*dt + l(i-1);
    
    u_hat = -k(i)*x_hat(i-1) + l(i)*r(i-1);
    x_hat_dot = a*x_hat(i-1) + b*u_hat;
    x_hat(i) = x_hat_dot * dt + x_hat(i-1);
    
    % X como controlador ideal
    u = -k_star*x(i-1) + l_star*r(i-1);
    x_dot = a*x(i-1) + b*u;
    x_i = x(i-1) + x_dot*dt;
    x(i) = x_i;
end

%% Experiência 7 (a) - Plots

figure(1)
hold on
grid on
plot(t, x, "LineWidth", 1)
plot(t, xm, '--', "LineWidth", 2)
plot(t, x_hat, "LineWidth", 1)
legend("Controlador Ideal", "Modelo de Referência", "Adaptativo")
xlabel("Tempo [s]")
ylabel("Amplitude")
title("Experiência 7(a) - Estado x")


figure(2)
hold on
grid on
plot(t,e)
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(a) - Erro de adaptação")

k_star_vector = k_star*ones(1, length(t));
figure(3)
hold on
grid on
plot(t,k)
plot(t, k_star_vector, '--')
legend("Adaptativo", "Ideal")
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(a) - Parâmetro de controle k")
ylim([0, 9])

l_star_vector = l_star*ones(1, length(t));
figure(4)
hold on
grid on
plot(t,l)
plot(t, l_star_vector, '--')
legend("Adaptativo", "Ideal")
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(a) - Parâmetro de controle l")
ylim([0, 1.5])

%% Experiência 7 - (b)

close all;
clc;
clear;
index_fig = 0;

a = -1;
b = 12;
am = 2;
bm = 2;

% Valores ideais de controle
k_star = (a+am)/b;
l_star = bm/b;

% Intervalo de tempo
dt = 0.01;
t = 0:dt:10;

% Valores de referência
r1 = 5*ones(1,length(t));
r2 = sin(2*t);
r3 = 1./(1+t);

r = [r1; r2; r3];

% Função de transferência modelo de referência
s = tf('s');
G = bm/(s+am);

xm_1 = lsim(G, r1, t);
xm_2 = lsim(G, r2, t);
xm_3 = lsim(G, r3, t);

xm = [xm_1 xm_2 xm_3]';

% Valores iniciais de adaptação
k0 = 1;
l0 = 1;
gamma1 = 1;
gamma2 = 20;

% Valores de x para o controlador adaptativo
x_hat = zeros(3,length(t));
x0_hat = 0;
x_hat(1) = x0_hat;

% Vetores para progressão dos ganhos e do erro
k = zeros(3,length(t));
k(1) = k0;
l = zeros(3,length(t));
l(1) = l0;
e = zeros(3,length(t));
e(1) = x_hat(1)-xm(1);

for j=[1,2,3]
    for i=2:length(t)
        % X com controlador adaptativo
        e(j,i-1) = x_hat(j,i-1) - xm(j,i-1);

        k_dot = gamma1*e(j,i-1)*x_hat(j,i-1)*sign(b);
        k(j,i) = k_dot*dt + k(j,i-1);

        l_dot = -gamma2*e(j,i-1)*r(j,i-1)*sign(b);
        l(j,i) = l_dot*dt + l(j,i-1);

        u_hat = -k(j,i)*x_hat(j,i-1) + l(j,i)*r(j,i-1);
        x_hat_dot = a*x_hat(j,i-1) + b*u_hat;
        x_hat(j,i) = x_hat_dot * dt + x_hat(j,i-1);
    end
end

%% Experiência 7 (b) - Plots

% Plot de x

index_fig = index_fig + 1;
figure(index_fig)
subplot(3,1,1)
hold on
grid on
plot(t, xm(1,:), '--', "LineWidth", 2)
plot(t, x_hat(1,:), "LineWidth", 1)
legend("Modelo de Referência", "Adaptativo")
xlabel("Tempo [s]")
ylabel("Amplitude")
title("Experiência 7(b): r=5")
subplot(3,1,2)
hold on
grid on
plot(t, xm(2,:), '--', "LineWidth", 2)
plot(t, x_hat(2,:), "LineWidth", 1)
legend("Modelo de Referência", "Adaptativo")
xlabel("Tempo [s]")
ylabel("Amplitude")
title("Experiência 7(b) - sin(2t)")
subplot(3,1,3)
hold on
grid on
plot(t, xm(3,:), '--', "LineWidth", 2)
plot(t, x_hat(3,:), "LineWidth", 1)
legend("Modelo de Referência", "Adaptativo")
xlabel("Tempo [s]")
ylabel("Amplitude")
title("Experiência 7(b)  r=1/(1+t)")

% Plot de e

index_fig = index_fig + 1;
figure(index_fig)
subplot(3,1,1)
hold on
grid on
plot(t,e(1,:))
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(b): r=5")
subplot(3,1,2)
hold on
grid on
plot(t,e(2,:))
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(b) - sin(2t)")
subplot(3,1,3)
hold on
grid on
plot(t,e(3,:))
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(b)  r=1/(1+t)")

k_star_vector = k_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,k)
plot(t, k_star_vector, '--')
legend("r=5", "r=sin(2t)", "r=1/(1+t)", "Ideal")
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(b) - Parâmetro de controle k")

l_star_vector = l_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,l)
plot(t, l_star_vector, '--')
legend("r=5", "r=sin(2t)", "r=1/(1+t)", "Ideal")
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(b) - Parâmetro de controle L")

%% Experiência 7 - (c)

close all;
clc;
clear;

%%
index_fig = 0;

a = 0.5;
b = 1.5;
d = 10;

am = 0.5;
bm = 0.5;

% Valores ideais de controle
k_star = (am-a)/b;
l_star = bm/b;
delta_star = -d*a/b;

% Intervalo de tempo
dt = 0.01;
t = 0:dt:10;

% Valores de referência
Vs = 55*ones(1,length(t));

% Função de transferência modelo de referência
s = tf('s');
G = bm/(s+am);

Vm = lsim(G, Vs, t);

% Valores iniciais de adaptação
k0 = 1;
l0 = 1;
delta0 = 1;

gamma1 = 2;
gamma2 = 1;
gamma3 = 1;

% Valores de V para o controlador adaptativo
V_hat = zeros(1,length(t));
V0_hat = 0;
V_hat(1) = V0_hat;

% Vetores para progressão dos ganhos e do erro
k = zeros(1,length(t));
k(1) = k0;
l = zeros(1,length(t));
l(1) = l0;
delta = zeros(1,length(t));
delta(1) = delta0;

e = zeros(1,length(t));
e(1) = V_hat(1)-Vm(1);

for i=2:length(t)
    % X com controlador adaptativo
    e(i-1) = V_hat(i-1) - Vm(i-1);

    k_dot = gamma1*e(i-1)*V_hat(i-1)*sign(b);
    k(i) = k_dot*dt + k(i-1);

    l_dot = -gamma2*e(i-1)*Vs(i-1)*sign(b);
    l(i) = l_dot*dt + l(i-1);
    
    delta_dot = -gamma3*e(i-1)*sign(b);
    delta(i) = delta_dot*dt + delta(i-1);

    u_hat = -k(i)*V_hat(i-1) + l(i)*Vs(i-1) + delta(i);
    V_hat_dot = - a*V_hat(i-1) + b*u_hat + d*a;
    V_hat(i) = V_hat_dot * dt + V_hat(i-1);
end

%% Plots

% Plot de x

index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t, Vm, '--', "LineWidth", 1)
plot(t, V_hat, "LineWidth", 1)
legend("Modelo de Referência", "Adaptativo")
xlabel("Tempo [s]")
ylabel("Velocidade")
title("Experiência 7(c) - Velocidade")

% Plot de e

index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,e)
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(c) - erro")

k_star_vector = k_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,k)
plot(t, k_star_vector, '--')
legend("Estimado", "Ideal")
xlabel("Tempo [s]")
title("Experiência 7(b) - Parâmetro de controle k")

l_star_vector = l_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,l)
plot(t, l_star_vector, '--')
legend("Estimado", "Ideal")
xlabel("Tempo [s]")
title("Experiência 7(c) - Parâmetro de controle L")

delta_star_vector = delta_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,l)
plot(t, delta_star_vector, '--')
legend("Estimado", "Ideal")
xlabel("Tempo [s]")
title("Experiência 7(c) - Parâmetro de controle delta")
clc
clear 
close all

index_fig = 0;

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
P0= 10*[100 0 0; 0 10 0; 0 0 10];

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

%% Experimento 5 - modelo parametrico - não funcionou

% Ra= 1.36;
% La= 0.0036;
% Km= 0.838;
% Kb=Km;
% J= 0.001;
% b= 0.268;
% t= linspace(0,10,10000);
% Va= sin(pi*t/3)+sin(pi*t/2)+sin(5*pi*t/3)+sin(100*pi*t/3)+sin(0.1*pi*t/3)+sin(400*pi*t/3)+sin(38*pi*t/3)+sin(0.001*pi*t/3)+sin(42*pi*t/3);
% lambda= 1;
% s= tf("s");
% 
% K1= Km/(J*La);
% K2= (J*Ra + b*La)/(J*La);
% K3= (Kb*Km)/(J*La); 
% 
% G_5= tf(K1, [1 K2 K3]);
% w= lsim(G_5, Va, t);
% 
% filtro= 1/(s+lambda)^2;
% 
% z= lsim(filtro*s^2,w,t);
% phi_1= lsim(filtro,Va,t);
% phi_2= lsim(-filtro*s,w,t);
% phi_3= lsim(-filtro,w,t);
% 
% phi= [phi_1 phi_2 phi_3]';
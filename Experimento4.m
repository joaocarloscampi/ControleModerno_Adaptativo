clc
clear 
close all

index_fig = 0;

%% Experiencia 4

% Sistema da questão 1

% Parametros do problema
M = 100;
f= 0.15;
k= 7;

lambda= 1; % Parametro do filtro

t = linspace(0,25,1000); % Vetor de tempo
u = 1+cos(pi*(t/3)); % Sinal de entrada

filtro = tf(1 , [1 2*lambda lambda^2]); % Implementacao do filtro
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